#!/usr/bin/env tcsh
#!/usr/bin/env python2.7
set cyan='\e[96m'
set white='\e[0m'


## Scritp to: warping (motion, fmap, coregidtation)
## -------------------------------------
## Set project specific parameters below
## Prepare slicetiming.1D
## -------------------------------------

# edit directory do dataset 
set bidsdir = "/media/sf_mri/bids"
set outdir = "/media/sf_mri/data"


# edit parameters for fieldmap unwarping, i.e. time echo difference (magnitude 1 vs magnitude 2) and effective echo spacing
set echodiff = 2.46
set effechospa = 0.000320001

# edit subject/session numbers
set subjects = (`count -digits 2 1 37`)

# edit name of functional tasks
set task_list = ('simon2' 'simon5' 'simple' 'stroop2' 'stroop5')


cd ${outdir}

# start subject loop
foreach sb ( ${subjects} )
echo "${cyan}\n >>> prepro part II for sb$sb starts <<< ${white}"

cd sb$sb/

# creating folder to fmap
if (! -d "fmap" ) then 
mkdir -p fmap
cd fmap
	
echo "${cyan}---- preparing fmap files for sb$sb  ${white}"
## copying fmap files
3dcopy ${bidsdir}/sub-$sb/fmap/sub-${sb}_magnitude_e2.nii.gz sb$sb.mag.00.nii
3dcopy ${bidsdir}/sub-$sb/fmap/sub-${sb}_phasediff_e2_ph.nii.gz sb$sb.phase.00.nii


## aligning fmap files to epi
3dZeropad -master ../simple/sb$sb.simple.02al.nii -prefix sb$sb.mag.01.zp.nii sb$sb.mag.00.nii
3dZeropad -master sb$sb.mag.01.zp.nii	-prefix sb$sb.phase.01.zp.nii sb$sb.phase.00.nii
align_epi_anat.py -deoblique on -giant_move -ex_mode quiet -anat_has_skull no -cost lpa \
	-epi2anat -epi_base 0 -tshift off -volreg off \
	-suffix .al -anat ../anat/anatSS.sb$sb.nii -epi sb$sb.mag.01.zp.nii

## creating 4D dummy file for Nwarp use
3dcalc -a sb$sb.mag.01.zp.nii -expr '(a*0)' -datum float -prefix sb$sb.fmap_dummyX.nii
3dbucket -abuc -prefix sb$sb.fmap_dummyXYZ.nii sb$sb.fmap_dummyX.nii \
	sb$sb.fmap_dummyX.nii sb$sb.fmap_dummyX.nii
3dNwarpApply -nwarp 'sb'$sb'.fmap_dummyXYZ.nii sb'$sb'.mag.01.zp.al_mat.aff12.1D' \
	-source sb$sb.mag.01.zp.nii -prefix sb$sb.mag.02.al.nii
3dNwarpApply -nwarp 'sb'$sb'.fmap_dummyXYZ.nii sb'$sb'.mag.01.zp.al_mat.aff12.1D' \
	-source sb$sb.phase.01.zp.nii -prefix sb$sb.phase.02.al.nii

## preparing aligned fmap files
bet sb$sb.mag.02.al.nii sb$sb.mag.03.ss.nii  -f 0.5 -g 0
fslmaths sb$sb.mag.03.ss.nii -ero sb$sb.mag.04.ero.nii
fsl_prepare_fieldmap SIEMENS sb$sb.phase.02.al.nii \
	sb$sb.mag.04.ero.nii sb$sb.rads.nii ${echodiff}

cd ../
endif



# start session loop
foreach task ( ${task_list} )
	
	cd $task/

	if (! -f sb$sb.$task.05.wr.nii) then
	echo "${cyan}---- running fugue and creating warp file for sb$sb $task ${white}"
	fugue -i sb$sb.$task.03.zp.nii  --dwell=${effechospa} --unwarpdir=y- \
		--saveshift=sb$sb.$task.04.shift.nii \
		--loadfmap=../fmap/sb$sb.rads.nii -u sb$sb.$task.04.uw.nii
	3dcalc -a sb$sb.$task.04.shift.nii -expr '(a*0)' -prefix sb$sb.$task.04.dummy.nii
	3dbucket -abuc -prefix sb$sb.$task.04.shiftXYZ.nii sb$sb.$task.04.dummy.nii \
	sb$sb.$task.04.shift.nii sb$sb.$task.04.dummy.nii

	echo "${cyan}---- all-in-one warping, i.e. align2anat, unwarp and motion correction for sb$sb $task ${white}"
	3dNwarpApply -nwarp 'sb'$sb'.'$task'.04.shiftXYZ.nii \
	 sb'$sb'.'$task'.00.motion.aff12.1D  sb'$sb'.'$task'.02.vol.al_mat.aff12.1D' \
		-source sb$sb.$task.03.zp.nii -prefix sb$sb.$task.05.wr.nii

	cp sb$sb.$task.05.wr.nii ../check

	endif

	cd ../

end # end task loop
echo "${cyan}>>> prepro part I for sb$sb ends <<< ${white}"

cd ../
end # end subject loop
echo "${cyan}\n ---- script done,check warping visually ----- \n ${white}"


exit


## Versions:
#18.0.1. > changed mask creation process and 3dBlurInMask instead of 3dmerge
