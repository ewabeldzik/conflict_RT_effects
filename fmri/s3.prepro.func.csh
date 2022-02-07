#!/usr/bin/env tcsh
#!/usr/bin/env python2.7
set cyan='\e[96m'
set white='\e[0m'


## Scritp to: reorient/despike/time shifting
## -------------------------------------
## Set project specific parameters below
## Prepare slicetiming.1D
## -------------------------------------

# EDIT directory do dataset 
set bidsdir = "/media/sf_mri/bids"
set outdir = "/media/sf_mri/data"

# EDIT subject/session numbers
set subjects = (`count -digits 2 1 37`)

# EDIT name of functional tasks
set task_list = ('simon2' 'simon5' 'simple' 'stroop2' 'stroop5')



cd ${outdir}

# start subject loop
foreach sb ( ${subjects} )
echo "${cyan}\n >>> prepro part I for sb$sb starts <<< ${white}"


# creating folder to check the aligment
if (! -d ${bidsdir}/check/ ) then 
mkdir -p sb$sb/check/
cp sb$sb/anat/anatSS.sb$sb.nii sb$sb/check/
endif

# start session loop
foreach task ( ${task_list} )

	if (! -d sb$sb/$task/ ) then 
	mkdir -p sb$sb/$task/
	endif

	cd sb$sb/$task/

	if (! -f sb$sb.$task.00.raw.nii ) then 
	echo "${cyan}---- copy data from BIDS to DATA directory for sb$sb.$task ${white}"
	3dcopy ${bidsdir}/sub-$sb/func/sub-${sb}_task-${task}_bold.nii.gz \
		${outdir}/sb$sb/$task/sb$sb.$task.00.raw.nii
	endif


	if (! -f sb$sb.$task.00.ro.nii ) then 
	echo "${cyan}---- reorienting for sb$sb $task ${white}"
	3dresample -orient rpi -prefix sb$sb.$task.00.ro.nii -input sb$sb.$task.00.raw.nii
	endif


	if (! -f  sb$sb.$task.00.motion.1D) then
	echo "${cyan}---- calculating motion estimates for sb$sb $task ${white}"
	3dvolreg -Fourier -base 0 -dfile sb$sb.$task.00.motion.1D \
		-1Dmatrix_save sb$sb.$task.00.motion -prefix NULL sb$sb.$task.00.ro.nii

	# compute de-meaned and derivatives of motion parameters (for use in regression)
	1d_tool.py -infile sb$sb.$task.00.motion.1D -set_nruns 1 -demean \
		-write sb$sb.$task.00.motion.demean.1D
	1d_tool.py -infile sb$sb.$task.00.motion.1D -set_nruns 1 -derivative -demean \
		-write sb$sb.$task.00.motion.deriv.1D 
	endif
	

	if (! -f sb$sb.$task.01.ds.nii) then
		echo "${cyan}---- despiking for sb$sb $task ${white}"
		3dDespike -NEW -nomask -localedit -ssave sb$sb.$task.01.ds.spikes.nii \
				-prefix sb$sb.$task.01.ds.nii sb$sb.$task.00.ro.nii
	endif
		


	if (! -f sb$sb.$task.02.ts.nii) then
		echo "${cyan}---- slice timing correction for sb$sb $task ${white}"
		3dTshift -tpattern @${outdir}/slicetiming.1D -prefix sb$sb.$task.02.ts.nii sb$sb.$task.01.ds.nii
	endif



	if (! -f sb$sb.$task.03.zp.nii) then
	echo "${cyan}---- calculating matrix for coregistration to anat for sb$sb $task ${white}" 
	3dbucket -fbuc -prefix sb$sb.$task.02.vol.nii sb$sb.$task.02.ts.nii'[0]'
	align_epi_anat.py -giant_move -ex_mode quiet -deoblique on -anat_has_skull no -cost lpc  \
		-epi_base 0 -epi2anat -tshift off -volreg off -suffix .al -anat ../anat/anatSS.sb$sb.nii \
		-epi sb$sb.$task.02.vol.nii
	3dAFNItoNIFTI -prefix sb$sb.$task.02al.nii sb$sb.$task.02.vol.al+orig 
	rm sb$sb.$task.02.vol.al+orig* sb$sb.$task.02.vol.nii
	3dZeropad -master sb$sb.$task.02al.nii \
		-prefix sb$sb.$task.03.zp.nii sb$sb.$task.02.ts.nii
	3dcopy sb$sb.$task.02al.nii ../check/
	endif
	

	cd ../../

end # end task loop
echo "${cyan}>>> prepro part I for sb$sb ends <<< ${white}"


end # end subject loop
echo "${cyan}\n ---- script done,check aligment visually ----- \n ${white}"


exit



## Versions:
#18.0.1. > changed mask creation process and 3dBlurInMask instead of 3dmerge
