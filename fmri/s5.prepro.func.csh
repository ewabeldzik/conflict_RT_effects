#!/usr/bin/env tcsh
#!/usr/bin/env python2.7
set cyan='\e[96m'
set white='\e[0m'



## Scritp to: scaling/CFS signal extraction/normalizing
## -------------------------------------
## set project specific parameters below
## -------------------------------------

# edit directory do dataset 
set outdir = "/media/sf_mri/data"
set maskfile = "/media/sf_mri/maps/maskMNI_cortex_3.5mm.nii"

# edit subject/session numbers
set subjects = (`count -digits 2 1 37`)

# edit name of functional tasks
set task_list = ('simon2' 'simon5' 'simple' 'stroop2' 'stroop5')

cd ${outdir}


# start subject loop
foreach sb ( ${subjects} )
echo "${cyan}\n >>> prepro part III for sb$sb starts <<< ${white}"


# start session loop
foreach task ( ${task_list} )

 
	cd sb$sb/$task/

	if (! -f sb$sb.$task.06.sc.nii ) then 
	echo "${cyan}---- creating mask and calculating % signal change for sb.$sb $task ${white}"
	3dAutomask -clfrac 0.4 -prefix sb$sb.$task.05.wr.mask.nii sb$sb.$task.05.wr.nii
	3dTstat -prefix sb$sb.$task.05.wr.mean.nii sb$sb.$task.05.wr.nii
	3dcalc -fscale -a sb$sb.$task.05.wr.nii -b sb$sb.$task.05.wr.mean.nii\
		 -c sb$sb.$task.05.wr.mask.nii \
		-expr "100*((a/b)-1)*c*step(b-1)" -prefix sb$sb.$task.06.sc.nii

	echo "${cyan}---- extracting tc of csf for sb.$sb $task ${white}"
	3dresample -prefix sb$sb.$task.mask_ventricles.nii -master sb$sb.$task.06.sc.nii\
		 -inset ../anat/sb$sb.mask_ventricles.nii
	3dmaskave -quiet -mask sb$sb.$task.mask_ventricles.nii \
		-mrange 1 1 sb$sb.$task.06.sc.nii > sb$sb.$task.CSFsignal.1D
	endif


	if (! -f sb$sb.$task.08.res.nii ) then 
	echo "${cyan}---- normalising to standard for sb.$sb $task ${white}"
	3dZeropad -S 7 -prefix sb$sb.$task.06.zp.nii sb$sb.$task.06.sc.nii
	3dNwarpApply -nwarp '../anat/anatQQ_total_WARP.nii' -source sb$sb.$task.06.zp.nii \
		-prefix sb$sb.$task.07.mni.nii
	3drefit -view tlrc -space MNI sb$sb.$task.07.mni.nii

	3dresample -master $maskfile -input sb$sb.$task.07.mni.nii \
		 -prefix sb$sb.$task.08.res.nii

	# compute MAD (median absolute deviation) in case you want to remove very bad voxels
	3dTstat -MAD -prefix sb$sb.$task.MAD.nii sb$sb.$task.08.res.nii 

	#copy one volume to check the coregistration to MNI
	3dbucket -fbuc -prefix sb$sb.$task.08.vol.nii sb$sb.$task.08.res.nii'[0]'
	mv sb$sb.$task.08.vol.nii ../../check/			
	endif


	cd ../../

end # end task loop
echo "${cyan}>>> prepro part III for sb$sb ends <<< ${white}"

end # end subject loop
