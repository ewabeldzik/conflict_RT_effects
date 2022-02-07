#!/usr/bin/env tcsh
set cyan='\e[96m'
set white='\e[0m'

## Scripts to regression analysis with nuisanse regressors only
## -------------------------------------
## set project specific parameters below
## -------------------------------------

# edit directory do dataset 
set datadir = "/media/sf_mri/data"
set outdir = "/media/sf_mri/ica/data"
set maskfile = "/media/sf_mri/maps/maskMNI_cortex_3.5mm.nii"

# edit subject/session numbers
set subjects = (`count -digits 2 1 37`)

# edit name of functional tasks
set task_list = ('simon2' 'simon5' 'simple' 'stroop2' 'stroop5')


cd ${datadir}

# start subject loop
foreach sb ( ${subjects} )
echo "${cyan}\n >>> regression for sb$sb starts <<< ${white}"

# start session loop
foreach task ( ${task_list} )

	cd sb$sb/$task/
	if (! -f sb$sb.$task.Rwherr.nii ) then 

	set nb_vol = `3dinfo -nt sb$sb.$task.08.res.nii`
	set reptime = `3dinfo -tr sb$sb.$task.08.res.nii`

	echo "${cyan}---- creating matrix for sb$sb.$task nb_vol=${nb_vol} TR=${reptime} ${white}"
	3dDeconvolve -polort A -num_stimts 13 -jobs 6 \
	-nodata ${nb_vol} ${reptime} \
	-stim_file 1 sb$sb.$task.00.motion.demean.1D'[1]' -stim_base 1 -stim_label 1 roll_01 \
	-stim_file 2 sb$sb.$task.00.motion.demean.1D'[2]' -stim_base 2 -stim_label 2 pitch_01 \
	-stim_file 3 sb$sb.$task.00.motion.demean.1D'[3]' -stim_base 3 -stim_label 3 yaw_01 \
	-stim_file 4 sb$sb.$task.00.motion.demean.1D'[4]' -stim_base 4 -stim_label 4 dS_01 \
	-stim_file 5 sb$sb.$task.00.motion.demean.1D'[5]' -stim_base 5 -stim_label 5 dL_01 \
	-stim_file 6 sb$sb.$task.00.motion.demean.1D'[6]' -stim_base 6 -stim_label 6 dP_01 \
	-stim_file 7 sb$sb.$task.00.motion.deriv.1D'[1]' -stim_base 7 -stim_label 7 roll_02 \
	-stim_file 8 sb$sb.$task.00.motion.deriv.1D'[2]' -stim_base 8 -stim_label 8 pitch_02 \
	-stim_file 9 sb$sb.$task.00.motion.deriv.1D'[3]' -stim_base 9 -stim_label 9 yaw_02 \
	-stim_file 10 sb$sb.$task.00.motion.deriv.1D'[4]' -stim_base 10 -stim_label 10 dS_02 \
	-stim_file 11 sb$sb.$task.00.motion.deriv.1D'[5]' -stim_base 11 -stim_label 11 dL_02 \
	-stim_file 12 sb$sb.$task.00.motion.deriv.1D'[6]' -stim_base 12 -stim_label 12 dP_02 \
	-stim_file 13 sb$sb.$task.CSFsignal.1D -stim_base 13 -stim_label 13 csf    \
	-x1D sb$sb.$task.xmat.denoise.1D -xjpeg sb$sb.$task.xmat.denoise.jpg

	mkdir -p $outdir/sb$sb/$task/

	echo "${cyan}---- running regression with 3dREMLfit for sb$sb.$task ${white}"
	3dREMLfit -matrix sb$sb.$task.xmat.denoise.1D -input sb$sb.$task.08.res.nii \
		-mask $maskfile -Rwherr $outdir/sb$sb/$task/sb$sb.$task.Rwherr.nii 
	

	endif # if file does not exist
	cd ../../

end # end task loop

echo "${cyan}\n >>> regression for sb$sb ends <<< ${white}"
end # end subject loop
