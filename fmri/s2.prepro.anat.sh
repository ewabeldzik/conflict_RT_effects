#!/usr/bin/env tcsh
set cyan='\e[96m'
set white='\e[0m'


## scripts deals with structrural data, i.e. prepare it for fMRI analysis
## it assumes that input directory (bidsdir) has BIDS structure

# EDIT directory do dataset 
set bidsdir = "/media/sf_mri/bids"
set outdir = "/media/sf_mri/data"
set atlaspath = "/usr/share/fsl/5.0/data/atlases/HarvardOxford/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz"

cd ${outdir}

if ( -f ventricles.nii ) then 
else
	echo "${cyan} >>> create map of ventricles ${white}"
	3dcalc -a $atlaspath -expr 'equals(a,3)+equals(a,14)' -prefix ventricles.nii
endif

## # EDIT subject numbers 
set subjects = (`count -digits 2 32 37`)

foreach sb ( $subjects )

mkdir -p sb$sb/anat

echo "${cyan} >> copy data from BIDS to DATA directory for sb$sb ${white}"
3dcopy ${bidsdir}/sub-${sb}/anat/sub-${sb}_T1w.nii.gz ${outdir}/sb$sb/anat/sb$sb.T1w.nii.gz

	cd sb$sb/anat

	### Cutting anat
	3dZeropad -I -75 -prefix sb$sb.T1w_cut.nii sb$sb.T1w.nii.gz

	echo "${cyan} >> reorienting for sb$sb ${white}"
	3dresample -orient rpi -prefix sb$sb.T1w_RO.nii -input sb$sb.T1w_cut.nii


	echo "${cyan} >> skull stripping and warping to MNI for sb$sb ${white}"
	@SSwarper -input sb$sb.T1w_RO.nii -base MNI152_2009_template_SSW.nii.gz -subid sb$sb 
	3dAutomask -dilate 1  -prefix anatSS.sb$sb.automask.nii anatSS.sb$sb.nii

	echo "${cyan} >> creating ventricles mask for CSF signal extraction for sb$sb ${white}"
	## inverse warping of ventricles masks from HarvardOxford-sub-maxprob-thr25-1mm.nii.gz
	3dNwarpCat -prefix anatQQ_total_WARP.nii -warp1 anatQQ.sb{$sb}_WARP.nii -warp2 anatQQ.sb$sb.aff12.1D 
	3dNwarpApply -nwarp 'INV(anatQQ_total_WARP.nii)' -source ../../ventricles.nii -master anatSS.sb$sb.nii -prefix sb$sb.ventricles.nii

	3dmask_tool -dilate_input -1 -input sb$sb.ventricles.nii -prefix sb$sb.mask_ventricles.nii 

	cd ../../

echo "${cyan} \n ---- Check warping visually for every subject ----- \n ${white}"

end




