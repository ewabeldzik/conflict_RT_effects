#!/bin/bash

## ver 1.02.
## scripts transforms dicom files to niifti and created bids structure toghether with defacing of structural scans
## to run it just EDIT project-dependent values and write in terminal "./our_dicom2bids.sh"
## it assumes that in dicomdir folder you have separate folders for each sequence in which there are folders for sbjects etc.

# EDIT directory to dicom files ordered in folder structure: project > seqname > sbject > session 
dicomdir="/media/sf_mri"
#mkdir ${dicomdir}/bids
outdir="${dicomdir}/bids"
#mkdir ${dicomdir}/anat_wholehead_nifti
wholehead_outdir="${dicomdir}/anat_wholehead_nifti"
# EDIT directory to mri_deface templates
brain_template=/home/neuron/abin/deface/talairach_mixed_with_skull.gca
face_template=/home/neuron/abin/deface/face.gca

cd ${dicomdir}/

# sbject loop
for sb in {01..37}; do
# sequence loop
for seqname in 'anat' 'fieldmap_m' 'fieldmap_p' 'simon2' 'simon5' 'simple' 'stroop2' 'stroop5' ; do
# 
	## checking if data exist
    	if [ -d raw/sb$sb/$seqname/ ]; then
		echo " >> Folder with raw data for sub-${sb}  ${seqname} exist! <<"

	## dcm2niix for each sequence
	if [ ${seqname} = 'anat' ]; then
		echo "processing anat for sub-${sb} "
	
	mkdir -p ${outdir}/sub-${sb}/anat/ 
	mkdir -p ${wholehead_outdir}/sub-${sb}/ 
	cd raw/sb${sb}/${seqname}/

	dcm2niix -z y -f sub-${sb}_T1w -o ${wholehead_outdir}/sub-${sb}/ .
	mri_deface ${wholehead_outdir}/sub-${sb}/sub-${sb}_T1w.nii.gz ${brain_template} ${face_template} ${outdir}/sub-${sb}/anat/sub-${sb}_T1w.nii.gz
	mv ${wholehead_outdir}/sub-${sb}/sub-${sb}_T1w.json ${outdir}/sub-${sb}/anat/

		if [ -f ${outdir}/sub-${sb}/anat/sub-${sb}_T1w.nii.gz ]; then echo " >> Data for sub-${sb}  anat exist! <<"
		else echo "Missing anat for sub-${sb} " >> ${outdir}/missing_data.log
		fi

	cd ../../../
	
	# EDIT functional sequence name
	elif [ ${seqname} = 'simon2' ] || [ ${seqname} = 'simon5' ] || [ ${seqname} = 'simple' ] || [ ${seqname} = 'stroop5' ] || [ ${seqname} = 'stroop2' ] ; then


	mkdir -p ${outdir}/sub-${sb}/func
	cd raw/sb${sb}/${seqname}/
		echo "processing func for sub-${sb} "

	dcm2niix -z y -f sub-${sb}_task-${seqname}_bold -o ${outdir}/sub-${sb}/func/ .

		if [ -f ${outdir}/sub-${sb}/func/sub-${sb}_task-${seqname}_bold.nii.gz ]; then echo " >> Data for sub-${sb}  func exist! <<"
		else echo "Missing func for sub-${sb} " >> ${outdir}/missing_data.log
		fi

	cd ../../../
	elif [ ${seqname} = 'fieldmap_m' ]; then
	mkdir -p ${outdir}/sub-${sb}/fmap
	cd raw/sb${sb}/${seqname}/
		echo "processing fmap_m for sub-${sb} "

	dcm2niix -z y -f sub-${sb}_magnitude -o ${outdir}/sub-${sb}/fmap/ .

		if [ -f ${outdir}/sub-${sb}/fmap/sub-${sb}_magnitude_e1.nii.gz ]; then echo " >> Data for sub-${sb}  fmap_mag_e1 exist! <<"
		elif  [ -f ${outdir}/sub-${sb}/fmap/sub-${sb}_magnitude.nii.gz ]; then echo " >> Data for sub-${sb}  fmap_mag exist! <<"
		else echo "Missing fmap_mag_e1 for sub-${sb} " >> ${outdir}/missing_data.log
		fi
		if [ -f ${outdir}/sub-${sb}/fmap/sub-${sb}_magnitude_e2.nii.gz ]; then echo " >> Data for sub-${sb}  fmap_mag_e2 exist! <<"
		else echo "Missing fmap_mag_e2 for sub-${sb} " >> ${outdir}/missing_data.log
		fi

	cd ../../../
	elif [ ${seqname} = 'fieldmap_p' ]; then
	mkdir -p ${outdir}/sub-${sb}/fmap
	cd raw/sb${sb}/${seqname}/
		echo "processing fmap_p for sub-${sb} "

	dcm2niix -z y -f sub-${sb}_phasediff -o ${outdir}/sub-${sb}/fmap/ .

		if [ -f ${outdir}/sub-${sb}/fmap/sub-${sb}_phasediff_e2_ph.nii.gz ]; then echo " >> Data for sub-${sb}  map_phase exist! <<"
		else echo "Missing fmap_phase for sub-${sb} " >> ${outdir}/missing_data.log
		fi

	cd ../../../
	else 
		echo ">>> SOMETHING WITH THIS LOOP IS WRONG!!! <<<"
	fi

	else 
		echo "Missing ${seqname} for sub-${sb} " >> ${outdir}/missing_folders_with_rawdata.log
	fi
done
done

exit

