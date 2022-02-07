#!/usr/bin/env tcsh


cd /media/sf_mri/ica/res_frontal25


set labels = ('FPNr' 'ACC' 'DMN' 'FPNl' 'preSMA' 'Motor' 'Aud' 'Ins')
set ic_list = ('13'  '17'  '15'  '20'    '24'    '14'    '21'   '4') #(`count -digits 1 0 35`)
set l = 1
foreach i ( $ic_list)

3dttest++ -setA\
sef_sub001_component_ica_s1_.nii'['$i']'\
sef_sub002_component_ica_s1_.nii'['$i']'\
sef_sub003_component_ica_s1_.nii'['$i']'\
sef_sub004_component_ica_s1_.nii'['$i']'\
sef_sub005_component_ica_s1_.nii'['$i']'\
sef_sub006_component_ica_s1_.nii'['$i']'\
sef_sub007_component_ica_s1_.nii'['$i']'\
sef_sub008_component_ica_s1_.nii'['$i']'\
sef_sub009_component_ica_s1_.nii'['$i']'\
sef_sub010_component_ica_s1_.nii'['$i']'\
sef_sub011_component_ica_s1_.nii'['$i']'\
sef_sub012_component_ica_s1_.nii'['$i']'\
sef_sub013_component_ica_s1_.nii'['$i']'\
sef_sub014_component_ica_s1_.nii'['$i']'\
sef_sub015_component_ica_s1_.nii'['$i']'\
sef_sub016_component_ica_s1_.nii'['$i']'\
sef_sub017_component_ica_s1_.nii'['$i']'\
sef_sub018_component_ica_s1_.nii'['$i']'\
sef_sub019_component_ica_s1_.nii'['$i']'\
sef_sub020_component_ica_s1_.nii'['$i']'\
sef_sub021_component_ica_s1_.nii'['$i']'\
sef_sub022_component_ica_s1_.nii'['$i']'\
sef_sub023_component_ica_s1_.nii'['$i']'\
sef_sub024_component_ica_s1_.nii'['$i']'\
sef_sub025_component_ica_s1_.nii'['$i']'\
sef_sub026_component_ica_s1_.nii'['$i']'\
sef_sub027_component_ica_s1_.nii'['$i']'\
sef_sub028_component_ica_s1_.nii'['$i']'\
sef_sub029_component_ica_s1_.nii'['$i']'\
sef_sub030_component_ica_s1_.nii'['$i']'\
sef_sub031_component_ica_s1_.nii'['$i']'\
sef_sub032_component_ica_s1_.nii'['$i']'\
sef_sub033_component_ica_s1_.nii'['$i']'\
sef_sub034_component_ica_s1_.nii'['$i']'\
sef_sub035_component_ica_s1_.nii'['$i']'\
sef_sub001_component_ica_s2_.nii'['$i']'\
sef_sub002_component_ica_s2_.nii'['$i']'\
sef_sub003_component_ica_s2_.nii'['$i']'\
sef_sub004_component_ica_s2_.nii'['$i']'\
sef_sub005_component_ica_s2_.nii'['$i']'\
sef_sub006_component_ica_s2_.nii'['$i']'\
sef_sub007_component_ica_s2_.nii'['$i']'\
sef_sub008_component_ica_s2_.nii'['$i']'\
sef_sub009_component_ica_s2_.nii'['$i']'\
sef_sub010_component_ica_s2_.nii'['$i']'\
sef_sub011_component_ica_s2_.nii'['$i']'\
sef_sub012_component_ica_s2_.nii'['$i']'\
sef_sub013_component_ica_s2_.nii'['$i']'\
sef_sub014_component_ica_s2_.nii'['$i']'\
sef_sub015_component_ica_s2_.nii'['$i']'\
sef_sub016_component_ica_s2_.nii'['$i']'\
sef_sub017_component_ica_s2_.nii'['$i']'\
sef_sub018_component_ica_s2_.nii'['$i']'\
sef_sub019_component_ica_s2_.nii'['$i']'\
sef_sub020_component_ica_s2_.nii'['$i']'\
sef_sub021_component_ica_s2_.nii'['$i']'\
sef_sub022_component_ica_s2_.nii'['$i']'\
sef_sub023_component_ica_s2_.nii'['$i']'\
sef_sub024_component_ica_s2_.nii'['$i']'\
sef_sub025_component_ica_s2_.nii'['$i']'\
sef_sub026_component_ica_s2_.nii'['$i']'\
sef_sub027_component_ica_s2_.nii'['$i']'\
sef_sub028_component_ica_s2_.nii'['$i']'\
sef_sub029_component_ica_s2_.nii'['$i']'\
sef_sub030_component_ica_s2_.nii'['$i']'\
sef_sub031_component_ica_s2_.nii'['$i']'\
sef_sub032_component_ica_s2_.nii'['$i']'\
sef_sub033_component_ica_s2_.nii'['$i']'\
sef_sub034_component_ica_s2_.nii'['$i']'\
sef_sub035_component_ica_s2_.nii'['$i']'\
sef_sub001_component_ica_s3_.nii'['$i']'\
sef_sub002_component_ica_s3_.nii'['$i']'\
sef_sub003_component_ica_s3_.nii'['$i']'\
sef_sub004_component_ica_s3_.nii'['$i']'\
sef_sub005_component_ica_s3_.nii'['$i']'\
sef_sub006_component_ica_s3_.nii'['$i']'\
sef_sub007_component_ica_s3_.nii'['$i']'\
sef_sub008_component_ica_s3_.nii'['$i']'\
sef_sub009_component_ica_s3_.nii'['$i']'\
sef_sub010_component_ica_s3_.nii'['$i']'\
sef_sub011_component_ica_s3_.nii'['$i']'\
sef_sub012_component_ica_s3_.nii'['$i']'\
sef_sub013_component_ica_s3_.nii'['$i']'\
sef_sub014_component_ica_s3_.nii'['$i']'\
sef_sub015_component_ica_s3_.nii'['$i']'\
sef_sub016_component_ica_s3_.nii'['$i']'\
sef_sub017_component_ica_s3_.nii'['$i']'\
sef_sub018_component_ica_s3_.nii'['$i']'\
sef_sub019_component_ica_s3_.nii'['$i']'\
sef_sub020_component_ica_s3_.nii'['$i']'\
sef_sub021_component_ica_s3_.nii'['$i']'\
sef_sub022_component_ica_s3_.nii'['$i']'\
sef_sub023_component_ica_s3_.nii'['$i']'\
sef_sub024_component_ica_s3_.nii'['$i']'\
sef_sub025_component_ica_s3_.nii'['$i']'\
sef_sub026_component_ica_s3_.nii'['$i']'\
sef_sub027_component_ica_s3_.nii'['$i']'\
sef_sub028_component_ica_s3_.nii'['$i']'\
sef_sub029_component_ica_s3_.nii'['$i']'\
sef_sub030_component_ica_s3_.nii'['$i']'\
sef_sub031_component_ica_s3_.nii'['$i']'\
sef_sub032_component_ica_s3_.nii'['$i']'\
sef_sub033_component_ica_s3_.nii'['$i']'\
sef_sub034_component_ica_s3_.nii'['$i']'\
sef_sub035_component_ica_s3_.nii'['$i']'\
sef_sub001_component_ica_s4_.nii'['$i']'\
sef_sub002_component_ica_s4_.nii'['$i']'\
sef_sub003_component_ica_s4_.nii'['$i']'\
sef_sub004_component_ica_s4_.nii'['$i']'\
sef_sub005_component_ica_s4_.nii'['$i']'\
sef_sub006_component_ica_s4_.nii'['$i']'\
sef_sub007_component_ica_s4_.nii'['$i']'\
sef_sub008_component_ica_s4_.nii'['$i']'\
sef_sub009_component_ica_s4_.nii'['$i']'\
sef_sub010_component_ica_s4_.nii'['$i']'\
sef_sub011_component_ica_s4_.nii'['$i']'\
sef_sub012_component_ica_s4_.nii'['$i']'\
sef_sub013_component_ica_s4_.nii'['$i']'\
sef_sub014_component_ica_s4_.nii'['$i']'\
sef_sub015_component_ica_s4_.nii'['$i']'\
sef_sub016_component_ica_s4_.nii'['$i']'\
sef_sub017_component_ica_s4_.nii'['$i']'\
sef_sub018_component_ica_s4_.nii'['$i']'\
sef_sub019_component_ica_s4_.nii'['$i']'\
sef_sub020_component_ica_s4_.nii'['$i']'\
sef_sub021_component_ica_s4_.nii'['$i']'\
sef_sub022_component_ica_s4_.nii'['$i']'\
sef_sub023_component_ica_s4_.nii'['$i']'\
sef_sub024_component_ica_s4_.nii'['$i']'\
sef_sub025_component_ica_s4_.nii'['$i']'\
sef_sub026_component_ica_s4_.nii'['$i']'\
sef_sub027_component_ica_s4_.nii'['$i']'\
sef_sub028_component_ica_s4_.nii'['$i']'\
sef_sub029_component_ica_s4_.nii'['$i']'\
sef_sub030_component_ica_s4_.nii'['$i']'\
sef_sub031_component_ica_s4_.nii'['$i']'\
sef_sub032_component_ica_s4_.nii'['$i']'\
sef_sub033_component_ica_s4_.nii'['$i']'\
sef_sub034_component_ica_s4_.nii'['$i']'\
sef_sub035_component_ica_s4_.nii'['$i']'\
-prefix tmap.ic$i.$labels[$l].nii

@ l = $l + 1
#3dbucket -session maps -prefix 

end

