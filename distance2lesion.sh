for sub in $(cat source/subjects.txt); do file=$(find NiftiLesions/  -type f | grep $sub); echo $file  ;done > masklist.txt


for file in $(cat masklist.txt); do id=$(echo  $file | cut  -d "/" -f 2 | cut -d "_" -f 1); 3dresample -input $file -master source/arterial_resampled.nii.gz -prefix  NiftiLesions/resampled/"$id"_lesion_resampled.nii.gz; done






for file in $(cat masklist.txt); do id=$(echo  $file | cut  -d "/" -f 2 | cut -d "_" -f 1);3dmaskdump  -mask NiftiLesions/resampled/"$id"_lesion_resampled.nii.gz -xyz -noijk  NiftiLesions/resampled/"$id"_lesion_resampled.nii.gz > 	distance2lesion/"$id"_lesion_coords.txt;done

for file in $(cat masklist.txt); do id=$(echo  $file | cut  -d "/" -f 2 | cut -d "_" -f 1);3dUndump -prefix "$id"_distance2lesion.nii.gz  -xyz -datum float -master /media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/yeo400_resampled.nii.gz distance2lesion/"$id"_distance2lesion.txt -overwrite;done
# Euclidean distance is calculated in matlab and saved back as nifti 

# This is after matlab code
3dbucket -prefix distance2lesion/all_distances.nii.gz $(find distance2lesion/ -type f | grep gz | sort)
