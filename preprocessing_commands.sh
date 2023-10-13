# Zeropad subjec ID numbers -- subs_nopad must be generated manually
for sub in $(cat subs_nopad.txt)
do 
    prefix=$(echo $sub | cut -c -7)
    suffix=$(echo $sub | cut -c 8-)
    new_name=$(echo "$prefix"0"$suffix")
    echo $new_name

    mv $sub $new_name
    rename "s/$sub/$new_name/g"  $(find . -type f | grep $sub)
    
done

# Get old names 
ls -d sub-* > subs.txt

for sub in $(cat subs.txt)
do 
    new_name=$sub
    old_name=$(cat "$sub"/*.tsv | awk '{print $3}' | tail -n -1)
    echo "$new_name  ----> $old_name"   
done > oldnames.txt

# Associate the current names with the ones in metadata and rename the lesion masks
grep study_id metadata.tsv > participants.tsv
for sub in $(cat subs.txt)
do 
    group=$(echo $sub | cut -c 5- | cut -c -3)
    id=$(echo $sub | cut -c 5- | cut -c 4-)

    if [[ "$group" == "CON" ]]; then
        associated_name=$(echo "AMC_"$id"_FCS")
        echo $sub $associated_name 
        #cat metadata.tsv | grep -w $associated_name | sed "s/$associated_name/$sub/g" >> participants.tsv
    elif [[ "$group" == "PAT" ]]; then
        associated_name=$(echo "FCS_"$id"")
        echo $sub $associated_name
        #cat metadata.tsv | grep -w $associated_name | sed "s/$associated_name/$sub/g" >> participants.tsv

        #lesion_newname=$(find LesionMaks_222/NiftiLesions/ -type f | grep FCS_027 | sed "s/$associated_name/$sub/g")
        
        filename=$(find LesionMaks_222/NiftiLesions/ -type f | grep $associated_name)
        if test -f "$filename";
        then
            rename  "s/$associated_name/$sub/g" $(find LesionMaks_222/NiftiLesions/ -type f | grep $associated_name)
        else
            echo "$filename has not been found"
        fi
        #rename  "s/$associated_name/$sub/g" $(find LesionMaks_222/NiftiLesions/ -type f | grep $associated_name)

    fi
done


# remove backup files and empty folders. 
rm $(find . -type f | grep backup)
find . -type d -empty -delete




# remove empty rows from participants.tsv
awk 'BEGIN{FS="\t";OFS="\t"}
    {
      for(i=1;i<=(NF-1);i++)
      {
        if($i == ""){
         $i="n/a"
        }
      }
      print
    }' ../stroke_after_Cemals_code/participants.tsv > temp.tsv

    mv temp.tsv participants.tsv

## Remove the dash character in between task name 

for file in $(find stroke_BIDS/ -type d |grep control-2); do  mv $file $(echo $file | sed "s/control-2/control2/g") ;done
for file in $(find stroke_BIDS/ -type f |grep control-2); do  mv $file $(echo $file | sed "s/control-2/control2/g") ;done

for file in $(find stroke_BIDS/ -type d |grep followup-2); do  mv $file $(echo $file | sed "s/followup-2/followup2/g") ;done
for file in $(find stroke_BIDS/ -type f |grep followup-2); do  mv $file $(echo $file | sed "s/followup-2/followup2/g") ;done


## Create a copy of the dataset that only includes the first runs (the same was repeated for second and third runs)
cp -r stroke_BIDS/ stroke_BIDS_firs_sessions
cd stroke_BIDS_firs_sessions/
find . -type  f | grep control2
rm $(find . -type  f | grep control2 )
find . -type  f | grep followup
rm $(find . -type  f | grep followup)
find . -type d -empty -print
find . -type d -empty -delete
# remove the folders that have no runs after the operation above 
rm -r $(find . -type d -links 2 |sort | grep -v anat |grep -v dwi |grep -v func )


# First we need to send the files to ares efficiently 

# Get the folder list 
ls -d -- */ > dir_list.txt

# Divide the folder list into batches
cd sources/parallelcopy 
split -l 17 ../../dir_list.txt
rm dir_list.txt
cd -
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xaa) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xab) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xac) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xad) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xae) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xaf) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xag) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xah) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xai) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xaj) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xak) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/
rsync -av -R $(cat /media/koba/MULTIBOOT/stroke_BIDS/source/parallelcopy/xal) plgkoba@ares.cyfronet.pl:/net/ascratch/people/plgkoba/stroke_BIDS/

## Create a singularity image for fmriprep on the local computer first 
sudo docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /home/koba/Desktop/Stroke/:/output \
    singularityware/docker2singularity \
    nipreps/fmriprep:latest

## Copy the image to Ares
scp /path/to/image plgkoba@ares.cyfronet.pl:/$SCRATCH

## Create a batch file for the singularity

#!/bin/bash -l 
## task_name 
#SBATCH -J test 
## num nodes 
#SBATCH -N 1 
## tasks per node 
#SBATCH --ntasks-per-node=1 
## cpus per taskq 
#SBATCH --cpus-per-task=16
## memory allocated per cpu 
#SBATCH --mem-per-cpu=4GB 
## max time 
#SBATCH --time=35:00:00 
## grant name 
#SBATCH -A plgsano4-cpu 
## partitionq 
#SBATCH --partition plgrid
## gpus allocation 
##SBATCH --gpus-per-node=2 
## output file path 
#SBATCH --output="/net/ascratch/people/plgkoba/output1.out" 
## error file path 
#SBATCH --error="/net/ascratch/people/plgkoba/error1.err" 
## loading modules that are not in the environment 
## module add cudnn/8.4.1.50-cuda-11.6.0 
## bash commands to execute 

#76 81 82 87 89

cd $SCRATCH
date

for sub in PAT084 PAT087 PAT088 PAT090 PAT092 PAT097 PAT098 PAT099 PAT100 PAT101 PAT102 PAT104 PAT106 PAT107 PAT108 PAT109 PAT112 PAT113 PAT115 PAT116 PAT117 PAT119 PAT120 PAT122
do
export APPTAINERENV_TEMPLATEFLOW_HOME=$SCRATCH/templateflow
mkdir $SCRATCH/temp/sub-"$sub"_temp
#sessions=$(cat sub-"$sub"/sub-"$sub"_sessions.tsv | awk '{print $1}' | grep -v session)

apptainer run -B $SCRATCH/ --cleanenv fmriprep.simg  \
 $SCRATCH/stroke_BIDS_second_sessions   $SCRATCH/stroke_BIDS_second_sessions/derivatives/fmriprep \
 participant --participant_label "$sub" --skip_bids_validation  \
 --fs-no-reconall --fs-license-file $SCRATCH/license.txt --nprocs 16 --work-dir $SCRATCH/temp/sub-"$sub"_temp
rm -r $SCRATCH/temp/sub-"$sub"_temp


done
date


## Submit the above code to the HPC
sbatch fmriprep_batch1.sh




### XCP commands
sublist=$(ls derivatives/fmriprep/*html | cut -d / -f 3 | cut -d . -f 1 | cut -d - -f 2 )
for sub in $(echo $sublist); do
sudo docker run --rm -it -v $fmriprep_dir:/data/ -v $out_dir:/out pennlinc/xcp_d:latest /data /out participant --nuisance-regressors $method --clean-workdir --smoothing 6  --fd-thresh 0.5 --min_coverage 0.001 --combineruns --participant-label sub-"$sub"
done
## Calculate the time lag maps 
for subject in $(cat /home/koba/Desktop/Stroke/scripts/subjects.txt)
do
	mkdir rapidtide/"$subject"
	rapidtide acompcor/xcp_d/"$subject"/ses-control/func/"$subject"_ses-control_task-rest_space-MNI152NLin2009cAsym_desc-denoisedSmoothed_bold.nii.gz rapidtide/"$subject"/"$subject" --searchrange -15 15
done
