#!/bin/bash


# create spherical rois placed on roi activations, save in 4D file  
# create_rois.sh FMRI_DIR MEAN_MAP RADIUS THR MASKS OUTPUT ROI_LIST
# ./create_rois.sh /home/share/LeftHand/fMRI/processed mean_maps_0.nii.gz 5 0.5 "/home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz" /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand5 /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand5.txt
# ./create_rois.sh /home/share/LeftHand/fMRI/processed mean_maps_0.nii.gz 8 0.5 "/home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz" /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand8 /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand8.txt
./create_rois.sh /home/share/LeftHand/fMRI/processed mean_maps_0.nii.gz 10 0.5 "/home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz" /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10 /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10.txt

FMRI_DIR=$1  #/home/share/LeftHand/fMRI/processed
MEAN_MAP=$2  # mean_maps_0.nii.gz
RADIUS=$3    # 10
THR=$4       # 0.5
MASKS=$5     # /home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz
OUTPUT=$6    # /home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10
ROI_LIST=$7  #/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10.txt


FMRI_DIR=/home/share/LeftHand/fMRI/processed
MEAN_MAP=mean_maps_0.nii.gz
RADIUS=10
THR=0.5
MASKS=/home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz
OUTPUT=/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10
ROI_LIST=/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand10.txt


cd $FMRI_DIR

rm -r $OUTPUT
mkdir $OUTPUT
rm max_coords.txt spheres.nii.gz

rm $ROI_LIST
COUNT=1
for mask in $MASKS; do

    NAME=`basename $mask .nii.gz`
    echo $NAME
    echo $NAME >> $ROI_LIST
     
    MASKMAX=`fslstats $mask -R | cut -f2 -d' '`
    fslmaths $mask -div $MASKMAX -thr $THR -bin -mul $MEAN_MAP mymap.nii.gz
    MAX=`fslstats mymap.nii.gz -x`
    echo $NAME $MAX >> max_coords.txt
    X=`echo $MAX | cut -d' ' -f1`
    Y=`echo $MAX | cut -d' ' -f2`
    Z=`echo $MAX | cut -d' ' -f3`

    fslmaths mymap.nii.gz -roi $X 1 $Y 1 $Z 1 0 1 -bin mymap.nii.gz
    
    fslmaths mymap.nii.gz -kernel sphere $RADIUS -fmean -bin mymask.nii.gz

    if [ -e spheres.nii.gz ]; then
        fslmaths mymask.nii.gz -mul $COUNT mymask.nii.gz 
        fslmerge -t spheres spheres mymask.nii.gz
    else 
        cp mymask.nii.gz spheres.nii.gz 
    fi
    
    mv mymask.nii.gz $OUTPUT/$NAME.nii.gz
    COUNT=$[$COUNT + 1]       
    fslmaths spheres.nii.gz -Tmean -mul $[$COUNT - 1] mean_spheres.nii.gz

done
fslmaths spheres.nii.gz -Tmax $OUTPUT.nii.gz
rm mymap.nii.gz mymask.nii.gz


