#!/bin/bash


# create spherical rois placed on roi activations, save in 4D file  
# create_rois.sh FMRI_DIR MEAN_MAP RADIUS THR MASKS OUTPUT

FMRI_DIR=$1  #/home/share/LeftHand/fMRI/processed
MEAN_MAP=$2  # mean_maps_0.nii.gz
RADIUS=$3    # 5
THR=$4       # 0.5
MASKS=$5     # /home/share/LeftHand/fMRI/processed/all_masks/fc_masks/*.nii.gz
OUTPUT=$6
cd $FMRI_DIR

rm max_coords.txt spheres.nii.gz

COUNT=1
for mask in $MASKS; do

    NAME=`basename $mask .nii.gz`
    echo $NAME
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
        fslmerge -t spheres_${i} spheres_${i} mymask.nii.gz
    else 
        cp mymask.nii.gz spheres_${i}.nii.gz 
    fi
    
    COUNT=$[$COUNT + 1]       
    fslmaths spheres.nii.gz -Tmean -mul $[$COUNT - 1] mean_spheres.nii.gz

done
mv spheres.nii.gz $OUTPUT
rm mymap.nii.ngz mymask.nii.gz
