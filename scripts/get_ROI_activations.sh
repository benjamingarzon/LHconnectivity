
#for i in 5 8; do for s in cope zstat; do ./get_ROI_activations.sh $i $s; done; done

RADIUS=$1 #5
ATLAS=/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/atlas_lefthand${RADIUS}.nii.gz
MAPNAME=$2 #zstat

INPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_lh_all.nii.gz
OUTPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_lh_all_ROI${RADIUS}.csv
mri_segstats --i $INPUT --seg $ATLAS --avgwf $OUTPUT --excludeid 0

INPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_rh_all.nii.gz
OUTPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_rh_all_ROI${RADIUS}.csv
mri_segstats --i $INPUT --seg $ATLAS --avgwf $OUTPUT --excludeid 0


INPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_lh-rh_all.nii.gz
OUTPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_lh-rh_all_ROI${RADIUS}.csv
mri_segstats --i $INPUT --seg $ATLAS --avgwf $OUTPUT --excludeid 0


INPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_rh-lh_all.nii.gz
OUTPUT=/home/share/LeftHand/LHconnectivity/Stockholm/${MAPNAME}_rh-lh_all_ROI${RADIUS}.csv
mri_segstats --i $INPUT --seg $ATLAS --avgwf $OUTPUT --excludeid 0

