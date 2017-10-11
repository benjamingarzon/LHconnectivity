# organize data for Stockholm

WD=/home/share/LeftHand/LHConnectivity/Stockholm
FUNC_DIR=cd /home/share/LeftHand/fMRI/processed
STRUCT_DIR=cd /home/share/LeftHand/fMRI/processed
NVOLS=199
TASK=localizer

cd $FUNC_DIR

for subject in LH*.*; do 

NAME=`echo $subject | cut -d'.' -f1`
SESS=`echo $subject | cut -d'.' -f2`

mkdir $WD/sub_$NAME
mkdir $WD/sub_$NAME/sess_$SESS

# fetch fMRI
for RUN in 1 2; do

if [ -e $FUNC_DIR/$subject/run$RUN ]; then
if [ -e $FUNC_DIR/$subject/run$RUN/fMRI.nii.gz ]; then
if [ -e fslval $FUNC_DIR/$subject/run$RUN/fMRI.nii.gz dim4 -eq "$NVOLS" ]; then

mkdir $WD/sub-$NAME/sess-$SESS/func/
ln -s $FUNC_DIR/$subject/run$RUN/fMRI.nii.gz $WD/sub-$NAME/sess-$SESS/func/sub-$NAME_task-$TASK_run-$RUN.nii.gz
mkdir $WD/sub-$NAME/sess-$SESS/func/sub-$NAME_task-$TASK_run-$RUN_events

ln -s $FUNC_DIR/$subject/run$RUN/EV1.csv $WD/sub-$NAME/sess-$SESS/func/sub-$NAME_task-$TASK_run-$RUN_events/EV1.csv
ln -s $FUNC_DIR/$subject/run$RUN/EV2.csv $WD/sub-$NAME/sess-$SESS/func/sub-$NAME_task-$TASK_run-$RUN_events/EV2.csv

fi
fi
fi

if [ -e $STRUCT_DIR/$subject/run$RUN ]; then


done



