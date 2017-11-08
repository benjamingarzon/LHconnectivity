#!/bin/sh
# organize data for Stockholm Experiment
# leading zeros and session 2 to 4 change

OVERWRITE=$1
WD=/home/share/LeftHand/LHconnectivity/Stockholm
FUNC_DIR=/home/share/LeftHand/fMRI/processed
STRUCT_DIR=/home/share/LeftHand/Structural/
LOG_FILE=$WD/logs/organize_data_Stockholm.log
NVOLS=199
TASK=localizer
# these subjects need to be renamed
RENAMELIST="LH2002.2 LH2003.4 LH2005.4 LH2007.2 LH2008.2 LH2009.4 LH2011.4 \ 
LH2012.2 LH2013.4 LH2014.2 LH2015.2 LH2017.4 LH2018.2"

# rename structural folders that were named differently for the structural
# and functional data 

if [ "$OVERWRITE" -eq "1" ]; then
rm -r $WD
mkdir $WD
mkdir $WD/logs
else
rm $WD/logs/* 
fi

cd $FUNC_DIR

for subject in LH*.*; do 

NAME=`echo ${subject} | cut -d'.' -f1`
SESS=`echo ${subject} | cut -d'.' -f2`
printf -v SESS "%02d" $SESS

mkdir $WD/sub-${NAME}
mkdir $WD/sub-${NAME}/sess-${SESS}

# fetch fMRI
for RUN in 1 2; do

if [ -e $FUNC_DIR/${subject}/run${RUN} ]; 
then
if [ -e $FUNC_DIR/${subject}/run${RUN}/fMRI.nii.gz ]; 
then
if [ `fslval $FUNC_DIR/${subject}/run${RUN}/fMRI.nii.gz dim4` -eq "$NVOLS" ]; 
then

mkdir $WD/sub-${NAME}/sess-${SESS}/func/
ln -s $FUNC_DIR/${subject}/run${RUN}/fMRI.nii.gz \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK}_run-${RUN}.nii.gz
mkdir $WD/sub-${NAME}/sess-${SESS}/func/\
sub-${NAME}_task-${TASK}_run-${RUN}_events

ln -s $FUNC_DIR/${subject}/run${RUN}/EV1.csv \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK}_run-${RUN}_events/\
EV1.csv
ln -s $FUNC_DIR/${subject}/run${RUN}/EV2.csv \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK}_run-${RUN}_events/\
EV2.csv

else
echo "${subject}/run${RUN}/fMRI.nii.gz: wrong number of vols." >> $LOG_FILE
fi
else
echo "${subject}/run${RUN}/fMRI.nii.gz does not exist." >> $LOG_FILE
fi
else
echo "${subject}/run${RUN} does not exist." >> $LOG_FILE
fi

done

# fetch anatomical


mkdir $WD/sub-${NAME}/sess-${SESS}/anat/

if [ "`echo $RENAMELIST| grep $subject|wc -l`" -eq "1" ]; 
then
subjectstruct="`echo $subject | cut -d'.' -f1`.3"
else
subjectstruct=$subject
fi


if [ -e $STRUCT_DIR/${subjectstruct}/unprocessed/3T/T1w_MPR1/\
${subjectstruct}_3T_T1w_MPR1.nii.gz ]; then
  ln -s $STRUCT_DIR/${subjectstruct}/unprocessed/3T/T1w_MPR1/\
${subjectstruct}_3T_T1w_MPR1.nii.gz $WD/sub-${NAME}/sess-${SESS}/\
anat/sub-${NAME}_T1w.nii.gz
else 
echo "${subjectstruct}/unprocessed/3T/T1w_MPR1/"\
"${subjectstruct}_3T_T1w_MPR1.nii.gz does not exist." >> $LOG_FILE
fi


if [ -e $STRUCT_DIR/${subjectstruct}/unprocessed/3T/T2w_SPC1/\
${subjectstruct}_3T_T2w_SPC1.nii.gz ]; then
  ln -s $STRUCT_DIR/${subjectstruct}/unprocessed/3T/T2w_SPC1/\
${subjectstruct}_3T_T2w_SPC1.nii.gz $WD/sub-${NAME}/sess-${SESS}/\
anat/sub-${NAME}_T2w.nii.gz
else 
echo "${subjectstruct}/unprocessed/3T/T2w_SPC1/"\
"${subjectstruct}_3T_T2w_SPC1.nii.gz does not exist." >> $LOG_FILE
fi

done



