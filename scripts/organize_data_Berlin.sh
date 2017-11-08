#!/bin/sh
# organize data for Berlin Experiment

OVERWRITE=$1
WD=/home/share/LeftHand/LHconnectivity/Berlin

SUBJECTS="Time1001A Time1001C Time1001E Time1001G Time1001I Time1001K Time1001M Time1001O Time1001post Time1001pre Time1002interim Time1002post Time1002pre Time1003A Time1003C Time1003E Time1003G Time1003I Time1003K Time1003M Time1003post Time1003pre Time1004A Time1004C Time1004G Time1004I Time1004K Time1004O Time1004post Time1004pre Time1005A Time1005C Time1005G Time1005I Time1005K Time1005M Time1005O Time1005post Time1005pre Time1006A Time1006C Time1006E Time1006G Time1006I Time1006K Time1006O Time1006post Time1006pre Time1007A Time1007C Time1007E Time1007G Time1007I Time1007K Time1007M Time1007O Time1007post Time1007pre Time1008A Time1008C Time1008E Time1008G Time1008I Time1008K Time1008M Time1008O Time1008post Time1008pre Time1009A Time1009C Time1009E Time1009G Time1009I Time1009K Time1009M Time1009O Time1009post Time1009pre Time1010A Time1010C Time1010E Time1010G Time1010I Time1010K Time1010M Time1010O Time1010post Time1010pre Time1011A Time1011C Time1011E Time1011G Time1011I Time1011K Time1011M Time1011O Time1011post Time1011pre Time1012A Time1012C Time1012E Time1012G Time1012I Time1012K Time1012M Time1012O Time1012post Time1012pre Time1013A Time1013C Time1013E Time1013G Time1013I Time1013K Time1013M Time1013O Time1013post Time1013pre Time1014A Time1014C Time1014E Time1014G Time1014I Time1014K Time1014M Time1014O Time1014post Time1014pre Time1015A Time1015C Time1015E Time1015G Time1015I Time1015K Time1015M Time1015O Time1015post Time1015pre Time2001interim Time2001post Time2001pre Time2002interim Time2002post Time2002pre Time2003interim Time2003post Time2003pre Time2004interim Time2004post Time2004pre Time2006interim Time2006post Time2006pre Time2007interim Time2007post Time2007pre Time2008A Time2008C Time2008E Time2008G Time2008I Time2008K Time2008M Time2008O Time2008post Time2008pre Time2009post Time2009pre Time2011interim Time2011post Time2011pre Time2013interim Time2013post Time2013pre Time2014interim Time2014post Time2014pre Time2015interim Time2015post Time2015pre Time2016interim Time2016post Time2016pre Time2017interim Time2017post Time2017pre Time2018pre Time2019interim Time2019post Time2019pre Time2020interim Time2020post Time2020pre"

FUNC_DIR1=/home/share/LeftHand/LHconnectivity/Berlin_raw/resting
FUNC_DIR2=/home/share/LeftHand/LHconnectivity/Berlin_raw/writing

STRUCT_DIR=/home/share/LeftHand/LHconnectivity/Berlin_raw/MPRAGE
LOG_FILE=$WD/logs/organize_data_Berlin.log
NVOLS1=140 #150
NVOLS2=260 #270

TASK1=resting
TASK2=writing

if [ "$OVERWRITE" -eq "1" ]; then
rm -r $WD
mkdir $WD
mkdir $WD/logs
fi

for subject in $SUBJECTS; do

NAME=`echo ${subject:4:4}` 
SESS_BER=`echo ${subject:8} `
echo $NAME $SESS_BER

if [ "$SESS_BER" == "pre" ];    then SESS=01; fi
if [ "$SESS_BER" == "A" ];       then SESS=02; fi
if [ "$SESS_BER" == "B" ];       then SESS=03; fi
if [ "$SESS_BER" == "C" ];       then SESS=04; fi
if [ "$SESS_BER" == "D" ];       then SESS=05; fi
if [ "$SESS_BER" == "E" ];       then SESS=06; fi
if [ "$SESS_BER" == "F" ];       then SESS=07; fi
if [ "$SESS_BER" == "G" ];       then SESS=08; fi
if [ "$SESS_BER" == "H" ];       then SESS=09; fi
if [ "$SESS_BER" == "I" ];       then SESS=10; fi
if [ "$SESS_BER" == "J" ];       then SESS=11; fi
if [ "$SESS_BER" == "K" ];       then SESS=12; fi
if [ "$SESS_BER" == "L" ];       then SESS=13; fi
if [ "$SESS_BER" == "M" ];       then SESS=14; fi
if [ "$SESS_BER" == "N" ];       then SESS=15; fi
if [ "$SESS_BER" == "O" ];       then SESS=16; fi
if [ "$SESS_BER" == "P" ];       then SESS=17; fi
if [ "$SESS_BER" == "post" ];   then SESS=18; fi

if [ "$SESS_BER" == "interim" ]; then SESS=09; fi

mkdir $WD/sub-${NAME}
mkdir $WD/sub-${NAME}/sess-${SESS}

# fetch fMRI - resting state
for RUN in 1 2; do
MYDIR=`echo $FUNC_DIR1/${subject}*/rest${RUN}`
echo $MYDIR
if [ -e "$MYDIR" ]; 
then
if [ ! -e "$MYDIR/fMRI.nii.gz" ]; then
mri_convert $MYDIR/*-0001.dcm $MYDIR/fMRI.nii.gz
fi
if [ -e "$MYDIR/fMRI.nii.gz" ]; 
then

MYVOLS=`fslval $FUNC_DIR1/${subject}*/rest${RUN}/fMRI.nii.gz dim4`
if [ "$MYVOLS" -le "$NVOLS1" ]; 
then
#echo "$FUNC_DIR1/${subject}/rest${RUN}/fMRI.nii.gz: right number of vols."    
#else
echo "$MYDIR/fMRI.nii.gz: too few vols: $MYVOLS" >> \
$LOG_FILE
fi


mkdir $WD/sub-${NAME}/sess-${SESS}/func/
ln -s $FUNC_DIR1/${subject}*/rest${RUN}/fMRI.nii.gz \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK1}_run-${RUN}.nii.gz
#mkdir $WD/sub-${NAME}/sess-${SESS}/func/\
#sub-${NAME}_task-${TASK}_run-${RUN}_events

else
echo "$MYDIR/fMRI.nii.gz does not exist." >> $LOG_FILE
fi
else
echo "$MYDIR does not exist." >> $LOG_FILE
fi

done


# fetch fMRI - writing task
RUN=1
MYDIR=`echo $FUNC_DIR2/${subject}*`
echo $MYDIR
if [ -e "$MYDIR" ]; 
then

if [ ! -e "$MYDIR/fMRI.nii.gz" ]; then
mri_convert $MYDIR/*-0001.dcm $MYDIR/fMRI.nii.gz
fi
if [ -e "$MYDIR/fMRI.nii.gz" ]; 
then

MYVOLS=`fslval $FUNC_DIR2/${subject}*/fMRI.nii.gz dim4`
if [ "$MYVOLS" -le "$NVOLS2" ]; 
then
#echo "$FUNC_DIR2/${subject}/fMRI.nii.gz: right number of vols."     
#else
echo "$MYDIR/fMRI.nii.gz: too few vols: $MYVOLS" >> \
$LOG_FILE
fi


mkdir $WD/sub-${NAME}/sess-${SESS}/func/

ln -s $MYDIR/fMRI.nii.gz \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK2}_run-${RUN}.nii.gz
mkdir $WD/sub-${NAME}/sess-${SESS}/func/\
sub-${NAME}_task-${TASK2}_run-${RUN}_events

python create_evs_Berlin.py $FUNC_DIR2/onsets/${subject}*.txt \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK2}_run-${RUN}_events/\
EV1.csv \
$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK2}_run-${RUN}_events/\
EV2.csv

#EV1.csv
#$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK2}_run-${RUN}_events/\
#EV1.csv
#ln -s $FUNC_DIR2/${subject}/run${RUN}/EV2.csv \
#$WD/sub-${NAME}/sess-${SESS}/func/sub-${NAME}_task-${TASK2}_run-${RUN}_events/\
#EV2.csv

else
echo "$MYDIR/fMRI.nii.gz does not exist." >> $LOG_FILE
fi
else
echo "$MYDIR/ does not exist." >> $LOG_FILE
fi


# fetch anatomical

mkdir $WD/sub-${NAME}/sess-${SESS}/anat/

#fslchfiletype NIFTI_GZ $STRUCT_DIR/${subject}*.img
MYFILE=`echo $STRUCT_DIR/${subject}*.nii.gz`
if [ -e "$MYFILE" ]; then
  ln -s $MYFILE $WD/sub-${NAME}/sess-${SESS}/\
anat/sub-${NAME}_T1w.nii.gz
else 
echo "$MYFILE does not exist." >> $LOG_FILE
fi

done
