cd /home/share/LeftHand/LHconnectivity

for FOLDER in Berlin Stockholm; do
rm $FOLDER/sub*/sess-*/funcproc_*/*/fMRI_mcf.nii.gz
rm -r $FOLDER/sub*/sess-*/funcproc_*/*/nilearn_cache
done

fslmerge -t all_mean ./sub-LH*/sess-*/funcproc_task-localizer/run-*/fMRI_mni_mean.nii.gz
