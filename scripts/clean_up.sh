cd /home/share/LeftHand/LHconnectivity

for FOLDER in Berlin Stockholm; do
rm $FOLDER/sub*/sess-*/funcproc_*/*/fMRI_mcf.nii.gz
#rm $FOLDER/sub*/sess-*/funcproc_*/*/fMRI_mcf*par
rm -r $FOLDER/sub*/sess-*/funcproc_*/*/nilearn_cache
rm $FOLDER/sub*/sess-*/funcproc_*/*/fMRI_mni.nii.gz
rm $FOLDER/sub*/sess-*/funcproc_*/*/fMRI_denoised_mni.nii.gz
#rm -r $FOLDER/sub*/sess-*/icaproc_*
#rm -r $FOLDER/sub*/sess-*/funcproc_*/*/analysis.feat*
rm -r $FOLDER/sub*/sess-*/funcproc_*/*/analysis.feat/featquery*

done

#fslmerge -t all_mean ./sub-LH*/sess-*/funcproc_task-localizer/run-*/fMRI_mni_mean.nii.gz
