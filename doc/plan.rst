PLAN

  WORKDIR = /home/share/LeftHand/LHConnectivity
# change events by data

** organize_data_Stockholm.sh
**  organize_data_Berlin.sh
 - Rename data as BIDS

** do_preprocessing.py  
 - Average all T1s for each subject longitudinal 
 - Run FIX or AROMA? 
 - Register 
 - Get confounds  
   - Friston 24
   - do I need to remove white matter and CSF signal ?  


** Estimate fMRI task data
 parallelize?
 clean intermediate files


 - Get striatum parcellation - Yeo --

 - Get individual ROIs and create parcellation file
 create_rois.sh
   - from localizer
   - from neurosynth? ICA, anatomical 
   - from contrast map - early vs late : to find e.g. PFC

 - Extract connectivity data matrices 
  - small smoothing
  - which estimator ? Ridge seems good! / Ledoit_Wolf
 
 - Cluster trajectories
 # do analysis

 - Extract substantia nigra, use MBST?


QUESTIONS
 - How does the connectivity pattern predict behaviour?  
  
  # remove nilearn_cache...
  
