PLAN

WORKDIR = /home/share/LeftHand/LHConnectivity

COMMBINATIONS FOR OUTPUT FILES
connectivity: full/lw
process: aromaaggr/aromanoaggr/noaroma
label: rest/all

# confounds in aroma or in niftimasker
** organize_data_Stockholm.sh
**  organize_data_Berlin.sh
x - Rename data as BIDS

** do_preprocessing.py  
x - Average all T1s for each subject longitudinal 
x - Run FIX or AROMA? 
x - Register 
x - Get confounds  
x - Friston 24
   - do I need to remove white matter and CSF signal ?  


** Estimate fMRI task data
x - parallelize?
 clean intermediate files


x - Get striatum parcellation - Yeo --

 - Get individual ROIs and create parcellation file
x create_rois.sh
   - from localizer
   - from neurosynth? ICA, anatomical 
   - from contrast map - early vs late : to find e.g. PFC

 - Extract connectivity data matrices 
  - small smoothing
x  - which estimator ? Ridge seems good! / Ledoit_Wolf
x  - output a few estimators : full / ledoit
 
 - Cluster trajectories
 # do analysis

 - Extract substantia nigra, use MBST?


QUESTIONS
 - How does the connectivity pattern predict behaviour?  
  
  # remove nilearn_cache...
  
