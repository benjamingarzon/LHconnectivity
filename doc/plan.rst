PLAN

MEETING
show comparison writing / resting in Berlin
possibly just weaker than Stockholm
show traces, flat
show clustering, show prediction

show clustering Stockholm
denoised/ not denoised
full vs lw correlations

CHECK FISHER transform

Prediction> Stockholm / Berlin / Stockholm - Berlin and interpret connections
OBSERVE AND CLEAN FITS
INSTALL CARET
X/VALIDATE PREDICTION



Add Prefrontal and midbrain

Do clustering
RUN*SESSION interaction
RUN difference

Compare writing and resting matrices
How many TRs to remove

Before clustering, remove motion, group difference etc...

Run denoised data
Compcor ?



--------------------------------------

WORKDIR = /home/share/LeftHand/LHConnectivity

COMBINATIONS FOR OUTPUT FILES
connectivity: full/lw
process: aromaaggr/aromanoaggr/noaroma
label: rest/all

/ change do preprocessing to get any file
am I using NVOLS / no
alter .fsf and check crop size


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
x - clean intermediate files


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


- Remove dicoms from repeated acquistitions
x - clean_dicom_Berlin. 

QUESTIONS
 - How does the connectivity pattern predict behaviour?  
  
  # remove nilearn_cache...
  
