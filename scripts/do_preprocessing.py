import sys, glob, os
import numpy as np
from joblib import Parallel, delayed
import json

sys.path.append('..')
import LHconnectivity.preprocess_funcs as pf
import LHconnectivity as lh

PARAMS_FILE = os.path.join(lh.__path__[0], 'data', sys.argv[1])
#PARAMS_FILE = os.path.join(lh.__path__[0], 'data', 'parameters_Stockholm.json')
with open(PARAMS_FILE) as data_file:    
    parameters = json.load(data_file)

WD = parameters['WD']
NJOBS = parameters['NJOBS']
CROP_SIZE = parameters['CROP_SIZE']
TR = parameters['TR']
mni = (parameters['MNI_FILENAME'], parameters['MNI_BRAIN_FILENAME'])
AROMA_PATH = parameters['AROMA_PATH']
FWHM = parameters['RS_FWHM']
LOW_PASS = parameters['LOW_PASS']
HIGH_PASS = parameters['HIGH_PASS']
FULL_FC_FILENAME = parameters['FULL_FC_FILENAME']
LW_FC_FILENAME = parameters['LW_FC_FILENAME']
AGGR_TYPE = parameters['AGGR_TYPE']
func_mni_filename = parameters['FUNC_MNI_FILENAME']
task = parameters['TASK']
rest = parameters['REST']
suffix = parameters['SUFFIX']
fsf_file = parameters['FSF_FILE']
atlas_filename = parameters['ATLAS_FILENAME']
atlas_names = parameters['ATLAS_NAMES']
mask_filename = parameters['MASK_FILENAME']
label = parameters['LABEL']
acq_order = parameters['ACQ_ORDER']
NREMOVE = parameters['NREMOVE']
MINVOLS = parameters['MINVOLS']
NVOLS = parameters['NVOLS']
ALFF_filename = parameters['ALFF_FILENAME']

os.chdir(WD)
#subjects = glob.glob('sub-LH100[1-5]')
subjects = glob.glob('sub-*')

# get average template
#Parallel(n_jobs = NJOBS)(delayed(pf.average_structurals)(WD, subject, CROP_SIZE) 
#    for subject in subjects) 

# register template
#Parallel(n_jobs = NJOBS)(delayed(pf.register_template)
#    (os.path.join(WD, subject), 'template.nii.gz', mni) 
#    for subject in subjects) 

# rest
# motion_correction
Parallel(n_jobs = NJOBS)(delayed(pf.motion_correction)(WD, subject, rest, TR,
    acq_order) for subject in subjects) 

# register functional data
Parallel(n_jobs = NJOBS)(delayed(pf.register_funcs)(WD, subject, rest, TR, mni) 
    for subject in subjects) 

#if rest != task:
# motion_correction
#    Parallel(n_jobs = NJOBS)(delayed(pf.motion_correction)(WD, subject, task, 
#        TR, acq_order) for subject in subjects) 

# register functional data
#    Parallel(n_jobs = NJOBS)(delayed(pf.register_funcs)(WD, subject, task, TR, 
#        mni) for subject in subjects) 

# run AROMA - only in rest
#Parallel(n_jobs = NJOBS)(delayed(pf.run_aroma)(WD, subject, rest, TR, 
#    AROMA_PATH, AGGR_TYPE, mni) 
#    for subject in subjects) 


#pf.get_ALFFs(WD, subjects, rest, mask_filename, TR, FWHM, 
#   LOW_PASS, HIGH_PASS, func_mni_filename, ALFF_filename + '.txt', 
#   ALFF_filename + '.nii.gz')

stophere

# process fMRI task data
Parallel(n_jobs = NJOBS)(delayed(pf.process_fMRI)(WD, subject, task, suffix, 
    fsf_file, mni)
    for subject in subjects) 

stophere

pf.get_percent_signal(WD, subjects, task, suffix, 'cope1', 
    'cope_ps_lh_all_ROI10.csv', atlas_filename, atlas_names, mni, NJOBS)

pf.get_percent_signal(WD, subjects, task, suffix, 'cope2', 
    'cope_ps_rh_all_ROI10.csv', atlas_filename, atlas_names, mni, NJOBS)

pf.get_percent_signal(WD, subjects, task, suffix, 'cope3', 
    'cope_ps_lh_rh_all_ROI10.csv', atlas_filename, atlas_names, mni, NJOBS)

pf.get_percent_signal(WD, subjects, task, suffix, 'cope4', 
    'cope_ps_rh_lh_all_ROI10.csv', atlas_filename, atlas_names, mni, NJOBS)


stophere


pf.get_activations(WD, subjects, task, suffix, 'zstat1.nii.gz', 
    'zstat_lh_all.txt', 'zstat_lh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'zstat2.nii.gz', 
    'zstat_rh_all.txt', 'zstat_rh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'zstat3.nii.gz', 
    'zstat_lh-rh_all.txt', 'zstat_lh-rh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'zstat4.nii.gz', 
    'zstat_rh-lh_all.txt', 'zstat_rh-lh_all.nii.gz')


pf.get_activations(WD, subjects, task, suffix, 'cope1.nii.gz', 
    'cope_lh_all.txt', 'cope_lh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'cope2.nii.gz', 
    'cope_rh_all.txt', 'cope_rh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'cope3.nii.gz', 
    'cope_lh-rh_all.txt', 'cope_lh-rh_all.nii.gz')

pf.get_activations(WD, subjects, task, suffix, 'cope4.nii.gz', 
    'cope_rh-lh_all.txt', 'cope_rh-lh_all.nii.gz')

#pf.get_ICA_lists(WD, subjects, rest, suffix, NVOLS)



# compute matrices
pf.get_fc_matrices(WD, subjects, rest, 
    atlas_filename, TR, label, FWHM, LOW_PASS, HIGH_PASS, 'lw', 
    func_mni_filename, LW_FC_FILENAME, NJOBS, NREMOVE, MINVOLS)

pf.get_fc_matrices(WD, subjects, rest, 
    atlas_filename, TR, label, FWHM, LOW_PASS, HIGH_PASS, 'full', 
    func_mni_filename, FULL_FC_FILENAME, NJOBS, NREMOVE, MINVOLS)

if task != rest:
    MINVOLS_TASK = parameters['MINVOLS_TASK']
    LOW_PASS_TASK = parameters['LOW_PASS_TASK']
    label_task = parameters['LABEL_TASK']
    FULL_FC_FILENAME = parameters['FULL_FC_FILENAME_TASK']
    LW_FC_FILENAME = parameters['LW_FC_FILENAME_TASK']

    pf.get_fc_matrices(WD, subjects, task, 
        atlas_filename, TR, label_task, FWHM, LOW_PASS_TASK, HIGH_PASS, 'lw', 
        func_mni_filename, LW_FC_FILENAME, NJOBS, NREMOVE, MINVOLS_TASK)

    pf.get_fc_matrices(WD, subjects, task, 
        atlas_filename, TR, label_task, FWHM, LOW_PASS_TASK, HIGH_PASS, 'full', 
        func_mni_filename, FULL_FC_FILENAME, NJOBS, NREMOVE, MINVOLS_TASK)




