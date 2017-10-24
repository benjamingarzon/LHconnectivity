import sys, glob, os
import numpy as np
from joblib import Parallel, delayed
import json

sys.path.append('..')
import LHconnectivity.preprocess_funcs as pf
import LHconnectivity as lh

PARAMS_FILE = os.path.join(lh.__path__[0], 'data', 'parameters_Stockholm.json')
with open(PARAMS_FILE) as data_file:    
    parameters = json.load(data_file)

WD = parameters['WD']
NJOBS = 5#parameters['NJOBS']
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
task = 'task-localizer'
suffix = 'lr'
fsf_file = parameters['FSF_FILE']
atlas_filename = parameters['ATLAS_FILENAME']
label = 'all'

os.chdir(WD)
subjects = glob.glob('sub-LH100[1-5]')

# get average template
Parallel(n_jobs = NJOBS)(delayed(pf.average_structurals)(WD, subject, CROP_SIZE) 
    for subject in subjects) 

# register template
Parallel(n_jobs = NJOBS)(delayed(pf.register_template)
    (os.path.join(WD, subject), 'template.nii.gz', mni) 
    for subject in subjects) 

# motion_correction
Parallel(n_jobs = NJOBS)(delayed(pf.motion_correction)(WD, subject, task, TR) 
    for subject in subjects) 

# register functional data
Parallel(n_jobs = NJOBS)(delayed(pf.register_funcs)(WD, subject, task, TR, mni) 
    for subject in subjects) 

# run AROMA
Parallel(n_jobs = NJOBS)(delayed(pf.run_aroma)(WD, subject, task, TR, 
    AROMA_PATH, AGGR_TYPE, mni) 
    for subject in subjects) 

# compute matrices
pf.get_fc_matrices(WD, subject, task, 
    atlas_filename, TR, label, FWHM, LOW_PASS, HIGH_PASS, 'lw', 
    func_mni_filename, LW_FC_FILENAME, NJOBS)

pf.get_fc_matrices(WD, subject, task, 
    atlas_filename, TR, label, FWHM, LOW_PASS, HIGH_PASS, 'full', 
    func_mni_filename, FULL_FC_FILENAME, NJOBS)

# process fMRI task data
Parallel(n_jobs = NJOBS)(delayed(pf.process_fMRI)(WD, subject, task, suffix, 
    fsf_file)
    for subject in subjects) 


