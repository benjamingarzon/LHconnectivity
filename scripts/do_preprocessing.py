import sys, glob, os
import numpy as np
from joblib import Parallel, delayed
import json

sys.path.append('..')
from LHconnectivity.preprocess_funcs import average_structurals, register

PARAMS_FILE = '../LHconnectivity/data/parameters_Stockholm.json'
with open(PARAMS_FILE) as data_file:    
    parameters = json.load(data_file)

NJOBS = parameters['NJOBS']
WD = parameters['WD']
subject = parameters['TEST_SUBJECT']
CROP_SIZE = parameters['CROP_SIZE']
TR = parameters['TR']
mni = (parameters['MNI_FILENAME'], parameters['MNI_BRAIN_FILENAME'])

os.chdir(WD)
subjects = glob.glob('sub-*')

# get average template
Parallel(n_jobs = NJOBS)(delayed(average_structurals)(WD, subject, CROP_SIZE)) 
    for subject in subjects) 

# register template
Parallel(n_jobs = NJOBS)(delayed(register_template))
    (os.path.join(WD, subject), 'template.nii.gz', mni) 
    for subject in subjects) 

# motion_correction
Parallel(n_jobs = NJOBS)(delayed(motion_correction))(WD, subject, TR) 
    for subject in subjects) 

# register subject
Parallel(n_jobs = NJOBS)(delayed(register_subject))(WD, subject, TR) 
    for subject in subjects) 

# run AROMA




# get connectivity matrices

label1 = sys.argv[1]
label2 = sys.argv[2]
mask_filename = sys.argv[3]
out_file = sys.argv[4]
weight_filename = sys.argv[5]

def mytask(scans, mysubject):
    results = []
 
    for (s, subject) in enumerate(scans["SUBJECT"]):
        if subject == mysubject:
            session = str(scans.SESSION[s])
            name = subject + '.' +  session
          
#            WD_subject = DataDir + subject + '/' 
            WDs = glob.glob(DataDir + name + '/run?/')
            for WD in WDs:
                run =  WD.split("/")[7][3]
        
                if os.path.exists(WD + 'EV1.csv') and os.path.exists(WD + 
                    'fMRI.nii.gz'):
                    print(WD)
                    structural_filename = glob.glob("%s/avg_%s.*_3T_T1w_MPR1.nii.gz"%
                        (VBM_DIR, subject))[0]
                    



results = Parallel(n_jobs = NJOBS)(delayed(average_structurals)(scans, subject, out_file) 
    for subject in subjects) 
