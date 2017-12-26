
# assume that data are registered in MNI
import sys, glob, itertools, os
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import json

sys.path.append('/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/LHconnectivity/')
sys.path.append('/home/ALDRECENTRUM/benjamin.garzon/Software/LeftHand/RSA/')
import LHconnectivity as lh
from LHconnectivity import preprocess_funcs as pf
from decoding_funcs import decode_inter_session_SVM
from nipype.interfaces import fsl

label1 = sys.argv[1]
label2 = sys.argv[2]
mask_filename = sys.argv[3]
out_file = sys.argv[4]
PARAMS_FILE = os.path.join(lh.__path__[0], 'data', sys.argv[5])
with open(PARAMS_FILE) as data_file:    
    parameters = json.load(data_file)

WD = parameters['WD']
NJOBS = parameters['NJOBS']
TR = parameters['TR']
task = parameters['TASK']
mni = (parameters['MNI_FILENAME'], parameters['MNI_BRAIN_FILENAME'])


def register(WD, WD_template, mni):

    os.chdir(WD)

    if not os.path.exists('fMRI_mni.nii.gz'):
        mni_filename = mni[0]   
        aw = fsl.ApplyWarp() 
        aw.run(in_file = 'fMRI_mcf.nii.gz',
            ref_file = mni_filename,
            field_file = os.path.join(WD_template, 'warps.nii.gz'), 
            out_file = 'fMRI_mni.nii.gz', 
            premat = 'func2struct.mat')

def mytask(WD, subject, out_file, task, mni):

    results = []

    os.chdir(WD)

    sessions, runs, paths = pf.iterate_sessions_runs(WD, subject, task)    

    WD_template = os.path.join(WD, subject, 'average/anat')
    items = zip(sessions, runs)

    for items1, items2 in \
        itertools.product(items, items): 
       
        session1, run1 = items1
        session2, run2 = items2
        os.chdir(os.path.join(WD, subject))
        WD1 = os.path.abspath(os.path.join(session1, 'funcproc_' + task, run1))
        WD2 = os.path.abspath(os.path.join(session2, 'funcproc_' + task, run2))
        events_dir1 = os.path.join(WD, subject, session1, 'func', subject + 
            '_' + task +  '_' + run1 + '_events') 
        events_dir2 = os.path.join(WD, subject, session2, 'func', subject + 
            '_' + task +  '_' + run2 + '_events')         

        register(WD1, WD_template, mni)        
        register(WD2, WD_template, mni)        

#                        session1 in ['1', '2', '3', '4', '5'] and \
        print(WD1, WD2)
        accuracy = decode_inter_session_SVM(WD1, WD2, events_dir1, events_dir2, 
            TR, label1, label2, mask_filename)
        fd1 = np.mean(pd.read_csv(os.path.join(WD1, 'fd.txt')).values)
        fd2 = np.mean(pd.read_csv(os.path.join(WD2, 'fd.txt')).values)
        results.append((subject, session1, run1, subject, session2, run2, accuracy, fd1, fd2))
        f = open(out_file, 'a')
        f.write('%s,%s,%s,%s,%s,%s,%f,%f,%f'%(subject.split('-')[1], 
            session1.split('-')[1], run1.split('-')[1], subject.split('-')[1], 
            session2.split('-')[1], run2.split('-')[1], accuracy, fd1, fd2))
        f.write('\n')
        f.close()

    return(results)                
    
os.chdir(WD)
#subjects = glob.glob('sub-LH100[1-5]')
subjects = glob.glob('sub-*')
f = open(out_file, 'w')
f.write('SUBJECT_TRAINING,SESSION_TRAINING,RUN_TRAINING,SUBJECT_TESTING,SESSION_TESTING,RUN_TESTING,ACCURACY,FD1,FD2\n')
f.close()

results = Parallel(n_jobs = NJOBS)(delayed(mytask)(WD, subject, out_file, task, mni) 
    for subject in subjects) 

#reorder and rewrite results
results = [item for sublist in results for item in sublist]
results = pd.DataFrame(results, columns=['SUBJECT_TRAINING', 'SESSION_TRAINING',
    'RUN_TRAINING','SUBJECT_TESTING','SESSION_TESTING','RUN_TESTING','ACCURACY','FD1','FD2'])
results.sort_values(by = ['SUBJECT_TRAINING', 'SESSION_TRAINING',
    'RUN_TRAINING','SUBJECT_TESTING','SESSION_TESTING','RUN_TESTING'], inplace = True)
results.to_csv(out_file, index = False)


