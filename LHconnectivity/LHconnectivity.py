from sklearn.covariance import ledoit_wolf
import os, shutil
from nilearn.input_data import NiftiLabelsMasker, NiftiMasker
from scipy import linalg
from scipy.spatial.distance import squareform, pdist
import numpy as np
import pandas as pd
import warnings

def find_labels(WD, TR, nvols):
    """
    Get labels from EV files.
    """

    left_onsets = pd.read_csv(os.path.join(WD, 'EV1.csv'), sep = ' ', 
        names = ['onset', 'duration', 'value' ] )
    right_onsets = pd.read_csv(os.path.join(WD, 'EV2.csv'), sep = ' ', 
        names = ['onset', 'duration', 'value' ] )

    times = np.arange(nvols)*TR
    labels = nvols*['rest']
    blocks = nvols*[0]
    thisblock = 0

    for (i, time) in enumerate(times):
        if any( (time >= np.array(left_onsets['onset'])) & 
            (time < np.array(left_onsets['onset'] + 
            left_onsets['duration']))):
            labels[i] = 'left'
        if any( (time >= np.array(right_onsets['onset'])) & 
            (time < np.array(right_onsets['onset'] + 
            right_onsets['duration']))):
            labels[i] = 'right'

        if (i > 0):
           if labels[i] != labels[i - 1]:
               thisblock = thisblock + 1
        blocks[i] = thisblock

    # restrict the analysis to specific conditions (left / right / ...)
    conditions = np.rec.array(zip(labels, blocks), [('label', 'S12'), 
        ('block', 'S12')])

    return(conditions)

def compute_fc(X, TYPE):
    """
    Compute different types of FC estimators. 
    """
    if TYPE == 'lw':
        fc_mat, _ = ledoit_wolf(X)
        np.fill_diagonal(fc_mat, 0)
        fc = squareform(fc_mat)

    if TYPE == 'full':
        fc = 1 - pdist(X.T, metric = 'correlation')
#        fc = squareform(np.corrcoef(X.T))
#        fc = np.arctanh(fc)
    
    return(fc) 

def get_fc(WD, func_mni_filename, TR, label, atlas_filename, confounds_filename, 
    FWHM, LOW_PASS, HIGH_PASS, TYPE, events_dir, NREMOVE, MINVOLS):
    """
    Extract connectivity matrix given atlas and processed fMRI data.
    """

    os.chdir(WD) 

    LOW_PASS = None if LOW_PASS == 0 else LOW_PASS
    HIGH_PASS = None if HIGH_PASS == 0 else HIGH_PASS

    masker = NiftiLabelsMasker(labels_img = atlas_filename, 
                     smoothing_fwhm = FWHM,
                     standardize = True, 
                     t_r = TR, 
                     detrend = True,
                     low_pass = LOW_PASS, 
                     high_pass = HIGH_PASS, 
                     memory_level = 1, 
                     memory = 'nilearn_cache')

    X = masker.fit_transform(func_mni_filename, confounds = confounds_filename)
    nvols = X.shape[0]

    if label == 'all':
        fc = compute_fc(X, TYPE)

    else:
        if os.path.exists(os.path.join(events_dir, 'EV1.csv')):
            conditions = find_labels(events_dir, TR, nvols)
            condition_mask = conditions['label'] == label            
            mask_diff = (np.diff(condition_mask*1)==1)
            # remove volumes 
            for x in np.where(mask_diff)[0]:
                condition_mask[x: np.min((x + NREMOVE + 1, 
                    len(condition_mask)))] = False
            
            if np.sum(condition_mask) > MINVOLS: 
                X = X[condition_mask]
                fc = compute_fc(X, TYPE) 
            else:
                print(os.path.abspath('.'))
                print('Too few vols! %d'%(np.sum(condition_mask))) 
                fc = ()
                print(X.shape)
        else: 
            fc = ()

    if os.path.exists('nilearn_cache'):
        shutil.rmtree('nilearn_cache')    

    return(fc)

def get_ALFF(WD, func_mni_filename, TR, mask_filename, FWHM, LOW_PASS, 
    HIGH_PASS, confounds_filename):
    """
    Compute ALFF.
    """

    os.chdir(WD) 

    LOW_PASS = None if LOW_PASS == 0 else LOW_PASS
    HIGH_PASS = None if HIGH_PASS == 0 else HIGH_PASS

    masker = NiftiMasker(mask_img = mask_filename, 
                     smoothing_fwhm = FWHM,
                     standardize = False, 
                     t_r = TR, 
                     detrend = True,
                     low_pass = LOW_PASS, 
                     high_pass = HIGH_PASS, 
                     memory_level = 1, 
                     memory = 'nilearn_cache')

    X = masker.fit_transform(func_mni_filename, confounds = confounds_filename)

    ALFF = np.std(X, axis = 0)
    m = np.mean(ALFF)
    sd = np.std(ALFF)

    # normalize?
    ALFF = (ALFF - m )/sd

    ALFF_img = masker.inverse_transform(ALFF)
#    ALFF_img.to_filename('ALFF.nii.gz')

    if os.path.exists('nilearn_cache'):
        shutil.rmtree('nilearn_cache')   

    return(ALFF_img) 




