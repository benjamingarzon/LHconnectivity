from sklearn.covariance import ledoit_wolf
import os, shutil
from nilearn.input_data import NiftiLabelsMasker
from scipy import linalg
from scipy.spatial.distance import squareform, pdist
import numpy as np
import pandas as pd

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


def get_fc(WD, func_mni_filename, TR, label, atlas_filename, confounds_filename, FWHM, LOW_PASS, 
    HIGH_PASS, TYPE, events_dir):
    """
    Extract connectivity matrix given atlas and processed fMRI data.
    """

    os.chdir(WD) 

    masker = NiftiLabelsMasker(labels_img=atlas_filename, 
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

    if os.path.exists(os.path.join(events_dir, 'EV1.csv')):

        if label != 'all':
            conditions = find_labels(events_dir, TR, nvols)
            condition_mask = conditions['label'] == label
            X = X[condition_mask]
        if TYPE == 'lw':
            fc_mat, _ = ledoit_wolf(X)
            np.fill_diagonal(fc_mat, 0)
            fc = squareform(fc_mat)

        if TYPE == 'full':
            fc = pdist(X.T, metric = 'correlation')
        shutil.rmtree('nilearn_cache')
    else: 
        fc = ()
    return(fc)


