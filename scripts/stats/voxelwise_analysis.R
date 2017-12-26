
rm(list = ls())
library(oro.nifti)
library(doParallel)
library(foreach)

setwd('~/Software/LeftHand/LHconnectivity/scripts/stats')
source('./analysis_tests.R')
source('./fc_analysis_func.R')
source('~/Software/LeftHand/VBM/stats_aux.R')
source('~/Software/LeftHand/VBM/get_session_data.R')

data.orig = data

NCLUSTERS = 15
REMOVE_OUTLIERS = F

data.imaging = read.table('/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh_all.txt', header=T, sep = ',')
data.imaging$SUBJECT = substr(data.imaging$SUBJECT, 5, 11)

data = inner_join(data.imaging, data.orig, by = c('SUBJECT','SESSION'))

OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Stockholm/activations'
#MASK_FILE = '/home/share/LeftHand/fMRI/RSA/motor_mask.nii.gz'
MASK_FILE = '/home/ALDRECENTRUM//benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/GMmask_MNI_2mm_Stockholm_masked.nii.gz'

#ACTIVATIONS

do_tests = do_tests_LME_TRAINING

OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Stockholm/activations/lh-rh'
IMAGING_FILE = '/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh-rh_all.nii.gz'
vbanalysis(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=F)

OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Stockholm/activations/rh-lh'
IMAGING_FILE = '/home//share/LeftHand/LHconnectivity/Stockholm/zstat_rh-lh_all.nii.gz'
vbanalysis(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=F)


OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Stockholm/activations/lh'
IMAGING_FILE = '/home//share/LeftHand/LHconnectivity/Stockholm/zstat_lh_all.nii.gz'
vbanalysis(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=F)

OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Stockholm/activations/rh'
IMAGING_FILE = '/home//share/LeftHand/LHconnectivity/Stockholm/zstat_rh_all.nii.gz'
vbanalysis(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=F)

MASK_FILE = '/home/ALDRECENTRUM//benjamin.garzon/Software/LeftHand/LHconnectivity/LHconnectivity/data/GMmask_MNI_2mm_Berlin_masked.nii.gz'
#MASK_FILE = '/home/share/LeftHand/fMRI/RSA/motor_mask.nii.gz'

#ALFF
data.imaging = read.table('/home//share/LeftHand/LHconnectivity/Berlin/ALFF_resting.txt', header=T, sep = ',')
data.imaging$SUBJECT = substr(data.imaging$SUBJECT, 5, 8)
data.imaging$GROUP = factor('none', levels = c('none','experimental','control'))
data.imaging$GROUP[data.imaging$SUBJECT %in% Berlin_experimental] = 'experimental'
data.imaging$GROUP[data.imaging$SUBJECT %in% Berlin_control] = 'control'

do_tests = do_tests_LME_SESSION_Berlin
data = data.imaging

OUTPUT_DIR = '/home//share/LeftHand/LHconnectivity/Berlin/ALFF'
IMAGING_FILE = '/home//share/LeftHand/LHconnectivity/Berlin/ALFF_resting.nii.gz'
#vbanalysis(IMAGING_FILE, OUTPUT_DIR, data, MASK_FILE, do_tests, excluded=NULL, to_gifti='', remove_outliers=F)

