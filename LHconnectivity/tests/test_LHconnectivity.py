
import sys, os, json
sys.path.append('..')
import LHconnectivity.preprocess_funcs as pf
import LHconnectivity as lh

PARAMS_FILE = os.path.join(lh.__path__[0], 'data', 'parameters_Stockholm.json')
with open(PARAMS_FILE) as data_file:    
    parameters = json.load(data_file)

WD = parameters['WD']
subject = parameters['TEST_SUBJECT']
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

task = 'task-localizer'
suffix = 'lr'
fsf_file = parameters['FSF_FILE']
atlas_filename = parameters['ATLAS_FILENAME']
label = 'all'

#def test_average_structurals():
#    pf.average_structurals(WD, subject, CROP_SIZE)

#def test_register_template():    
#    pf.register_template(os.path.join(WD, subject), 
#        'template.nii.gz', mni)

#def test_iterate_sessions_runs():
#   print(pf.iterate_sessions_runs(WD, subject, task))

#def test_motion_correction():
#   pf.motion_correction(WD, subject, task, TR)

#def test_register_funcs():
#  pf.register_funcs(WD, subject, task, TR, mni)

def test_run_aroma():
   pf.run_aroma(WD, subject, task, TR, AROMA_PATH, AGGR_TYPE, mni)

def test_get_fc_matrices():
   pf.get_fc_matrices(WD, subject, task, atlas_filename, TR, label, FWHM, LOW_PASS, 
   HIGH_PASS, 'lw', LW_FC_FILENAME)
   pf.get_fc_matrices(WD, subject, task, atlas_filename, TR, label, FWHM, LOW_PASS, 
   HIGH_PASS, 'full', FULL_FC_FILENAME)

def test_process_fMRI():
   pf.process_fMRI(WD, subject, task, suffix, fsf_file)

#test_average_structurals()
#test_iterate_sessions_runs()
#test_register_template()
#test_motion_correction()
#test_register_funcs()
#test_run_aroma()
#test_process_fMRI()
#test_get_fc_matrices()
