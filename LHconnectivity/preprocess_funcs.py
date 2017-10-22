from __future__ import division
import numpy as np
from nilearn.input_data import NiftiMasker
from nipype.interfaces import fsl
from nilearn import image
import os, glob, sys, random, itertools, shutil, errno
#from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform, pdist
from scipy.stats import pearsonr
from LHconnectivity import get_fc
import pandas as pd


def average_structurals(WD, subject, CROP_SIZE):
    """
    Downsample, crop and obtain an average template for all the structurals. 
    """
    
    os.chdir(os.path.join(WD, subject))

    if not os.path.exists('average/anat'):
        os.makedirs('average/anat')

    structural_files = glob.glob('sess*/anat/*T1w.nii.gz')
    for structural_file in structural_files: 
        resampled_structural_file = structural_file.split('.')[0] + \
            '_1mm.nii.gz'
        flt = fsl.FLIRT()
        flt.run(in_file = structural_file, 
                reference = structural_file,
                interp = 'spline', 
                apply_isoxfm = 1,
                out_file = resampled_structural_file,
                out_matrix_file = 'tmp.mat')
        os.remove('tmp.mat')
        os.system('fslroi ' + resampled_structural_file + ' ' +
            resampled_structural_file + ' ' + CROP_SIZE)
    
    os.system('mri_robust_template --mov sess*/anat/*T1w*_1mm.nii.gz ' + 
        '--template average/anat/template.nii.gz --satit')
    
    # clean up
    os.system('rm sess*/anat/*T1w*_1mm.nii.gz ')    

def motion_correction(WD, subject, task, TR):

    """
    Do motion correction and get parameters. 
    """

    def do_motion_correction(WD, TR):
        os.chdir(WD)
        if not os.path.exists('fMRI_mcf.nii.gz'):
            os.system('slicetimer -i fMRI.nii.gz -o fMRI_st.nii.gz --odd -r ' + 
                str(TR))
            os.system('mcflirt -in fMRI_st.nii.gz -o fMRI_mcf -plots')
            os.system('rm fMRI_st.nii.gz')

        if not os.path.exists('mean.nii.gz'):
            os.system('fslmaths fMRI_mcf.nii.gz -Tmean mean.nii.gz')

        if not os.path.exists('brain_mask.nii.gz'):
            btr = fsl.BET()
            btr.run(in_file = 'mean.nii.gz',
                    robust = True,
                    out_file = 'brain.nii.gz',
                    mask = True, 
                    frac = 0.3, 
                    no_output = True)

        if not os.path.exists('fd.txt'):
            os.system('fsl_motion_outliers -i fMRI.nii.gz ' + 
                '-m brain_mask.nii.gz -o confounds.txt -s fd.txt --fd')

        # Calculate Friston 24 parameters
        if not os.path.exists('fMRI_mcf_fr24.par'):
            pars = pd.read_csv('fMRI_mcf.par', delim_whitespace = True).values
            pars_diff = np.diff(pars, axis = 0)
            pars_diff = np.concatenate((np.zeros(( 1, pars.shape[1])), 
                pars_diff), axis = 0)
            fr24 = np.concatenate((pars, pars_diff, pars**2, pars_diff**2), 
                axis = 1)
            np.savetxt('fMRI_mcf_fr24.par', fr24, delimiter = ' ')

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    
    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))

        proc_dir = os.path.join(session, 'funcproc_' + task, run)
        if not os.path.exists(proc_dir):
            os.makedirs(proc_dir)            
        link = os.path.join(proc_dir, 'fMRI.nii.gz')
        if os.path.exists(link):
            os.remove(link)        
        os.symlink(os.path.abspath(path), link) 
        do_motion_correction(proc_dir, TR)


def register_template(WD, structural_filename, mni):
    
    os.chdir(os.path.join(WD, 'average/anat'))
    
    mni_filename = mni[0]
    mni_brain_filename = mni[1]
    

    # segment and make nice
    if not os.path.exists('structural_wmseg.nii.gz'):

        shutil.copyfile(structural_filename, 'structural.nii.gz')
        os.system('fsl_anat --clobber --weakbias --nosubcortseg ' + 
            '-i structural.nii.gz')
        shutil.copyfile('structural.anat/T1_biascorr.nii.gz', 
            'structural.nii.gz')
        shutil.copyfile('structural.anat/T1_biascorr_brain.nii.gz', 
            'structural_brain.nii.gz')
        shutil.copyfile('structural.anat/T1_biascorr_brain_mask.nii.gz', 
            'structural_brain_mask.nii.gz')
        os.system('fslmaths structural.anat/T1_fast_pve_2.nii.gz -thr 0.5 ' + 
            '-bin structural_wmseg.nii.gz')
        shutil.rmtree('structural.anat')

    # register structural
    if not os.path.exists('warps.nii.gz'):
        flt = fsl.FLIRT()
        flt.run(in_file = 'structural_brain.nii.gz', 
                reference = mni_brain_filename,
                dof = 12, 
                out_matrix_file = 'structural.mat')

        fnt = fsl.FNIRT()
        fnt.run(in_file = 'structural.nii.gz', 
                ref_file = mni_filename, 
                affine_file = 'structural.mat', 
                inmask_file = 'structural_brain_mask.nii.gz', 
                field_file = 'warps.nii.gz')


def register_funcs(WD, subject, task, TR, mni):

    """
    Compute and apply registrations to fMRI data. 
    """

    mni_filename = mni[0]

    def do_register_funcs(WD, WD_template):
        os.chdir(WD)

        # register fMRI
        if not os.path.exists('func2struct.mat'):
            epr = fsl.EpiReg()
            epr.run(epi = 'mean.nii.gz', 
                    t1_head = os.path.join(WD_template, 'structural.nii.gz'),
                    t1_brain = os.path.join(WD_template, 
                        'structural_brain.nii.gz'),
                    out_base = 'epireg',
                    wmseg = os.path.join(WD_template, 
                        'structural_wmseg.nii.gz'))
            shutil.copyfile('epireg.mat', 'func2struct.mat')
            os.system('rm epireg*')

    
        if not os.path.exists('fMRI_mni.nii.gz'):
            aw = fsl.ApplyWarp() 
            aw.run(in_file = 'fMRI_mcf.nii.gz',
                   ref_file = mni_filename,
                   field_file = os.path.join(WD_template, 'warps.nii.gz'), 
                   out_file = 'fMRI_mni.nii.gz', 
                   premat = 'func2struct.mat')

            #get average volume to check that everything is alrighty
            os.system('fslmaths ' + 'fMRI_mni.nii.gz' + 
                ' -Tmean fMRI_mni_mean.nii.gz')

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

    WD_template = os.path.join(WD, subject, 'average/anat')

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))

        proc_dir = os.path.join(session, 'funcproc_' + task, run)
        do_register_funcs(proc_dir, WD_template)


def iterate_sessions_runs(WD, subject, task):
   
    """
    Return a list with scans for all sessions and runs of the subject.
    """

    paths_list = []
    runs_list = []
    sessions_list = []

    os.chdir(os.path.join(WD, subject))

    sessions = glob.glob('sess-*')
    for session in sessions:
        paths = glob.glob(os.path.join(session, 'func', subject + '_' +
        task + '_run*.nii.gz'))
        paths_list += paths
        runs = [os.path.basename(path).split('.')[0].split('_')[2] 
            for path in paths]
        runs_list += runs
        sessions_list += [session]*len(runs)

    return((sessions_list, runs_list, paths_list))


def run_aroma(WD, subject, task, TR, AROMA_PATH, AGGR_TYPE, mni):

    """
    Run aroma for denoising. 
    """

    mni_filename = mni[0]

    def do_run_aroma(WD, WD_template):
        os.chdir(WD)
        command = 'python ' + AROMA_PATH + \
            ' -i ' + os.path.abspath('fMRI_mcf.nii.gz') + \
            ' -o ' + os.path.join(WD, 'aroma') + \
            ' -a ' + os.path.abspath('func2struct.mat') + \
            ' -w ' + os.path.join(WD_template, 'warps.nii.gz') + \
            ' -mc ' + os.path.abspath('fMRI_mcf.par') + \
            ' -tr ' + str(TR) + \
            ' -m ' + os.path.abspath('brain_mask.nii.gz') + \
            ' -den both'
        os.system(command)

        link = os.path.join(WD, 'fMRI_denoised.nii.gz')
        if os.path.exists(link):
            os.remove(link)        
        fmri_filename = os.path.join(WD, 'aroma', 'denoised_func_data_' + 
            AGGR_TYPE + '.nii.gz')
        os.symlink(fmri_filename, link) 

        aw = fsl.ApplyWarp() 
        aw.run(in_file = fmri_filename,
            ref_file = mni_filename,
            field_file = os.path.join(WD_template, 'warps.nii.gz'), 
            out_file = 'fMRI_denoised_mni.nii.gz', 
            premat = 'func2struct.mat')
        
#        aroma = fsl.ICA_AROMA()
#        aroma.inputs.in_file = 'fMRI_mcf.nii.gz'
#        aroma.inputs.mat_file = 'func2struct.mat'
#        aroma.inputs.fnirt_warp_file = os.path.join(WD_template, 'warps.nii.gz')
#        aroma.inputs.motion_parameters = 'fMRI_mcf.par'
#        aroma.inputs.mask = 'brain_mask.nii.gz'
#        aroma.inputs.denoise_type = 'both'
#        aroma.inputs.out_dir = 'ICA_testout'

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    WD_template = os.path.join(WD, subject, 'average/anat')

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.abspath(os.path.join(session, 'funcproc_' + task, 
            run))        
        do_run_aroma(proc_dir, WD_template)


def process_fMRI(WD, subject, task, suffix, fsf_file):
    """
    Do the processing for fMRI data.
    """

    def do_process_fMRI(WD, WD_template, fsf_file):

        # loop over EVs instead
        os.chdir(WD)        
        shutil.copyfile(fsf_file, 'fMRI.fsf')
        os.system("sed -i 's @fMRI %s ' %s"%(os.path.abspath('fMRI.nii.gz'),
            'fMRI.fsf'))
        os.system("sed -i 's @STRUCTURAL_BRAIN %s ' %s"%
            (os.path.join(WD_template, 'structural_brain.nii.gz'), 'fMRI.fsf'))
        os.system("sed -i 's @EV1 %s %s"%
            (os.path.abspath('EV1.csv'), 'fMRI.fsf'))
        os.system("sed -i 's @EV2 %s %s"%
            (os.path.abspath('EV2.csv'), 'fMRI.fsf'))
        os.system("sed -i 's @ANALYSIS %s %s"%
            (os.path.join(os.getcwd(), 'analysis'), 'fMRI.fsf'))

        os.system('feat fMRI.fsf')

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)    
    WD_template = os.path.join(WD, subject, 'average/anat')

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.join(session, 'funcproc_' + task + '-' + suffix, run)
        if not os.path.exists(proc_dir):
            os.makedirs(proc_dir)            
        link = os.path.join(proc_dir, 'fMRI.nii.gz')
        if os.path.exists(link):
            os.remove(link)        
        os.symlink(os.path.abspath(path), link) 
        do_process_fMRI(proc_dir, WD_template, fsf_file)


def get_fc_matrices(WD, subject, task, atlas_filename, TR, label, FWHM, LOW_PASS, 
    HIGH_PASS, TYPE, output_filename):
    """
    Compute FC matrices for each run. 
    """
    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    results = []
    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))

        proc_dir = os.path.join(session, 'funcproc_' + task, run)
        confounds_filename = 'fMRI_mcf_fr24.par'
        fc = get_fc(proc_dir, TR, label, atlas_filename, confounds_filename, 
        FWHM, LOW_PASS, HIGH_PASS, TYPE)
        fd = np.mean(pd.read_csv(os.path.join(proc_dir, 'fd.txt')).values)

        results.append((subject, int(session), int(run), fd) + tuple(fc))
  
    results = pd.DataFrame(results, columns=['SUBJECT','SESSION','RUN','FD'] + 
         [ 'fc_%03d'%(i+1)  for i in range(len(fc)) ])
    results.sort_values(by = ['SUBJECT', 'SESSION', 'RUN'], inplace = True)
    results.to_csv(output_filename, index = False)


