from __future__ import division
import numpy as np
from nilearn.input_data import NiftiMasker
from nipype.interfaces import fsl
from nipype.algorithms.confounds import CompCor
from nilearn import image
import os, glob, sys, random, itertools, shutil, errno
#from scipy.spatial import distance_matrix
from scipy.spatial.distance import squareform, pdist
from scipy.stats import pearsonr
from LHconnectivity import get_fc, get_ALFF
import pandas as pd
from joblib import Parallel, delayed


def average_structurals(WD, subject, CROP_SIZE):
    """
    Downsample, crop and obtain an average template for all the structurals. 
    """
    
    os.chdir(os.path.join(WD, subject))

    if not os.path.exists('average/anat'):
        os.makedirs('average/anat')

    if not os.path.exists('average/anat/template.nii.gz'):

        structural_files = glob.glob('sess*/anat/*T1w.nii.gz')
        for structural_file in structural_files: 
            resampled_structural_file = structural_file.split('.')[0] + \
                '_1mm.nii.gz'
            flt = fsl.FLIRT()
            flt.run(in_file = structural_file, 
                    reference = structural_file,
                    interp = 'trilinear', 
                    apply_isoxfm = 1,
                    out_file = resampled_structural_file,
                    out_matrix_file = 'tmp.mat')
            os.remove('tmp.mat')
            os.system('fslroi ' + resampled_structural_file + ' ' +
                resampled_structural_file + ' ' + CROP_SIZE)
            os.system('fslmaths ' + resampled_structural_file + ' -thr 0 ' +
                resampled_structural_file)
    
        os.system('mri_robust_template --mov sess*/anat/*T1w*_1mm.nii.gz ' + 
            '--template average/anat/template.nii.gz --satit --iscale')
    
        # clean up
        os.system('rm sess*/anat/*T1w*_1mm.nii.gz')    

def motion_correction(WD, subject, task, TR, acq_order):
    """
    Do motion correction and get parameters. 
    """

    if acq_order == 'interleaved':
        acq_order_flag = ' --odd'
    if acq_order == 'top-down':
        acq_order_flag = ' --down'
    if acq_order == 'bottom-up':
        acq_order_flag = ''

    def do_motion_correction(WD, TR):
        os.chdir(WD)
        if not os.path.exists('fMRI_mcf.nii.gz'):
            if acq_order != 'none':
                os.system('slicetimer -i fMRI.nii.gz -o fMRI_st.nii.gz -r ' + 
                    str(TR) + acq_order_flag)
                os.system('mcflirt -in fMRI_st.nii.gz -o fMRI_mcf -plots')
                os.system('rm fMRI_st.nii.gz')
            else:
                os.system('mcflirt -in fMRI.nii.gz -o fMRI_mcf -plots')
        

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
            pars = pd.read_csv('fMRI_mcf.par', delim_whitespace = True, 
                header = None).values
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
        if not os.path.exists('fMRI_denoised_mni.nii.gz'):

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
            shutil.rmtree('aroma')
            os.unlink(link)
        
    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    WD_template = os.path.join(WD, subject, 'average/anat')

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.abspath(os.path.join(session, 'funcproc_' + task, 
            run))        
        do_run_aroma(proc_dir, WD_template)



def run_compcor(func_mni_filename, mask_filename):

    def do_run_compcor(WD):
        os.chdir(WD)
        ccinterface = CompCor()
        ccinterface.inputs.realigned_file = func_mni_filename
        ccinterface.inputs.mask_files = mask_filename
        ccinterface.inputs.num_components = 1
        ccinterface.inputs.pre_filter = 'polynomial'
        ccinterface.inputs.regress_poly_degree = 2
        ccinterface.run()

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

    os.chdir(os.path.join(WD, subject))

    for session, run, path in zip(sessions, runs, paths):
        proc_dir = os.path.join(session, 'funcproc_' + task, run)

        do_run_compcor(WD)


def process_fMRI(WD, subject, task, suffix, fsf_file, mni):
    """
    Do the processing for fMRI data.
    """

    def do_process_fMRI(WD, events_dir, WD_template, fsf_file, mni):

        mni_filename = mni[0]

        # loop over EVs instead
        os.chdir(WD)        
        shutil.copyfile(fsf_file, 'fMRI.fsf')
        os.system("sed -i 's @fMRI %s ' %s"%(os.path.abspath('fMRI.nii.gz'),
            'fMRI.fsf'))
        os.system("sed -i 's @STRUCTURAL_BRAIN %s ' %s"%
            (os.path.join(WD_template, 'structural_brain.nii.gz'), 'fMRI.fsf'))
        os.system("sed -i 's @EV1 %s ' %s"%
            (os.path.join(events_dir, 'EV1.csv'), 'fMRI.fsf'))
        os.system("sed -i 's @EV2 %s ' %s"%
            (os.path.join(events_dir, 'EV2.csv'), 'fMRI.fsf'))
        os.system("sed -i 's @ANALYSIS %s ' %s"%
            (os.path.join(WD, 'analysis'), 'fMRI.fsf'))
        
        if os.path.exists('analysis.feat'):
            if not os.path.exists(
            'analysis.feat/rendered_thresh_zstat1.nii.gz'):
                shutil.rmtree('analysis.feat')
                os.system('feat fMRI.fsf')
          
        else:
            os.system('feat fMRI.fsf')

        # warp the contrasts to MNI
        os.system('cp analysis.feat/stats/zstat* .')
        zstats = glob.glob('zstat*')
        maps = zstats

#        tstats = glob.glob('analysis.feat/stats/tstat*')
#        cope = glob.glob('analysis.feat/stats/cope*')
#        varcope = glob.glob('analysis.feat/stats/varcope*')
#        maps = tstats + zstats + cope + varcope


        for mymap in maps:
            aw = fsl.ApplyWarp() 
            aw.run(in_file = mymap,
                ref_file = mni_filename,
                field_file = os.path.join(WD_template, 'warps.nii.gz'), 
                out_file = mymap, 
                premat = 'analysis.feat/reg/example_func2highres.mat')

        # remove analysis dir
#        shutil.rmtree('analysis.feat')

    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)    
    WD_template = os.path.join(WD, subject, 'average/anat')

    for session, run, path in zip(sessions, runs, paths):

        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.join(WD, subject, session, 'funcproc_' + task + 
            '-' + suffix, run)  
        if not os.path.exists(proc_dir):
            os.makedirs(proc_dir)            
        link = os.path.join(proc_dir, 'fMRI.nii.gz')
        if os.path.exists(link):
            os.remove(link)        
        os.symlink(os.path.abspath(path), link) 
        events_dir = os.path.join(WD, subject, session, 'func', subject + 
            '_' + task +  '_' + run + '_events') 
        do_process_fMRI(proc_dir, events_dir, WD_template, fsf_file, mni)


def do_get_fc_matrices(WD, subject, task, atlas_filename, TR, label, FWHM, 
    LOW_PASS, HIGH_PASS, TYPE, func_mni_filename, NREMOVE, MINVOLS):
    """
    Compute FC matrices for each run. 
    """
    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    results = []

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.join(session, 'funcproc_' + task, run)
        confounds_filename = 'fMRI_mcf_fr24.par'
        events_dir = os.path.join(WD, subject, session, 'func', subject + 
            '_' + task +  '_' + run + '_events') 
        #print(os.path.join(WD, subject, proc_dir))
        fc = get_fc(proc_dir, func_mni_filename, TR, label, atlas_filename, 
            confounds_filename, FWHM, LOW_PASS, HIGH_PASS, TYPE, events_dir,
            NREMOVE, MINVOLS)
        fd = np.mean(pd.read_csv(os.path.join(WD, subject, proc_dir, 
            'fd.txt')).values)
        if len(fc) > 0:
            results.append((subject, int(session.split('-')[1]), 
                int(run.split('-')[1]), fd) + tuple(fc))

    return(results)


def get_fc_matrices(WD, subjects, task, atlas_filename, TR, label, FWHM, 
    LOW_PASS, HIGH_PASS, TYPE, func_mni_filename, output_filename, NJOBS, 
    NREMOVE, MINVOLS):
    """
    Compute FC matrices for each run for all subjects. 
    """

    results = Parallel(n_jobs = NJOBS)(delayed(do_get_fc_matrices)(WD, subject, 
    task, atlas_filename, TR, label, FWHM, LOW_PASS, HIGH_PASS, TYPE, 
    func_mni_filename, NREMOVE, MINVOLS)
    for subject in subjects) 

    results = [item for sublist in results for item in sublist]
    results = pd.DataFrame(results)
    results.columns = ['SUBJECT','SESSION','RUN','FD'] + \
        [ 'fc_%03d'%(i+1)  for i in range(results.shape[1] - 4) ]
    results.sort_values(by = ['SUBJECT', 'SESSION', 'RUN'], inplace = True)
    results.to_csv(output_filename, index = False)

def do_get_ALFF(WD, subject, task, mask_filename, TR, FWHM, 
    LOW_PASS, HIGH_PASS, func_mni_filename):
    
    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)
    results = []

    for session, run, path in zip(sessions, runs, paths):
        os.chdir(os.path.join(WD, subject))
        proc_dir = os.path.join(session, 'funcproc_' + task, run)
        confounds_filename = 'fMRI_mcf_fr24.par'
        ALFF = get_ALFF(proc_dir, func_mni_filename, TR, mask_filename, 
            FWHM, LOW_PASS, HIGH_PASS, confounds_filename)
        fd = np.mean(pd.read_csv(os.path.join(WD, subject, proc_dir, 
                'fd.txt')).values)
    if len(ALFF) > 0:
        results.append((subject, int(session.split('-')[1]), 
        int(run.split('-')[1]), fd) + tuple(ALFF))

    return(results)

def get_ALFFs(WD, subjects, task, mask_filename, TR, FWHM, 
    LOW_PASS, HIGH_PASS, func_mni_filename, output_filename, 
    output_img_filename):
    """
    Compute ALFF for each run for all subjects. 
    """
    indices = []
    maps = []
    for subject in subjects:
        os.chdir(os.path.join(WD, subject))
        sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

        for session, run, path in zip(sessions, runs, paths): 
            proc_dir = os.path.join(WD, subject, session, 'funcproc_' + task, 
            run) 

            confounds_filename = 'fMRI_mcf_fr24.par'
            ALFF = get_ALFF(proc_dir, func_mni_filename, TR, mask_filename, 
            FWHM, LOW_PASS, HIGH_PASS, confounds_filename)
            fd = np.mean(pd.read_csv(os.path.join(WD, subject, session, 
                'funcproc_' + task, run, 'fd.txt')).values)
            indices.append((subject, int(session.split('-')[1]), 
            int(run.split('-')[1]), fd) )
            maps.append(ALFF)

    os.chdir(WD)

    indices = pd.DataFrame(indices)
    indices.columns = ['SUBJECT','SESSION','RUN','FD']
    indices.to_csv(output_filename, index = False)

    # concatenate the list of images
    output_img = image.concat_imgs(maps)
    output_img.to_filename(output_img_filename)


def get_activations(WD, subjects, task, suffix, stat_filename, output_filename, 
    output_img_filename):
    """
    Join all activation maps. 
    """

    indices = []
    filenames = []
    for subject in subjects:
        os.chdir(os.path.join(WD, subject))
        sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

        for session, run, path in zip(sessions, runs, paths): 
            proc_dir = os.path.join(WD, subject, session, 'funcproc_' + task + 
                '-' + suffix, run) 

            abs_stat_filename = os.path.join(proc_dir, stat_filename)
            fd = np.mean(pd.read_csv(os.path.join(WD, subject, session, 
                'funcproc_' + task, run, 'fd.txt')).values)
            if os.path.exists(abs_stat_filename):  
                indices.append((subject, int(session.split('-')[1]), 
                    int(run.split('-')[1]), fd) )
                filenames.append(abs_stat_filename)

    os.chdir(WD)

#    indices = [item for sublist in indices for item in sublist]

    indices = pd.DataFrame(indices)
    indices.columns = ['SUBJECT','SESSION','RUN','FD']
    #indices.sort_values(by = ['SUBJECT', 'SESSION', 'RUN'], inplace = True)
    indices.to_csv(output_filename, index = False)

    # concatenate the list of images
    output_img = image.concat_imgs(filenames)
    output_img.to_filename(output_img_filename)


def do_percent_signal(WD, WD_template, mni, atlas_filename, stat_filename):
    print(WD)
    os.chdir(WD)
    mni_filename = mni[0]
    warp_filename = os.path.join(WD_template, 'warps.nii.gz') 
    inv_warp_filename = os.path.join(WD_template, 'inv_warps.nii.gz') 
    warped_atlas_filename = 'analysis.feat/warped_atlas.nii.gz'

    # invert warp
    if not os.path.exists(inv_warp_filename):
        os.system('invwarp -w %s -r %s -o %s'%(warp_filename, mni_filename, 
            inv_warp_filename))
            
    # register atlas 
    aw = fsl.ApplyWarp() 
    aw.run(in_file = atlas_filename,
    ref_file = 'analysis.feat/example_func.nii.gz',
        field_file = inv_warp_filename, 
        out_file = warped_atlas_filename, 
        postmat = 'analysis.feat/reg/highres2example_func.mat',
        interp = 'nn')

    # featquery and read output
    atlas = image.load_img(warped_atlas_filename)
    nvols = np.max(atlas.get_data())
    mask_filename = os.path.abspath('querymask.nii.gz') 
    ps = []
    for i in range(nvols):
        mask =  image.new_img_like(atlas, atlas.get_data() == i + 1) 
        mask.to_filename(mask_filename)
        os.system(('featquery 1 analysis.feat 1 stats/%s featquery -p %s')%
            (stat_filename, mask_filename))
        values = pd.read_csv('analysis.feat/featquery/report.txt', 
            sep=' ', header = None)
        ps.append(values[6][0]) 
        shutil.rmtree('analysis.feat/featquery')
        os.system('rm ' + mask_filename)
    return(ps) 
       


def loop_percent_signal(WD, subject, task, suffix, stat_filename, 
    atlas_filename, mni):

    os.chdir(os.path.join(WD, subject))
    sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

    indices = []
    for session, run, path in zip(sessions, runs, paths): 
        proc_dir = os.path.join(WD, subject, session, 'funcproc_' + task + 
            '-' + suffix, run) 
        WD_template = os.path.join(WD, subject, 'average/anat')                
        abs_stat_filename = os.path.join(proc_dir, 'analysis.feat/stats',
            stat_filename + '.nii.gz')

        if os.path.exists(abs_stat_filename):  
            fd = np.mean(pd.read_csv(os.path.join(WD, subject, session, 
                'funcproc_' + task, run, 'fd.txt')).values)
            ps = do_percent_signal(proc_dir, WD_template, mni, atlas_filename, 
                stat_filename)
            indices.append((subject, int(session.split('-')[1]), 
                int(run.split('-')[1]), fd) + tuple(ps))
    return(indices)

def get_percent_signal(WD, subjects, task, suffix, stat_filename, 
    output_filename, atlas_filename, atlas_names, mni, NJOBS):
    """
    Get percent signal with featquery. 
    """

    names = pd.read_csv(atlas_names, header = None).values.tolist()

    indices = Parallel(n_jobs = NJOBS)(delayed(loop_percent_signal)(WD, subject, 
        task, suffix, stat_filename, atlas_filename, mni)
    for subject in subjects) 

    os.chdir(WD)

    indices = [item for sublist in indices for item in sublist]
    indices = pd.DataFrame(indices)
    indices.columns = ['SUBJECT','SESSION','RUN','FD'] + \
        [ x[0] for x in names ]
    #'ROI_%03d'%(i+1)  for i in range(indices.shape[1] - 4)] # #
    indices.to_csv(output_filename, index = False)


def get_ICA_lists(WD, subjects, task, suffix, NVOLS):
    """
    Produce image lists for ICA analysis etc
    """
    indices = []
    func_filenames = []
    struct_filenames = []

    for subject in subjects:
        os.chdir(os.path.join(WD, subject))
        sessions, runs, paths = iterate_sessions_runs(WD, subject, task)

        for session, run, path in zip(sessions, runs, paths): 
            print(path)
            if NVOLS > 0:
                img = image.load_img(path)
                print(img.shape[3])

                if img.shape[3] != NVOLS: 
                    print('No')
                    continue
            proc_dir = os.path.join(WD, subject, session, 'icaproc_' + task, 
                run)
            if not os.path.exists(proc_dir):
                os.makedirs(proc_dir)
            link = os.path.join(proc_dir, 'fMRI.nii.gz')
            if os.path.exists(link):
                os.remove(link)        
            os.symlink(os.path.abspath(path), link) 

            func_filenames.append(link)
            struct_filenames.append(os.path.join(WD, subject, 'average/anat',
                 'structural_brain.nii.gz'))

            fd = np.mean(pd.read_csv(os.path.join(WD, subject, session, 
                'funcproc_' + task, run, 'fd.txt')).values)
            indices.append((subject, int(session.split('-')[1]), 
                int(run.split('-')[1]), fd) )


    os.chdir(WD)

    indices = pd.DataFrame(indices)
    indices.columns = ['SUBJECT','SESSION','RUN','FD']
    indices.to_csv('ICA_table.csv', index = False)

    with open('ICA_func.txt','wb') as f:
        for name in func_filenames:
            f.write(name + '\n')

    with open('ICA_struct.txt','wb') as f:
        for name in struct_filenames:
            f.write(name + '\n')



