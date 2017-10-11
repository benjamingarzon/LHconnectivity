
def register(WD, subject, session, run, structural_filename, TR):

    """
    Compute and apply registrations. 
    """

    WD_subject = WD + subject + '/' 
    os.chdir(WD_subject)

    # brain extraction
    if not os.path.exists('structural_wmseg.nii.gz'):
#        T1images = '/home/share/LeftHand/Structural/' + subject + '.*/unprocessed/3T/T1w_MPR1/LH*.nii.gz' 
#        os.system('mri_robust_template --mov `echo ' + T1images + '` --template structural.nii.gz --satit' )
        os.system('fslroi ' + structural_filename + 
            ' structural.nii.gz 0 -1 0 -1 50 450' )
        os.system('fslmaths structural.nii.gz -thr 0 structural.nii.gz' )
        os.system('flirt -in structural.nii.gz -ref structural.nii.gz ' + 
            '-out structural.nii.gz -applyisoxfm 1')
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

    #fMRI
    # motion correction and parameters
    os.chdir(WD + subject + '.' + session + '/run' + run)

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
                mask = True)

    if not os.path.exists('fd.txt'):
        os.system('fsl_motion_outliers -i fMRI.nii.gz -m brain_mask.nii.gz ' +
            '-o confounds.txt -s fd.txt --fd')

    # register fMRI
    if not os.path.exists('BOLD.mat'):
        epr = fsl.EpiReg()
        epr.run(epi = 'mean.nii.gz', 
                t1_head = WD_subject + 'structural.nii.gz',
                t1_brain = WD_subject + 'structural_brain.nii.gz',
                out_base = 'epireg',
                wmseg = WD_subject + 'structural_wmseg.nii.gz')
        shutil.copyfile('epireg.mat', 'BOLD.mat')
        os.system('rm epireg*')

    
    if not os.path.exists(func_mni_filename):
        aw = fsl.ApplyWarp() 
        aw.run(in_file = 'fMRI_mcf.nii.gz',
               ref_file = mni_filename,
               field_file = WD_subject + 'warps.nii.gz', 
               out_file = func_mni_filename, 
               premat = 'BOLD.mat')

        #get average volume to check that everything is alrighty
        os.system('fslmaths ' + func_mni_filename + 
            ' -Tmean fMRI_mni_mean.nii.gz')

