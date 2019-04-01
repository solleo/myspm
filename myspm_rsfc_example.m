function EXP = myspm_rsfc_example(EXP)


EXP.subjID = fsss_subjID(EXP.subjID);

%% 1. fmri/t1w preprocessing: STC + unwarp/realignment (rigid-motion param)
%  .fwhm_mm
%  .fname_epi
%  .fname_t1w
%  .fname_vdm
EXP = myspm_fmriprep12_func(EXP); % 'a' for STC; 'u' for unwarping/realign

%% 2. Scrubbing regressor by ART

EXP=[];
EXP.dir_base='/scr/vatikan3/APConn/rest12.410/';
EXP.subjID = subjID17;
EXP.name_epi='uarest410.nii';
EXP.name_rp ='rp_arest410.txt';
EXP.dir_figure='/scr/vatikan3/APConn/rest12.410/fig_denoising';
EXP.global_threshold=3;
EXP.motion_threshold=0.5;
EXP = myspm_art(EXP); % to make sure that you don't need scrubbing

%% 3. Compcor regressor
EXP.bpf1=[0 inf];
EXP=myy_compcor(EXP);

%% 4. compare results from various regressions

EXP.param_cc='wmcsf99_n16d1v1b0.00-Inf';
EXP.name_cc=['cc_',EXP.param_cc,'_eigenvec.txt'];
EXP.param_art='3.0std_0.5mm';
EXP.name_art=['art_regression_outliers_and_movement_uarest410_',EXP.param_art,'.mat'];
EXP.name_rp='rp_arest410.txt';
EXP.fname_gmmask='oc1t1w_99.nii';
EXP.bpf2=[0.01 0.10];
EXP.covset=[1 2 3 4];
EXP.cov_idx=4; % i'll go with global-signal!
EXP = myspm_denoise(EXP);  % 'r' for residual; 'f' for filtering

%% 5. bring parcellation into function

EXP=[];
EXP.dir_base = '/scr/vatikan3/APConn/rest12.410/';
EXP.dir_fs   = '/scr/vatikan3/APConn/FSspm12';
EXP.name_t1w = 'Brain.nii';
EXP.name_epi = 'meanuarest410.nii';
EXP.subjID = subjID;
myspm_aparc(EXP)

%% 6. compute correlation without smoothing
% smoothing (doesn't seem to be a good idea for ROIs. Averaging already
% cancels out tiny noise. Even 2-voxel-fhwm increases skewness and drops
% gof from ks-test.... (but i'll use it for seed-based correlation!@
EXP=[];
EXP.dir_base = '/scr/vatikan3/APConn/rest12.410/';
EXP.subjID = subjID;
EXP.dir_figure = '/scr/vatikan3/APConn/figures.local/crrmtx_sg';
EXP.name_res='fruarest410.nii';
myspm_crrmtx(EXP)




%% also FS-projection...


%% 7. register only 3D image! 
% (registering 4-D double images TAKES TOO TOO LONG! 10 mins per subject)
% coreg 2.3mm,410scans,single func => strc (2mm): takes 3 min
% deform




end
