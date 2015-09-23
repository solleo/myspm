function EXP = myspm_rsfc(EXP)

EXP.subjID = fsss_subjID(EXP.subjID);

%% 1. fmri preprocessing: STC + unwarp/realignment
EXP = myspm_fmriprep12(EXP);

EXP1=[];
EXP1.subjID = subjID;
EXP1.fname_proc={'rest410.nii','arest410.nii','uarest410.nii'};
EXP1.name_proc={'orig','stc','u+a'};
EXP1.dir_base='/scr/vatikan3/APConn/rest12';
EXP1.dir_png='/scr/vatikan3/APConn/rest12/fig_ua';
myspm_check_timeseries(EXP1);

%% 2. motion assessment
EXP = myspm_art(EXP); % to make sure that you don't need scrubbing

%% 3. compcor & filtering using Yan's implementation (too long?)
% find slice timing in msec and repetition time in sec from a example DICOM

EXP = myy_compcor_pbf(EXP);

end