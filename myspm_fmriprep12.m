function job = myspm_fmriprep12 (job)
% job = myspm_fmriprep12 (job)
%
% For preprocessing of EPI and T1w images. This script is based on a batch
% script included in SPM12.
% It does slice-timing-correction before realignment/unwarping.
%
% CURRENTLY, THIS SCRIPT REQUIRES FSL AND FREESURFER TOO FOR TINY BUT
% VERY CONVENIENT UTILITIES.
%
% This does:
%  (1) STRUCTURAL DATA (fname_t1w)
%      -> unified segmentation+normalization -> brain-masking
%      -> normalize into mni
%  (2) FIELDMAP DATA (fname_pha/fname_mag)
%      -> construct VDM
%  (3) FUNCTIONAL DATA (fname_epi)
%      -> slice timing correction -> unified unwarp+realign
%      -> coregister into t1w (only header)
%      -> normalize into mni (resampling) -> smoothing (if defined)
%
% JOB requires:
% [data]
%  .fname_epi   : a string for a single session
% (.fname_json) : (or with the same filename as the .nii file to be found)
%  .fname_t1w   : a string for anatomical image (full path)
%
% [Cropping]
% (.dummy_sec)  : removing first K sec for unsteady state (adding suffix in
%                 TR)
%
% [Slice timing correction]
% * TR_sec, slice_oder_msec, ref_slices_msec can be read from a paired json
% file of fname_epi
%
% (.slice_order_msec) : can be 1-based indices or actual time from the
%                       pulse (msec)
% (.ref_slice_msec)   : reference slice (in the same unit as
%                       .slice_order_msec)
% (.fname_dcm)  : DICOM file of a volume to read slice time (only for
%                 SIEMENS scanners)
% (.stc_fsl)    : use FSL instead of SPM for a very big file that SPM can't
%                 handle! (But are you sure you're doing right with the
%                 nortorious FSL slice-timing syntax?)
% (.noSTC)      : for a sparse sampling (e.g. TR>TA) STC would make no
%                 effect or could introduce a strange interpolation between
%                 distant time points?
% (.TR_sec)     : can be read from a paired .json file.
%
% [Fieldmap: geometric distortion correction]
% (.fname_mag)  : magnitude a short TE
% (.fname_pha)  : phase difference
% * TEs_fmap and totalreadout_msec can be read from paired json files of
% fname_mag and fname_pha
% (.TE_fmap)
% (.totalreadout_msec)
% (.fname_vdm)  : if it's already computed
%
% [Unwarping+realignment]
%
%
%
% [Spatial normalization]
% (.vox_mm)     : resampling resolution for EPI in MNI152 space
% (.bbox_mm)
%
% [Smoothing]
% (.fwhm_mm)    : smoothing in MNI152 space (default = none)
%
% [Misc.]
% (.useants)    : 0 [default] | 1 (takes LOOONG: need MYANTS wrappers)
% (.overwrite)  : 0 [default] or 1
%
% This script uses original and modified MATLAB codes by others:
% - NIFTI toolbox by Jimmy SHEN:
% - DPARSF by Chao-Gan YAN: http://rfmri.org/DPARSF
%
% (cc) 2015-2020. sgKIM. mailto:solleo@gmail.com  https://ggooo.wordpress.com

if nargin==0, help(mfilename); return; end
if isempty(getCurrentWorker) % CHECK ONLY WHEN NOT RUNNING IN PAR
  a = spm('Ver');
  if ~strcmp(a(4:5),'12')
    error('Run this script on SPM12!');
  end
  
  % check NIFTI toolbox
  if isempty(which('load_untouch_nii'))
    web(['https://mathworks.com/matlabcentral/fileexchange/',...
      '8797-tools-for-nifti-and-analyze-image']);
    error(['This script uses NIFTI MATLAB toolbox because it is much',...
      'faster than spm_vol or MRIread. Please download and path ',...
      'the toolbox.']);
  end
end
spm('Defaults','fmri');
spm_jobman('initcfg');
if ~isfield(job,'useants'), job.useants = 0; end
if ~isfield(job,'fwhm_mm'), job.fwhm_mm = 0; end
if job.useants
  nprocs = 6 + job.fwhm_mm;
else
  nprocs = 5 + job.fwhm_mm;
end

if ~isfield(job,'overwrite'), job.overwrite =0; end
if ~isfield(job,'dname_data')
  job.dname_data = myfileparts(job.fname_epi);
end
cd(job.dname_data)

%% KEEP DIARY
diary('myspm_fmriprep12.log') % OPEN
disp(repmat('=',[1 72]))
fprintf('[%s:%s] START\n', mfilename, datestr(now,31))
disp(job)
disp(repmat('=',[1 72]))  

%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  m${t1w}.nii     : bias-corrected t1w
% [5]  bm${t1w}.nii    : skull-stripped brain image
% [6]  wbm${t1w}.nii   : skull-stripped brain image registered in mni152
pid = 1;
[p2, f2, e2] = myfileparts(job.fname_t1w);
job.fname_t1w = [p2,'/',f2,e2];
fname_in = job.fname_t1w;
ls(fname_in);
fname_out = [p2,'/bm',f2,e2];
fname_def = [p2,'/y_',f2,e2];
if ~exist(fname_out,'file') || job.overwrite
  fprintf('[%s:%i/%i] Unified segmentation of T1w..\n',...
    mfilename, pid, nprocs)
  job1 = job;
  job1.norm = 1;
  if ~isfield(job,'seg')
    job.seg = 'spm12';
  end
  try
    if strcmp(job.seg,'spm12')
      myspm_seg12(job1,'ss');
    elseif strcmp(job.seg,'cat12')
      myspm_seg_cat12(job1,'ss');
    end
  end
end
isdone(fname_out, pid, nprocs);
pid = pid+1;


%% ANTS registration on skull-stripped T1w to skull-stripped MNI template
if job.useants
  % outputs:
  % [1]  ${t1w}_to_mni_brain_stage2_Composite.h5: deformation field
  % [2]  ${t1w}_to_mni_brain_stage2_fwdwarped.nii.gz: %{t1w} in MNI
  % [3]  ${t1w}_to_mni_brain_stage2_invwarped.nii.gz: MNI in ${t1w}
  
  if ~exist('mni_brain.nii.gz','file')
    system(['ln -sf ',fullfile(getenv('FSLDIR'),'data','standard',...
      'MNI152_T1_1mm_brain.nii.gz'),...
      ' mni_brain.nii.gz']);
  end
  job1 = [];
  job1.fname_fixed = 'mni_brain.nii.gz';
  job1.fname_moving = [p2,'/bm',f2,e2];
  job1.reg_stages = 2;
  fname_out = [p2,'/bm',f2,'_to_mni_brain_stage2_Composite.h5'];
  if ~isfile(fname_out)
    fprintf('[%s:%i/%i] skull-stripped T1w to MNI152-brain using ANTs..\n',...
      mfilename, pid, nprocs)
    try
      job1 = myants_antsRegistration(job1);
    end
  end
  isdone(fname_out, pid, nprocs);
  pid = pid+1;
end


%% ---EPI processing---

%% Read parameters from JSON, DICOM, NIFTI headers or user inputs:
if ~isfield(job,'fname_json')
  [path1,name1,~] = myfileparts(job.fname_epi);
  job.fname_json = [path1,'/',name1,'.json'];
end
fname_epi = job.fname_epi;
[path1,name1,ext1] = myfileparts(fname_epi);
fname_epi = [path1,'/',name1,ext1];
ls(fname_epi)

%% find TR
fn_json = job.fname_json;
if ~isfield(job,'TR_sec') % manually given?
  if exist(fn_json,'file') % using json files
    disp(['Found json file: ',fn_json]);
    disp(['Reading TR & SliceTiming from json']);
    json = jsondecode(fileread(fn_json));
    job.TR_sec = json.RepetitionTime; % *** it's SEC!! ***
  elseif isfield(job,'fname_dcm')
    if exist(fn_dcm,'file')
      disp(['Found dicom file: ',fname_dcm]);
      hdr2 = spm_dicom_headers(job.fname_dcm);
      job.TR_sec  = hdr2{1}.RepetitionTime/1000; % it's msec!!!
    end
  else
    error('TR cannot be determined!')
  end
end
fprintf('> TR = %i msec\n',job.TR_sec*1000);

%% CROP dummy scans (for old scanners or new protocols)
if isfield(job,'dummy_sec')
  fname_epi_old = fname_epi;
  n = ceil(job.dummy_sec / job.TR_sec);
  fprintf('Removing first %i sec (%i volumes)\n',job.dummy_sec, n);
  fname_epi = [path1,'/',name1,'_skip',num2str(n),ext1];
  if ~exist(fname_epi,'file')
    unix(['mri_convert --nskip ',num2str(n),' ',fname_epi_old,' ',fname_epi]);
  end
  ls(fname_epi)
end

%% FIND slice timing
hdr = load_nii_hdr(fname_epi);
NumFrames = hdr.dime.dim(5);
fprintf('> Number of frames = %i\n',NumFrames);

if exist(fn_json,'file')
  if ~contains(json.Manufacturer,'GE')
    slice_order_msec = json.SliceTiming*1000;
    ref_slice_msec = job.TR_sec*1000/2;
  else
    disp('Is this GE data?')
    % NOTE: dcm2niix has no access to GE machine, thus still reading it
    % strangely. Use BIAC header instead:
    if exist(job.fname_bxh,'file')
      disp(['Found BIAC header file: ',job.fname_bxh])
    end
    bxh = parsebxh(job.fname_bxh);
    slice_order_msec = bxh.acquisitiontime_msec;
    ref_slice_msec = job.TR_sec*1000/2;
  end
else
  % if no json
  if job.TR_sec >= 6 && ~isfield(job,'noSTC') && ~isfield(job,'fname_dcm')
    warning(...
      ['TR is long (>= 6 s) and found no dicom file (.fname_dcm)',...
      'to read actual slice time. Thus setting .noSTC=1']);
    disp(['[!] It may possible to use vectors in .slice_order_msec ',...
      'and .ref_slice_msec']);
    job.noSTC=1;
  end
  % find slice timing in msec and repetition time in sec from a example
  % DICOM
  if isfield(job,'fname_dcm')
    hdr2 = spm_dicom_headers(job.fname_dcm);
    if isfield(hdr2{1},'Private_0019_1029') % recent Siemens scanners
      slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
    else
      if ~isfield(job,'slice_order')
        error('No slice timing information in the DICOM header!')
      end
    end
    ref_slice_msec = job.TR_sec*1000/2;
  elseif isfield(job,'slice_order_msec')
    disp('Applying manually entered .slice_order_msec');
    slice_order_msec = job.slice_order_msec;
    ref_slice_msec = job.TR_sec*1000/2;
    job.noSTC = 0;
  else
    warning(['NO slice timing information is given. Skipping STC ', ...
      'CREATING a link...']);
    job.noSTC = 1;
  end
end
% override of manually given:
if isfield(job,'slice_order_msec')
  slice_order_msec = job.slice_order_msec;
end
if isfield(job,'ref_slice_msec')
  ref_slice_msec = job.ref_slice_msec;
end


%% 2. slice timing correction
% output:
% [1]  a${epi}.nii   : slice-timing corrected EPI
ls(fname_epi);
[p1,f1,e1] = myfileparts(fname_epi);
fname_out = [p1,'/a',f1,e1];
if ~exist(fname_out,'file') || job.overwrite
  fprintf('[%s:%i/%i] slice timing correction..\n',...
    mfilename, pid, nprocs)
  disp(['> Slice order (msec) = ',num2str(reshape(slice_order_msec,1,[]))]);
  disp(['> Reference slice = ',num2str(ref_slice_msec),' msec']);
  if isfield(job,'noSTC') && job.noSTC
    disp(['Create a link for ',fname_epi,' as ',fname_out,'..']);
    unix(['ln -sf ',fname_epi,' ',fname_out])
  else
    try
      if isfield(job,'stc_fsl')&&job.stc_fsl
        if isfield(job,'fname_dcm')
          hdr2 = spm_dicom_headers(job.fname_dcm);
          slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
          T = slice_order_msec'./(job.TR_sec*1000); % in TR
        else
          T = (job.slice_order_msec-1)./max(job.slice_order_msec);
        end
        fname_stc = [p1,'/',f1,'_timing.txt'];
        dlmwrite(fname_stc,T);
        unix(['FSLOUTPUTTYPE=NIFTI; slicetimer.fsl -i ',fname_epi,...
          ' -o ',fname_out,...
          ' --tcustom=',fname_stc,' -r ',num2str(job.TR_sec),' -v -d 3'])
      else
        job1 = struct('fname_epi',fname_epi, 'NumFrames',NumFrames, ...
          'TR_sec',job.TR_sec, 'slice_order_msec',slice_order_msec, ...
          'ref_slice_msec',ref_slice_msec);
        myspm_stc(job1);
      end
    end
  end
end
isdone(fname_out, pid, nprocs);
pid = pid+1;

%% 3. unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_a${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  a${epi}.mat      : [4x4xT] realign transform
% [3]  a${epi}_uw.mat   : unwarping meta data
% [4]  ua${epi}.nii     : unwarped/realigned image
% [5]  meanua${epi}.nii : mean image of [4]
job1 = job;
job1.fname_epi = [p1,'/a',f1,e1];
fname_out = [p1,'/meanua',f1,e1];
if ~isfile(fname_out)
  fprintf('[%s:%i/%i] unwarping and realignment..\n',...
    mfilename, pid, nprocs)
  try
    myspm_unwarp(job1);
  end
end
isdone(fname_out, pid, nprocs);
pid = pid+1;

%% 4. Coregistration of EPI to native T1w
% outputs:
% - if SPM12 used:
% [1]  ua${epi}.nii       : header-modified EPI images
% [2]  rmmeanua${epi}.nii : resampled mean EPI image (for quality check)
%
% - if ANTs used:
% [1]  mmeanuaa${epi}_to_bm${t1w}_stage0_Composite.h5
%      : rigid-body transform (by ANTs)

if ~job.useants
  % OUTPUT? resampled volume#1 for sanity check
  fname_out = [p1,'/rmmeanua',f1,'_in_',f2,e1];
  job1 = job;
  job1.fname_epi = [p1,'/ua',f1,e1];
  if ~isfile(fname_out)
    fprintf('[%s:%i/%i] coreg EPI to T1w using SPM12..\n',...
      mfilename, pid, nprocs)
    myspm_coreg_hdr(job1);
  end
    
else
  % Intensity-bias correction of EPI (for stable coregistration)
  fname_epi_unbiased = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
  if ~exist(fname_epi_unbiased,'file')
    unix(['mri_nu_correct.mni --i ',p1,'/meanua',f1,e1,...
      ' --o ',fname_epi_unbiased]);
  end
  
  % RIGID from EPI to T1w
  job1 = [];
  job1.fname_moving = [p1,'/mmeanua',f1,e1];  % bias-corrected EPI
  job1.fname_fixed  = [p2,'/bm',f2,e2];  % bias-corrected T1w
  job1.reg_stages = 0;
  fname_out = ['mmeanua',f1,'_to_bm',f2,'_stage0_Composite.h5'];
  if ~isfile(fname_out)
    fprintf('[%s:%i/%i] coreg EPI to T1w using ANTs..\n',...
      mfilename, pid, nprocs)
    myants_antsRegistration(job1);
  end

end
isdone(fname_out, pid, nprocs);
pid = pid+1;

%% 5. Apply forward deformation on EPI for registration into MNI152
% outputs:
% - if SPM12 used:
% [1]  wua${epi}.nii    : resampled EPI images in MNI152 (by SPM12)
%
% - if ANTs used:
% [1]  ua${epi}_mni.nii : resampled EPI images in MNI152 (by ANTs)

if ~isfield(job,'vox_mm')
  job.vox_mm = hdr.dime.pixdim(2:4);
  disp(['Resampling at [',num2str(job.vox_mm),'] mm in MNI152 space..'])
end
if numel(job.vox_mm) == 1
  job.vox_mm = [1 1 1]*job.vox_mm;
end
if ~isfield(job,'bbox_mm')
  job.bbox_mm = [];
end

if ~job.useants
  fname_out=[p1,'/wua',f1,e1];
  if ~exist(fname_out,'file') || job.overwrite
    fprintf('[%s:%i/%i] resampling EPI in MNI152 using SPM12..\n',...
      mfilename, pid, nprocs)
    myspm_norm(struct('fname_moving',[p1,'/ua',f1,e1], ...
      'fname_deform',fname_def, 'vox_mm',job.vox_mm, 'interp',4, ...
      'bbox_mm',job.bbox_mm));
  end
  
  % quality check:
  fname_png = [p1,'/wua',f1,'_in_mni152.png'];
  if ~isfile(fname_png)
    fname_epi_in_mni1 = [p1,'/wua',f1,'1.nii.gz'];
    setenv('FSLOUTPUTTYPE','NIFTI_GZ');
    unix(['fslroi ',p1,'/wua',f1,e1,' ',fname_epi_in_mni1,' 0 1']);
    ls(fname_epi_in_mni1);
    fname_mni = [getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm.nii.gz'];
    slices(fname_epi_in_mni1, [], ...
      struct('fname_png',fname_png,'contour',fname_mni));
  end
else
  % CREATE FUNCTIONAL REFERENCE AT GIVEN VOXE RESOLUTION
  if ~isfile('mni_funcref.nii')
    % reslice
    system(['mri_convert -vs ',num2str(job.vox_mm),...
      ' mni_brain.nii.gz mni_funcref.nii'])
    % bounding-box
    if isfield(job,'bbox_mm') && ~isempty(job.bbox_mm)
      myspm_boundingbox('mni_funcref.nii', job.bbox_mm,...
        'mni_funcref.nii')
    else
      myspm_boundingbox('mni_funcref.nii', 'canon',...
        'mni_funcref.nii')
    end
  end
  
  % COMBINE TRANSFORMS
  fn_warp = ['mmeanua',f1,'_to_mni_Warping.nii.gz'];
  fn_reg_epi_to_t1w = ['mmeanua',f1,'_to_bm',f2,'_stage0_Composite.h5'];
  fn_reg_t1w_to_mni = ['bm',f2,'_to_mni_brain_stage2_Composite.h5'];
  job1 = struct('fname_out',fn_warp, ...
    'fname_fixed','mni_funcref.nii',...
    'transforms',{{fn_reg_epi_to_t1w,0;fn_reg_t1w_to_mni,0}});
  if ~isfile(fn_warp)
    myants_combinetransforms(job1)
  end
  
  % APPLY ON TIMESERIES
  fname_out = [p1,'/xua',f1,e1];
  job1 = struct('fname_moving',[p1,'/ua',f1,e1], ...
    'fname_fixed','mni_funcref.nii', 'transforms',{{fn_warp,0}}, ...
    'fname_out',fname_out);
  if ~isfile(fname_out)
    fprintf('[%s:%i/%i] resampling EPI in MNI152 using ANTs..\n',...
      mfilename, pid, nprocs)
    myants_antsApplyTransformsTimeseries(job1)
  end
    
end
isdone(fname_out, pid, nprocs);
pid = pid+1;


%% (6). Smoothing
if job.fwhm_mm
  if ~job.useants
    fname_in = [p1,'/wua',f1,e1];
    fname_out = [p1,'/s',num2str(job.fwhm_mm(1)),'wua',f1,e1];
  else
    fname_in = [p1,'/xua',f1,e1];
    fname_out = [p1,'/s',num2str(job.fwhm_mm(1)),'xua',f1,e1];
  end
  if ~isfile(fname_out)
    fprintf('[%s:%i/%i] smoothing..\n',mfilename, pid, nprocs);
    myspm_smooth(struct('fname',fname_in,'fwhm_mm',job.fwhm_mm));
  end
  isdone(fname_out, pid, nprocs);
  pid = pid+1;
end

diary off % CLOSE (UPDATE)
end

%% Subroutines:
function isdone(fname_out, pid, nprocs)
if exist(fname_out,'file')
  fprintf('[myspm_fmriprep12:%i/%i] DONE:', pid, nprocs);
  ls(fname_out)
  diary('myspm_fmriprep12.log') % CLOSE (UPDATE)
  diary('myspm_fmriprep12.log') % open (wait)
else
  fprintf('[myspm_fmriprep12:%i/%i] FAILED to create: %s', ...
    pid, nprocs, fname_out);
  warning('PROCESS ABORTED!');
  diary off % CLOSE (UPDATE)
  return
end
end
