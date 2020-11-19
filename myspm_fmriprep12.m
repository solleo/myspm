function job = myspm_fmriprep12 (job)
% job = myspm_fmriprep12 (job)
%
% For preprocessing of EPI and T1w images. This script is based on a batch
% script included in SPM12.
% It does slice-timing-correction before realignment/unwarping.
%
% This does:
%  (1) fname_t1w -> unified segmentation+normalization -> brain-masking
%      (skull-stripping) -> normalize into mni
%  (2) fname_pha/fname_mag -> construct VDM
%  (3) fname_epi -> slice timing correction -> unified unwarp+realign
%      -> coregister into t1w (only header) -> denoising (+compcor)
%      -> normalize into mni (resampling) -> smoothing (default: fwhm=2.5 vox)
%
% JOB requires:
% [data]
%  .fname_epi   : a string for a single session
% (.fname_json) : (or with the same filename as the .nii file to be found)
%  .fname_t1w   : a string for anatomical image (full path)
%
% [Cropping]
% (.dummy_sec)  : removing first K sec for unsteady state (adding suffix in TR)
%
% [Slice timing correction]
% * TR_sec, slice_oder_msec, ref_slices_msec can be read from a paired json
% file of fname_epi
% (.slice_order_msec) : can be 1-based indices or actual time from the pulse (msec)
% (.ref_slice_msec)   : reference slice (in the same unit as .slice_order_msec)
% (.fname_dcm)  : DICOM file of a volume to read slice time (only for SIEMENS scanners)
% (.stc_fsl)    : use FSL instead of SPM for a very big file
% (.noSTC)      : for a sparse sampling (e.g. TR>6) STC would make no effect
%                 or could introduce a strange interpolation between distant time points
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
%
% [Smoothing]
% (.fwhm_mm)    : smoothing in MNI152 space
%
% [Misc.]
% (.dir_fig)    : a directory to save head motion plots of multiple subjects together
% (.useants)    : 0 [default] | 1 (takes LOOONG: need MYANTS wrappers)
% (.overwrite)  : 0 [default] or 1
%
% This script uses original and modified MATLAB codes by others:
% - NIFTI toolbox by Jimmy SHEN:
% - DPARSF by Chao-Gan YAN: http://rfmri.org/DPARSF
%
% (cc) 2015-2020. sgKIM. mailto:solleo@gmail.com  https://ggooo.wordpress.com

if nargin==0, help(mfilename); return; end

if isempty(getCurrentWorker)
  a = spm('Ver');
  if ~strcmp(a(4:5),'12')
    error('Run this script on SPM12!');
  end
  % check NIFTI toolbox
  if isempty(which('load_untouch_nii'))
    web('https://mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image')
    error('This script uses NIFTI MATLAB toolbox because it is much faster than spm_vol or MRIread. Please download and path the toolbox.')
  end
end

if ~isfield(job,'overwrite'), job.overwrite =0; end
if isfield(job,'dname_data'), cd(job.dname_data); end
spm('Defaults','fmri');
spm_jobman('initcfg');


%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  m${t1w}.nii     : bias-corrected t1w
% [5]  bm${t1w}.nii    : skull-stripped brain image
% [6]  wbm${t1w}.nii   : skull-stripped brain image registered in mni152
[p2,f2,e2]=myfileparts(job.fname_t1w);
job.fname_t1w=[p2,'/',f2,e2];
fname_in = job.fname_t1w;
ls(fname_in);
fname_out = [p2,'/bm',f2,e2];
fname_def = [p2,'/y_',f2,e2];
if ~exist(fname_out,'file') || job.overwrite
  disp('[1] Unified segmentation of T1w..');
  job1 = job;
  job1.norm = 1;
  if ~isfield(job,'seg')
    job.seg = 'spm12';
  end
  if strcmp(job.seg,'spm12')
    myspm_seg12(job1,'ss')
  elseif strcmp(job.seg,'cat12')
    myspm_seg_cat12(job1,'ss')
  end
  isdone(fname_out,1);
end


%% ANTS registration on skull-stripped T1w to skull-stripped MNI template
if job.useants
  if ~exist('mni_brain.nii.gz','file')
    system(['ln -sf ',fullfile(getenv('FSLDIR'),'data','standard',...
      'MNI152_T1_1mm_brain.nii.gz'),...
      ' mni_brain.nii.gz']);
  end
  cfg = [];
  cfg.fname_fixed = 'mni_brain.nii.gz';
  cfg.fname_moving = [p2,'/bm',f2,e2];
  cfg.reg_stages = 2;
  cfg = myants_antsRegistration(cfg);
  fn_reg_t1w_to_mni = cfg.fname_reg;
end


%% ---EPI processing---
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

%% Cropping dummy scans (for old scanners)
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

%% find slice timing
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
    warning(['TR is long (>= 6 s) and found no dicom file (.fname_dcm) to read actual slice time. Thus setting .noSTC=1']);
    disp(['[!] It may possible to use vectors in .slice_order_msec and .ref_slice_msec']);
    job.noSTC=1;
  end
  % find slice timing in msec and repetition time in sec from a example DICOM
  if isfield(job,'fname_dcm')
    hdr2 = spm_dicom_headers(job.fname_dcm);
    if isfield(hdr2{1},'Private_0019_1029') % recent Siemens scanners
      slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
    else
      if ~isfield(job,'slice_order')
        error('Slice timing information cannot be found in the DICOM header!')
      end
    end
    ref_slice_msec = job.TR_sec*1000/2;
  elseif isfield(job,'slice_order_msec')
    disp(['Found .slice_order_msec; applying those values']);
    slice_order_msec = job.slice_order_msec;
    ref_slice_msec = job.TR_sec*1000/2;
    job.noSTC = 0;
  else
    warning('NO slice timing information is given. Creating a link instead of STC...');
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


%% 3. slice timing correction
% output:
% [1]  a${epi}.nii   : slice-timing corrected EPI
ls(fname_epi);
[p1,f1,e1] = myfileparts(fname_epi);
fname_output = [p1,'/a',f1,e1];
if ~exist(fname_output,'file') || job.overwrite
  disp('[3] slice timing correction..');
  disp(['> Slice order (msec) = ',num2str(reshape(slice_order_msec,1,[]))]);
  disp(['> Reference slice = ',num2str(ref_slice_msec),' msec']);
  if isfield(job,'noSTC') && job.noSTC
    disp(['Create a link for ',fname_epi,' as ',fname_output,'..']);
    unix(['ln -sf ',fname_epi,' ',fname_output])
  else
    if isfield(job,'stc_fsl')&&job.stc_fsl
      if isfield(job,'fname_dcm')
        hdr2 = spm_dicom_headers(job.fname_dcm);
        slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
        T = slice_order_msec'./(job.TR_sec*1000); % in TR
      else
        T = (job.slice_order_msec-1)./max(job.slice_order_msec);
      end
      fname_stc = [p1,'/',f1,'_timing.txt'];
      dlmwrite(fname_stc,T)
      unix(['FSLOUTPUTTYPE=NIFTI; slicetimer.fsl -i ',fname_epi,...
        ' -o ',fname_output,...
        ' --tcustom=',fname_stc,' -r ',num2str(job.TR_sec),' -v -d 3'])
    else
      st1.scans={};
      for t = 1:NumFrames
        st1.scans{1}{t,1} = [fname_epi,',',num2str(t)];
      end
      st1.nslices = numel(slice_order_msec);
      st1.tr = job.TR_sec;
      st1.ta = job.TR_sec-(job.TR_sec/st1.nslices); % NOTE: this value will
      % not be used because .so and .refslice are in msec (see SPM Batch
      % Editor help)
      st1.so = slice_order_msec;
      st1.refslice = ref_slice_msec;
      st1.prefix = 'a';
      matlabbatch = {};
      matlabbatch{1}.spm.temporal.st = st1;
      fname_matlabbatch=[p1,'/myspm_fmriprep3_',f1,'.mat'];
      save(fname_matlabbatch,'matlabbatch');
      spm_jobman('run', matlabbatch);
    end
    isdone(fname_out,3);
  end
end


%% 4a. preparing VDM
if isfield(job,'fname_mag') && isfield(job,'fname_pha')
  if ~isfield(job,'TEs_fmap')
    [p5,f5,~] = myfileparts(job.fname_mag);
    try
      json1 = readjson([p5,'/',f5,'.json']);
      [p5,f5,~] = myfileparts(job.fname_pha);
      json2 = readjson([p5,'/',f5,'.json']);
      job.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      %         JOB.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  job.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  [~,f_ph,e_ph] = myfileparts(job.fname_pha);
  job.fname_vdm = [p1,'/vdm5_sc',f_ph,e_ph];
  if ~exist(job.fname_vdm,'file')
    myspm_prepare_vdm(job.fname_mag, job.fname_pha, job.TEs_fmap, fname_epi, ...
      job.totalreadout_msec); % coreg of t1w here does not help much...
  end
  ls(job.fname_vdm)
end


%% 4. unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_a${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  a${epi}.mat      : [4x4xT] realign transform
% [3]  a${epi}_uw.mat   : unwarping meta data
% [4]  ua${epi}.nii     : unwarped/realigned image
% [5]  meanua${epi}.nii : mean image of [4]
realignunwarp1 = [];
for t = 1:NumFrames
  realignunwarp1.data.scans{t,1} = [p1,'/a',f1,e1,',',num2str(t)];
end
if isfield(job,'fname_vdm')
  ls(job.fname_vdm)
  realignunwarp1.data.pmscan = {[job.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(job,'realign_sep')
  realignunwarp1.eoptions.sep = job.realign_sep;
end
realignunwarp1.eoptions.fwhm = 5;
realignunwarp1.eoptions.rtm = 1; % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
realignunwarp1.eoptions.einterp = 2;
realignunwarp1.eoptions.ewrap = [0 0 0];
realignunwarp1.eoptions.weight = '';
realignunwarp1.uweoptions.basfcn = [12 12];
realignunwarp1.uweoptions.regorder = 1;
realignunwarp1.uweoptions.lambda = 100000;
realignunwarp1.uweoptions.jm = 0;
realignunwarp1.uweoptions.fot = [4 5];
realignunwarp1.uweoptions.sot = [];
realignunwarp1.uweoptions.uwfwhm = 4;
realignunwarp1.uweoptions.rem = 1;
realignunwarp1.uweoptions.noi = 5;
realignunwarp1.uweoptions.expround = 'Average';
realignunwarp1.uwroptions.uwwhich = [2 1];
realignunwarp1.uwroptions.rinterp = 4;
realignunwarp1.uwroptions.wrap = [0 0 0];
realignunwarp1.uwroptions.mask = 1;
realignunwarp1.uwroptions.prefix = 'u';
matlabbatch = {};
matlabbatch{1}.spm.spatial.realignunwarp = realignunwarp1;
fname_in = [p1,'/a',f1,e1];
ls(fname_in);
fname_output = [p1,'/ua',f1,e1];
if ~exist(fname_output,'file') || job.overwrite
  disp('[4] Unwarp & realign..');
  fname_matlabbatch = [p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  
  % visualize unwarping results:
  [p2,f2,e2] = myfileparts(job.fname_t1w);
  if strcmp(f2,'uni')
    fn_t1w = [p2,'/bm',f2,e2];
  else
    fn_t1w = [p2,'/',f2,e2];
  end
  compare_unwarped([p1,'/a',f1,e1], [p1,'/ua',f1,e1], fn_t1w)
  isdone(fname_out,4);
end


%% 5. Intensity-bias correction of EPI (for stable coregistration)
fname_epi_unbiased = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
if ~exist(fname_epi_unbiased,'file')
  unix(['mri_nu_correct.mni --i ',p1,'/meanua',f1,e1,' --o ',fname_epi_unbiased]);
end


%% 6. Coregistration of EPI to native T1w
if ~job.useants
  
  % outputs: <modifying transform matrices in headers>
  % [1]  mmeanua${epi}.nii
  % [2]  ua${epi}.nii        : "other" images, header modified
  % [3]  ua${epi}.mat        : "other" images, header modified (like
  %                            realignment, affine transform matrices for all
  %                            volumes (4x4x#vols)
  estimate1=[];
  fname_t1w_brain     = [p2,'/bm',f2,e2]; % bias-corrected T1w
  estimate1.ref{1}    = fname_t1w_brain;
  ls(fname_t1w_brain)
  fname_epi_unbiased  = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
  estimate1.source{1} = fname_epi_unbiased ;
  ls(fname_epi_unbiased )
  % ------------------here you specified OTHER imgaes-----------------
  fname_epi_ua = [p1,'/ua',f1,e1];
  ls(fname_epi_ua);
  estimate1.other{1} = fname_epi_ua; % NIFTI-format have only ONE xfm for all volumes (no need to put all frames here)
  % ------------------------------------------------------------------
  estimate1.eoptions.cost_fun = 'nmi';
  estimate1.eoptions.sep = [4 2];
  estimate1.eoptions.tol = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
  estimate1.eoptions.fwhm = [7 7];
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.estimate = estimate1;
  fname_out=[p1,'/ua',f1,'.mat'];
  if ~exist(fname_out,'file') || job.overwrite
    disp('[6] Coregistration of EPI to native T1w..');
    fname_matlabbatch=[p1,'/myspm_fmriprep6_',f1,'.mat'];
    save(fname_matlabbatch,'matlabbatch');
    spm_jobman('run', matlabbatch);
    isdone(fname_out,6);
  end
  
  % resample the first volume (to check coregistration quality)
  [p4,f4,e4]=myfileparts(fname_epi_unbiased);
  if ~exist([p4,'/r',f4,e4],'file')  || job.overwrite
    matlabbatch={};
    matlabbatch{1}.spm.spatial.coreg.write.ref{1}    = fname_t1w_brain;
    matlabbatch{1}.spm.spatial.coreg.write.source{1} = [fname_epi_unbiased,',1'];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
    fname_epi_in_t1w = [p1,'/rmmeanua',f1,e1];
    slices(fname_epi_in_t1w, [], ...
      struct('fname_png',[p1,'/rmmeanua',f1,'_in_',f2,'.png'],...
      'contour',fname_t1w_brain))
  end
else
  % RIGID from EPI to T1w
  cfg = [];
  cfg.fname_moving = fname_epi_unbiased;
  cfg.fname_fixed  = fname_t1w_brain;
  cfg.reg_stages = 0;
  cfg = myants_antsRegistration(cfg);
  fn_reg_epi_to_t1w = cfg.fname_reg;
end


%% 8. Apply forward deformation on EPI for registration into MNI152
if ~isfield(job,'vox_mm')
  job.vox_mm = hdr.dime.pixdim(2:4);
  disp(['Resampling at [',num2str(job.vox_mm),'] mm in MNI152 space..'])
end
if numel(job.vox_mm) == 1
  job.vox_mm=[1 1 1]*job.vox_mm;
end

if ~job.useants
  matlabbatch={};
  write1=[];
  write1.subj.def{1} = fname_def;
  write1.subj.resample{1} = fname_epi_ua;
  write1.woptions.bb = [-78 -112 -70; 78 76 85];
  
  write1.woptions.vox = job.vox_mm;
  write1.woptions.interp = 4;
  matlabbatch{1}.spm.spatial.normalise.write = write1;
  fname_out=[p1,'/wua',f1,e1];
  if ~exist(fname_out,'file') || job.overwrite
    disp('[8] Registration of EPI to MNI152..');
    fname_matlabbatch=[p1,'/myspm_fmriprep8_',f1,'.mat'];
    save(fname_matlabbatch,'matlabbatch');
    spm_jobman('run', matlabbatch);
    isdone(fname_out,8);
  end
  % quality check:
  fname_png=[p1,'/wua',f1,'_in_mni152.png'];
  fname_epi_in_mni1=[p1,'/wua',f1,'1.nii.gz'];
  setenv('FSLOUTPUTTYPE','NIFTI_GZ');
  unix(['fslroi ',p1,'/wua',f1,e1,' ',fname_epi_in_mni1,' 0 1']);
  ls(fname_epi_in_mni1);
  fname_mni = [getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm.nii.gz'];
  slices(fname_epi_in_mni1, [], ...
    struct('fname_png',fname_png,'contour',fname_mni));
else
  % CREATE FUNCTIONAL REFERENCE AT GIVEN VOXE RESOLUTION
  system(['mri_convert -vs ',num2str(job.vox_mm),...
    ' mni_brain.nii.gz mni_funcref.nii.gz'])
  
  % COMBINE TRANSFORMS
  fn_warp = ['mmeanua',f1,'_to_mni_Warping.nii.gz'];
  cfg = struct('fname_out',fn_reg_epi_to_mni, ...
    'fname_fixed','mni_funcref.nii.gz',...
    'transforms',{{fn_reg_epi_to_t1w,0;fn_reg_t1w_to_mni,0}});
  myants_combinetransforms(cfg)
  
  % APPLY ON TIMESERIES
  cfg = struct('fname_moving',[p1,'/ua',f1,e1], ...
  'fname_fixed','mni_funcref.nii.gz', 'transforms',{{fn_warp,0}}, ...
  'fname_out',[p1,'/ua',f1,'_mni',e1]);
  myants_antsApplyTransformsTimeseries(cfg)
end

%% (10). Smoothing
if isfield(job,'fwhm_mm')
  if ~job.useants
    myspm_smooth(struct('fname',[p1,'/wua',f1,e1],'fwhm_mm',job.fwhm_mm));
  else
    myspm_smooth(struct('fname',[p1,'/wua',f1,'_mni,',e1],...
      'fwhm_mm',job.fwhm_mm));
  end
end


end

%% Subroutines:
function isdone(fname_out,proci)
if exist(fname_out,'file')
  disp(['[',num2str(proci),'] done: ',fname_out])
else
  error(['[',num2str(proci),'] failed: ',fname_out,' was not created'])
end
end
