function JOB = myspm_fmriprep12withants (JOB)
% JOB = myspm_fmriprep12withants (JOB)
%
% FOR DIFFICULT REGISTRATIONS, this pipeline uses ANTs instead of SPM12:
% (1) coregistration between mean realigned/unwarped EPI to T1w: ants-rigid
% (2) skull-stripped T1w to MNI-brain: ants-SyN
% (3) apply same registration for all volumes in a EPI timeseries
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
% ["Denoising" or motion artifacts suppression]
% (.bpf)        : (default=[1/128 inf] Hz)
% (.restbpf)    : 1=[0.009 0.08] Hz from Satterthwaite.2013
% (.num_pcs)    : 3 [default] # of CompCor regressors to extract
% (.subjid)     : for motion plots and compcor if FS directory is available
% (.cov_idx)    : for motion-artifact removal: rp+cc+sc+gs
% (.out_prefix) : to take a note on threshold selection etc.
%
% [Spatial normalization]
% (.vox_mm)     : resampling resolution for EPI in MNI152 space
%
% [Smoothing]
% (.fwhm_mm)    : smoothing in MNI152 space
%
% [Misc.]
% (.dir_fig)    : a directory to save head motion plots of multiple subjects together
%
%
% (.overwrite)  : 0 [default] or 1
%
% This script uses original and modified MATLAB codes by others:
% - NIFTI toolbox by Jimmy SHEN:
% - DPARSF by Chao-Gan YAN: http://rfmri.org/DPARSF
% - SCREEN2JPEG by Sean P. McCARTHY (MALTAB builtin)
%
% (cc) 2015,2017,2018. sgKIM. mailto:solleo@gmail.com  https://ggooo.wordpress.com

a=spm('version');
if ~strcmp(a(4:5),'12')
  error('Run this script on SPM12!');
end
if nargin==0, help(mfilename); return; end
% check NIFTI toolbox
if isempty(which('load_untouch_nii'))
  web('https://mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image')
  error('This script uses NIFTI MATLAB toolbox because it is much faster than spm_vol or MRIread. Please download and path the toolbox.')
end
if ~isfield(JOB,'overwrite'), JOB.overwrite =0; end
if isfield(JOB,'dname_data'), cd(JOB.dname_data); end
spm('Defaults','fmri');
spm_jobman('initcfg');

% try
%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  m${t1w}.nii     : bias-corrected t1w
% [5]  bm${t1w}.nii    : skull-stripped brain image
% [6]  wbm${t1w}.nii   : skull-stripped brain image registered in mni152
[p2,f2,e2]=myfileparts(JOB.fname_t1w);
JOB.fname_t1w=[p2,'/',f2,e2];
fname_in = JOB.fname_t1w;
ls(fname_in);
fname_out = [p2,'/bm',f2,e2];
fname_def = [p2,'/y_',f2,e2];
if ~exist(fname_out,'file') || JOB.overwrite
  disp('[1] Unified segmentation of T1w..');
  job1 = JOB;
  job1.norm = 1;
  myspm_seg12(job1,'ss')
  isdone(fname_out,1);
end
%% ANTS registration on skull-stripped T1w to skull-stripped MNI template
if ~exist('mni_brain.nii','file')
  copyfile(...
    fullfile(getenv('FSLDIR'),'data','standard','MNI152_T1_1mm_brain.nii.gz'),...
    'mni_brain.nii.gz');
  unix('gunzip mni_brain.nii.gz')
end
job = [];
job.fname_fixed = 'mni_brain.nii';
job.fname_moving = [p2,'/bm',f2,e2];
job.reg_stages = 2;
job = myants_antsRegistration(job);
fn_reg_t1w_to_mni = job.fname_reg;



%% ---EPI processing---
if ~isfield(JOB,'fname_json')
  [path1,name1,~] = myfileparts(JOB.fname_epi);
  JOB.fname_json = [path1,'/',name1,'.json'];
end
fname_epi = JOB.fname_epi;
[path1,name1,ext1] = myfileparts(fname_epi);
fname_epi = [path1,'/',name1,ext1];
ls(fname_epi)
%% find TR
fn_json = JOB.fname_json;
if ~isfield(JOB,'TR_sec') % manually given?
  if exist(fn_json,'file') % using json files
    disp(['Found json file: ',fn_json]);
    disp(['Reading TR & SliceTiming from json']);
    json = jsondecode(fileread(fn_json));
    JOB.TR_sec = json.RepetitionTime; % *** it's SEC!! ***
  end
elseif isfield(JOB,'fname_dcm')
  if exist(fn_dcm,'file')
    disp(['Found dicom file: ',fname_dcm]);
    hdr2 = spm_dicom_headers(JOB.fname_dcm);
    JOB.TR_sec  = hdr2{1}.RepetitionTime/1000; % it's msec!!!
  end
else
  error('TR cannot be determined!')
end
fprintf('> TR = %i msec\n',JOB.TR_sec*1000);

%% Cropping dummy scans (for old scanners)
if isfield(JOB,'dummy_sec')
  fname_epi_old = fname_epi;
  n = ceil(JOB.dummy_sec / JOB.TR_sec);
  fprintf('Removing first %i sec (%i volumes)\n',JOB.dummy_sec, n);
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
  if ~strcmp(json.Manufacturer,'GE')
    slice_order_msec = json.SliceTiming*1000;
    ref_slice_msec = JOB.TR_sec*1000/2;
  else
    bxh = parsebxh(JOB.fname_bxh);
    slice_order_msec = bxh.acquisitiontime_msec;
    ref_slice_msec = JOB.TR_sec*1000/2;
  end
else
  % if no json
  if JOB.TR_sec >= 6 && ~isfield(JOB,'noSTC') && ~isfield(JOB,'fname_dcm')
    warning(['TR is long (>= 6 s) and found no dicom file (.fname_dcm) to read actual slice time. Thus setting .noSTC=1']);
    disp(['[!] It may possible to use vectors in .slice_order_msec and .ref_slice_msec']);
    JOB.noSTC=1;
  end
  % find slice timing in msec and repetition time in sec from a example DICOM
  if isfield(JOB,'fname_dcm')
    hdr2 = spm_dicom_headers(JOB.fname_dcm);
    if isfield(hdr2{1},'Private_0019_1029') % recent Siemens scanners
      slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
    else
      if ~isfield(JOB,'slice_order')
        error('Slice timing information cannot be found in the DICOM header!')
      end
    end
    ref_slice_msec = JOB.TR_sec*1000/2;
  elseif isfield(JOB,'slice_order_msec')
    disp(['Found .slice_order_msec; applying those values']);
    slice_order_msec = JOB.slice_order_msec;
    ref_slice_msec = JOB.TR_sec*1000/2;
    JOB.noSTC = 0;
  else
    warning('NO slice timing information is given. Creating a link instead of STC...');
    JOB.noSTC=1;
  end
end
% override of manually given:
if isfield(JOB,'slice_order_msec')
  slice_order_msec = JOB.slice_order_msec;
end
if isfield(JOB,'ref_slice_msec')
  ref_slice_msec = JOB.ref_slice_msec;
end


%% 3. slice timing correction
% output:
% [1]  a${epi}.nii   : slice-timing corrected EPI
ls(fname_epi);
[p1,f1,e1] = myfileparts(fname_epi);
fname_output = [p1,'/a',f1,e1];
if ~exist(fname_output,'file') || JOB.overwrite
  disp('[3] slice timing correction..');
  disp(['> Slice order (msec) = ',num2str(reshape(slice_order_msec,1,[]))]);
  disp(['> Reference slice = ',num2str(ref_slice_msec),' msec']);
  if isfield(JOB,'noSTC') && JOB.noSTC
    disp(['Create a link for ',fname_epi,' as ',fname_output,'..']);
    unix(['ln -sf ',fname_epi,' ',fname_output])
  else
    if isfield(JOB,'stc_fsl')&&JOB.stc_fsl
      if isfield(JOB,'fname_dcm')
        hdr2 = spm_dicom_headers(JOB.fname_dcm);
        slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
        T=slice_order_msec'./(JOB.TR_sec*1000); % in TR
      else
        T=(JOB.slice_order_msec-1)./max(JOB.slice_order_msec);
      end
      fname_stc=[p1,'/',f1,'_timing.txt'];
      dlmwrite(fname_stc,T)
      unix(['FSLOUTPUTTYPE=NIFTI; slicetimer.fsl -i ',fname_epi,' -o ',fname_output,...
        ' --tcustom=',fname_stc,' -r ',num2str(JOB.TR_sec),' -v -d 3'])
    else
      st1.scans={};
      for t=1:NumFrames
        st1.scans{1}{t,1} = [fname_epi,',',num2str(t)];
      end
      st1.nslices = numel(slice_order_msec);
      st1.tr = JOB.TR_sec;
      if min(slice_order_msec)<1
        st1.ta = 0; %OR timing = [0 TR] when previous inputs are specified in milliseconds
      else
        st1.ta = JOB.TR_sec-(JOB.TR_sec/st1.nslices);
      end
      st1.so = slice_order_msec;
      st1.refslice = ref_slice_msec;
      st1.prefix = 'a';
      matlabbatch={};
      matlabbatch{1}.spm.temporal.st = st1;
      fname_matlabbatch=[p1,'/myspm_fmriprep3_',f1,'.mat'];
      save(fname_matlabbatch,'matlabbatch');
      spm_jobman('run', matlabbatch);
    end
    isdone(fname_out,3);
  end
end

%% 4a. preparing VDM
if isfield(JOB,'fname_mag') && isfield(JOB,'fname_pha')
  if ~isfield(JOB,'TEs_fmap')
    [p5,f5,~] = myfileparts(JOB.fname_mag);
    try
      json1 = readjson([p5,'/',f5,'.json']);
      [p5,f5,~] = myfileparts(JOB.fname_pha);
      json2 = readjson([p5,'/',f5,'.json']);
      JOB.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      %         JOB.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  JOB.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  [~,f_ph,e_ph] = myfileparts(JOB.fname_pha);
  JOB.fname_vdm = [p1,'/vdm5_sc',f_ph,e_ph];
  if ~exist(JOB.fname_vdm,'file')
    myspm_prepare_vdm(JOB.fname_mag, JOB.fname_pha, JOB.TEs_fmap, fname_epi, ...
      JOB.totalreadout_msec); % coreg of t1w here does not help much...
  end
  ls(JOB.fname_vdm)
end

%% 4. unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_a${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  a${epi}.mat      : [4x4xT] realign transform
% [3]  a${epi}_uw.mat   : unwarping meta data
% [4]  ua${epi}.nii     : unwarped/realigned image
% [5]  meanua${epi}.nii : mean image of [4]
realignunwarp1=[];
for t=1:NumFrames
  realignunwarp1.data.scans{t,1} = [p1,'/a',f1,e1,',',num2str(t)];
end
if isfield(JOB,'fname_vdm')
  ls(JOB.fname_vdm)
  realignunwarp1.data.pmscan = {[JOB.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(JOB,'realign_sep')
  realignunwarp1.eoptions.sep = JOB.realign_sep;
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
matlabbatch={};
matlabbatch{1}.spm.spatial.realignunwarp = realignunwarp1;
fname_in=[p1,'/a',f1,e1];
ls(fname_in);
fname_output = [p1,'/ua',f1,e1];
if ~exist(fname_output,'file') || JOB.overwrite
  disp('[4] Unwarp & realign..');
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  
  % visualize unwarping results:
  [p2,f2,e2] = myfileparts(JOB.fname_t1w);
  if strcmp(f2,'uni')
    fn_t1w = [p2,'/bm',f2,e2];
  else
    fn_t1w = [p2,'/',f2,e2];
  end
  compare_unwarped([p1,'/a',f1,e1], [p1,'/ua',f1,e1], fn_t1w)
  isdone(fname_out,4);
end
% %% 5. Intensity-bias correction of EPI (for stable coregistration)
fname_epi_unbiased = [p1,'/meanua',f1,e1]; % unprocessed EPI
% if ~exist(fname_epi_unbiased,'file')
%   unix(['mri_nu_correct.mni --i ',p1,'/meanua',f1,e1,' --o ',fname_epi_unbiased]);
% end
%% 6. Coregistration of EPI to native T1w (un-skull-stripped)
fname_t1w = [p2,'/m',f2,e2]; % bias-corrected T1w
% fname_epi_unbiased  = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
fname_spinecho = JOB.fname_spe;

job = [];
job.fname_moving = fname_epi_unbiased;
job.fname_fixed  = fname_spinecho;
job.reg_stages = 0;
job = myants_antsRegistration(job);
fn_reg_epi_to_spe = job.fname_reg;

job = [];
job.fname_moving = fname_spinecho;
job.fname_fixed  = fname_t1w;
job.reg_stages = 0;
job = myants_antsRegistration(job);
fn_reg_spe_to_t1w = job.fname_reg;
% %%
% fn_reg_epi_to_t1w = [p1,'/meanua',f1,'_to_m',f2,'_Composite.h5'];
% unix(['antsApplyTransforms -i ',fname_epi_unbiased,...
%   ' -r ',fname_t1w,...
%   ' -o [',fn_reg_epi_to_t1w,',1]']);
% 
%%
job=[];
job.fname_fixed = 'mni_brain.nii';
job.fname_moving = fname_epi_unbiased;
job.transforms = {
  fn_reg_epi_to_spe,0;
  fn_reg_spe_to_t1w,0;
  fn_reg_t1w_to_mni,0};
myants_antsApplyTransforms(job)

%% 7. Compcor
if ~isfield(JOB,'num_pcs'), JOB.num_pcs=3; end
if isfield(JOB,'restbpf'), JOB.bpf=[0.009 0.08]; end
if ~isfield(JOB,'bpf'), JOB.bpf=[0 inf]; end
job1=struct('path1',p1, 't1w_suffix',f2, 'bpf', JOB.bpf, ...
  'TR_sec', JOB.TR_sec, 'num_pcs', JOB.num_pcs, ...
  'name_epi', ['ua',f1,e1], 'name_t1w', ['bm',f2,e2], ...
  'dir_data', p1, 'name_rp',['rp_a',f1,'.txt']);
if isfield(JOB,'nofigure'), job1.nofigure=JOB.nofigure; end
cc_suffix= sprintf('n%db%0.2f-%0.2f',JOB.num_pcs, JOB.bpf);
fname_epi_ua = [p1,'/ua',f1,e1];
[~,name1,~]= myfileparts(fname_epi_ua);
fname_out= [p1,'/',name1,'_',cc_suffix,'_eigenvec.txt'];
if isfield(JOB,'cov_idx'), job1.cov_idx=JOB.cov_idx; end
if isfield(JOB,'dir_fig'), job1.dir_fig=JOB.dir_fig; end
if isfield(JOB,'out_prefix'), job1.out_prefix=JOB.out_prefix; end
if isfield(JOB,'dir_fs'), job1.dir_fs=JOB.dir_fs; end
if isfield(JOB,'subjid'), job1.subjid=JOB.subjid; end
if ~exist(fname_out,'file') || JOB.overwrite || isfield(JOB,'cov_idx')
  disp(['[7] Creating ',num2str(JOB.num_pcs),' anatomical CompCor regressors..']);
  job1 = myspm_denoise(job1);
  if isfield(job1,'fname_cov')
    fname_cov = job1.fname_cov;
    JOB.fname_cov = fname_cov;
    fname_cc = job1.fname_out;
  end
  isdone(fname_out,7);
end

%% 8. Apply forward deformation on EPI for registration into MNI152

fname_out=[p1,'/wua',f1,e1];
fn_mni_ref=[p1,'/mni_brain_ref.nii'];
if numel(JOB.vox_mm) == 1
  JOB.vox_mm = ones(1,3) * JOB.vox_mm;
end
unix(sprintf('mri_convert %s/mni_brain.nii -vs %i %i %i %s/mni_brain_ref.nii',...
  p1, JOB.vox_mm, p1))
if ~exist(fname_out,'file') || JOB.overwrite
  disp('[8] Applying antsRegistration to EPI..');
  dn=tempname;
  mkdir(dn);
  unix(['FSLOUTPUTTYPE=NIFTI; fslsplit ',fname_epi_ua,' ',dn,'/ -t']);
  FN=dir([dn,'/*.nii']);
  for v=1:numel(FN)
    unix(['antsApplyTransforms -d 3',...
      ' --input ',FN(v).folder,'/',FN(v).name,...
      ' --reference-image ',fn_mni_ref,...
      ' --output ',FN(v).folder,'/w',FN(v).name,...
      ' --float 1 ',...
      ' --interpolation BSpline[3] ',...
      ' --transform [',fn_reg_epi_to_spe,',0] ',...
      ' --transform [',fn_reg_spe_to_t1w,',0] ',...
      ' --transform [',fn_reg_t1w_to_mni,',0] ']);
  end
  unix(['FSLOUTPUTTYPE=NIFTI; cd ',dn,'; fslmerge -t ',fname_out,' w*.nii']);
  rmdir(dn,'s')
  isdone(fname_out,8);
  slices(fname_out,...
    struct('fname_png',[p1,'/wua',f1,'_in_mni152.png'],'contournum',1),...
    fn_mni_ref)
end

%% (9). Transforming denoised image:
if exist('fname_cov','var')
  [p3,f3,e3]=myfileparts(fname_cov);
  fname_out=[p3,'/w',f3,e3];
  if ~exist(fname_out,'file') || JOB.overwrite
    disp('[9] Registration of denoised EPI to MNI152..');
    dn=tempname;
    mkdir(dn);
    unix(['FSLOUTPUTTYPE=NIFTI; fslsplit ',fname_cov,' ',dn,'/ -t']);
    FN=dir([dn,'/*.nii']);
    for v=1:numel(FN)
      unix(['antsApplyTransforms -d 3',...
        ' --input ',FN(v).folder,'/',FN(v).name,...
        ' --reference-image ',fn_mni_ref,...
        ' --output ',FN(v).folder,'/w',FN(v).name,...
        ' --float 1 ',...
        ' --interpolation BSpline[3] ',...
      ' --transform [',fn_reg_epi_to_spe,',0] ',...
      ' --transform [',fn_reg_spe_to_t1w,',0] ',...
      ' --transform [',fn_reg_t1w_to_mni,',0] ']);
    end
    unix(['FSLOUTPUTTYPE=NIFTI; cd ',dn,'; fslmerge -t ',fname_out,' w*.nii']);
    rmdir(dn,'s')
    isdone(fname_out,9);
    slices(fname_out,...
      struct('fname_png',[p1,'/wua',f1,'_in_mni152.png'],'contournum',1),...
      fn_mni_ref)
  end
end

%% (10). Smoothing
if isfield(JOB,'fwhm_mm')
  myspm_smooth(struct('fname',[p1,'/wua',f1,e1],'fwhm_mm',JOB.fwhm_mm));
  if exist('fname_cov','var')
    [p3,f3,e3]=myfileparts(fname_cov);
    myspm_smooth(struct('fname',[p3,'/wua',f3,e3],'fwhm_mm',JOB.fwhm_mm));
  end
end

%   catch ME
%     warning(ME.identifier, '##### ERROR DURING PROCESSING: %s\n',JOB.fname_epi{j});
%     warning(ME.identifier, '%s', ME.message);
%     for iStack=1:numel(ME.stack)
%       fprintf('#####  In %s (line %i):\n', ...
%         ME.stack(iStack).file, ME.stack(iStack).line);
%     end
%     JOB.ME = ME;
%     JOB.error_with = JOB.fname_epi{j};
%   end
end

%% Subroutines:
function isdone(fname_out,proci)
if exist(fname_out,'file')
  disp(['[',num2str(proci),'] done: ',fname_out])
else
  error(['[',num2str(proci),'] failed: ',fname_out,' was not created'])
end
end
