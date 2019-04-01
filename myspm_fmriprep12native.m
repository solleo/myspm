function EXP = myspm_fmriprep12native (EXP)
% EXP = myspm_fmriprep12native (EXP)
%
% NOTE: for analysis in native space, register EPI into t1w
% For preprocessing of EPI and T1w images. This script is based on a batch
% script included in SPM12.
% It does slice-timing-correction before realignment/unwarping.
%
% This does:
%  (1) fname_epi -> slice timing correction -> unified unwarp+realign
%  (2) fname_t1w -> unified segmentation+normalization -> brain-masking (skull-stripping) -> normalize into mni
%  (3) fname_epi -> coregister into t1w (resampled) -> smoothing (default: fwhm=2.5 vox)
%
% EXP requires:
% [data]
% (.dname_data): not to specify for all data files
%  .fname_epi  : a string for a single session
% (.fname_json)
%  .fname_t1w  : a string for anatomical image (full path)
%
% [slice timing correction]
% (.slice_order_msec) : can be 1-based indices or actual time from the pulse (msec)
%                  (can be read from a paired json file) *WARNING: (sec)!*
% (.ref_slice_msec)   : reference slice (in the same unit as .slice_order_msec)
% (.fname_dcm)   : DICOM file of a volume to read slice time (only for
%                  SIEMENS scanners)
% (.stc_fsl)     : use FSL instead of SPM for a very big file
% (.noSTC)       : for a sparse sampling (e.g. TR>6) STC would make no effect
%                  or could introduce a strange interpolation between distant time points
%
% [temporal process]
% (.TR_sec)    : can be read from a paired .json file.
%                (can be read from a paired json file)
%
% [geometric distortion correction]
% (.fname_vdm) : if a fieldmap is available
%                (or .fname_mag, .fname_pha, (.TEs_fmap), and (.totalreadout_msec))
%
% [Misc.]
% (.dir_fig)   : a directory to save head motion plots of multiple subjects together
% (.vox_mm)    : resampling resolution in MNI152
%
% ["Denoising" or motion artifacts suppression]
% (.bpf)       : (default=[1/128 inf] Hz)
% (.restbpf)   : 1=[0.009 0.08] Hz from Satterthwaite.2013
% (.num_pcs)   : 3 [default] # of CompCor regressors to extract
% (.subjid)    : for motion plots and compcor if FS directory is available
% (.cov_idx)   : for motion-artifact removal: rp+cc+sc+gs
% (.out_prefix)
%
% (.overwrite) : 0 [default] or 1
%
% This script uses original and modified MATLAB codes by others:
% - NIFTI toolbox by Jimmy SHEN:
% - DPARSF by Chao-Gan YAN: http://rfmri.org/DPARSF
% - SCREEN2JPEG by Sean P. McCARTHY (MALTAB builtin)
%
% (cc) 2015-2019. sgKIM. mailto:solleo@gmail.com  https://ggooo.wordpress.com

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
if ~isfield(EXP,'overwrite'), EXP.overwrite =0; end
if isfield(EXP,'dname_data'), cd(EXP.dname_data); end
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
[p2,f2,e2]=myfileparts(EXP.fname_t1w);
EXP.fname_t1w=[p2,'/',f2,e2];
fname_in = EXP.fname_t1w;
ls(fname_in);
fname_out = [p2,'/bm',f2,e2];
fname_def = [p2,'/y_',f2,e2];
if ~exist(fname_out,'file') || EXP.overwrite
  disp('[1] Unified segmentation of T1w..');
  exp1 = EXP;
  exp1.norm = 1;
  myspm_seg12(exp1,'ss')
  isdone(fname_out,1);
end

%% ---EPI processing---
if ~isfield(EXP,'fname_json')
  [path1,name1,~] = myfileparts(EXP.fname_epi);
  EXP.fname_json = [path1,'/',name1,'.json'];
end
fname_epi = EXP.fname_epi;
[path1,name1,ext1] = myfileparts(fname_epi);
fname_epi = [path1,'/',name1,ext1];
ls(fname_epi)
hdr = load_nii_hdr(fname_epi);
NumFrames = hdr.dime.dim(5);

% using json files *** it's SEC!! ***
fn_json = EXP.fname_json;
if exist(fn_json,'file')
  disp(['Found json file: ',fn_json]);
  disp(['Reading TR & SliceTiming from json']);
  json = jsondecode(fileread(fn_json));
  EXP.TR_sec = json.RepetitionTime;
  slice_order_msec = json.SliceTiming*1000;
  ref_slice_msec = EXP.TR_sec*1000/2;
else
  % if no json
  if ~isfield(EXP,'TR_sec')
    EXP.TR_sec=hdr.dime.pixdim(5);
    if ~EXP.TR_sec
      error(['TR=0 in the file header; Enter EXP.TR_sec !']);
    end
  end
  if EXP.TR_sec >= 6 && ~isfield(EXP,'noSTC') && ~isfield(EXP,'fname_dcm')
    warning(['TR is long (>= 6 s) and found no dicom file (.fname_dcm) to read actual slice time. Thus setting .noSTC=1']);
    disp(['[!] It may possible to use vectors in .slice_order_msec and .ref_slice_msec']);
    EXP.noSTC=1;
  end
  % find slice timing in msec and repetition time in sec from a example DICOM
  if isfield(EXP,'fname_dcm')
    hdr2 = spm_dicom_headers(EXP.fname_dcm);
    if isfield(hdr2{1},'Private_0019_1029') % recent Siemens scanners
      slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
    else
      if ~isfield(EXP,'slice_order')
        error('Slice timing information cannot be found in the DICOM header!')
      end
    end
    EXP.TR_sec  = hdr2{1}.RepetitionTime/1000;
    ref_slice_msec = EXP.TR_sec*1000/2;
    NumFrames = hdr.dime.dim(5);
  elseif isfield(EXP,'slice_order_msec') && isfield(EXP,'ref_slice_msec')
    disp(['Found .slice_order_msec and .ref_slice_msec; applying those values']);
    slice_order_msec = EXP.slice_order_msec;
    ref_slice_msec = EXP.ref_slice_msec;
    EXP.noSTC = 0;
  else
    warning('NO slice timing information is given. Creating a link instead of STC...');
    EXP.noSTC=1;
  end
end
if isfield(EXP,'slice_order_msec')
  slice_order_msec = EXP.slice_order_msec;
end
if isfield(EXP,'ref_slice_msec')
  ref_slice_msec = EXP.ref_slice_msec;
end

disp(['> Number of frames = ',num2str(NumFrames)]);
disp(['> TR = ',num2str(EXP.TR_sec),' sec']);

%% 3. slice timing correction
% output:
% [1]  a${epi}.nii   : slice-timing corrected EPI
ls(fname_epi);
[p1,f1,e1] = myfileparts(fname_epi);
fname_output = [p1,'/a',f1,e1];
if ~exist(fname_output,'file') || EXP.overwrite
  disp('[3] slice timing correction..');
  disp(['> Slice order (msec) = ',num2str(reshape(slice_order_msec,1,[])),' msec']);
  disp(['> Reference slice = ',num2str(ref_slice_msec),' msec']);
  if isfield(EXP,'noSTC') && EXP.noSTC
    disp(['Create a link for ',fname_epi,' as ',fname_output,'..']);
    unix(['ln -sf ',fname_epi,' ',fname_output])
  else
    if isfield(EXP,'stc_fsl')&&EXP.stc_fsl
      if isfield(EXP,'fname_dcm')
        hdr2 = spm_dicom_headers(EXP.fname_dcm);
        slice_order_msec = hdr2{1}.Private_0019_1029; % in msec
        T=slice_order_msec'./(EXP.TR_sec*1000); % in TR
      else
        T=(EXP.slice_order_msec-1)./max(EXP.slice_order_msec);
      end
      fname_stc=[p1,'/',f1,'_timing.txt'];
      dlmwrite(fname_stc,T)
      unix(['FSLOUTPUTTYPE=NIFTI; slicetimer.fsl -i ',fname_epi,' -o ',fname_output,...
        ' --tcustom=',fname_stc,' -r ',num2str(EXP.TR_sec),' -v -d 3'])
    else
      st1.scans={};
      for t=1:NumFrames
        st1.scans{1}{t,1} = [fname_epi,',',num2str(t)];
      end
      st1.nslices = numel(slice_order_msec);
      st1.tr = EXP.TR_sec;
      if min(slice_order_msec)<1
        st1.ta = 0; %OR timing = [0 TR] when previous inputs are specified in milliseconds
      else
        st1.ta = EXP.TR_sec-(EXP.TR_sec/st1.nslices);
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
if isfield(EXP,'fname_mag') && isfield(EXP,'fname_pha')
  if ~isfield(EXP,'TEs_fmap')
    [p5,f5,~] = myfileparts(EXP.fname_mag);
    try
      json1 = readjson([p5,'/',f5,'.json']);
      [p5,f5,~] = myfileparts(EXP.fname_pha);
      json2 = readjson([p5,'/',f5,'.json']);
      EXP.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      %         EXP.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  EXP.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  [~,f_ph,e_ph] = myfileparts(EXP.fname_pha);
  EXP.fname_vdm = [p1,'/vdm5_sc',f_ph,e_ph];
  if ~exist(EXP.fname_vdm,'file')
    myspm_prepare_vdm(EXP.fname_mag, EXP.fname_pha, EXP.TEs_fmap, fname_epi, ...
      EXP.totalreadout_msec); % coreg of t1w here does not help much...
  end
  ls(EXP.fname_vdm)
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
if isfield(EXP,'fname_vdm')
  ls(EXP.fname_vdm)
  realignunwarp1.data.pmscan = {[EXP.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(EXP,'realign_sep')
  realignunwarp1.eoptions.sep = EXP.realign_sep;
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
if ~exist(fname_output,'file') || EXP.overwrite
  disp('[4] Unwarp & realign..');
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  
  % visualize unwarping results:
  [p2,f2,e2] = myfileparts(EXP.fname_t1w);
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
if ~exist(fname_out,'file') || EXP.overwrite
  disp('[6] Coregistration of EPI to native T1w..');
  fname_matlabbatch=[p1,'/myspm_fmriprep6_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  isdone(fname_out,6);
end
%% resample the first volume (to check coregistration quality)
[p4,f4,e4]=myfileparts(fname_epi_unbiased);
if ~exist([p4,'/r',f4,e4],'file')  || EXP.overwrite
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.write.ref{1}    = fname_t1w_brain;
  matlabbatch{1}.spm.spatial.coreg.write.source{1} = [fname_epi_unbiased,',1'];
  matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
  spm_jobman('run', matlabbatch);
  fname_epi_in_t1w=[p1,'/rmmeanua',f1,e1];
  myspm_denan(fname_epi_in_t1w,1);
  slices(fname_epi_in_t1w,...
    struct('fname_png',[p1,'/rmmeanua',f1,'_in_',f2,'.png'],'contournum',1),...
    fname_t1w_brain)
end
%% and resample all the images:
if ~exist([p1,'/r',f1,e1],'file')  || EXP.overwrite
  if ~isfield(EXP,'vox_mm')
    EXP.vox_mm = hdr.dime.pixdim(2:4);
    disp(['Resampling at [',num2str(EXP.vox_mm),'] mm in native T1w space..'])
  end
  fname_t1w_brain_ref=[p2,'/bm',f2,'_',num2str(EXP.vox_mm(1)),'mm',e2];
  unix(['mri_convert -vs ',num2str(EXP.vox_mm),' ',fname_t1w_brain,' ',fname_t1w_brain_ref])
  
  job = [];
  job.fname_ref = fname_t1w_brain_ref;
	job.fname_src = fname_epi_ua; 
	job.interp = 4;
  job.prefix = 'r';
  job = myspm_reslice(job);
end


%% 7. Compcor
if ~isfield(EXP,'num_pcs'), EXP.num_pcs=3; end
if isfield(EXP,'restbpf'), EXP.bpf=[0.009 0.08]; end
if ~isfield(EXP,'bpf'), EXP.bpf=[0 inf]; end
exp1=struct('path1',p1, 't1w_suffix',f2, 'bpf', EXP.bpf, ...
  'TR_sec', EXP.TR_sec, 'num_pcs', EXP.num_pcs, ...
  'name_epi', ['ua',f1,e1], 'name_t1w', ['bm',f2,e2], ...
  'dir_data', p1, 'name_rp',['rp_a',f1,'.txt']);
if isfield(EXP,'nofigure'), exp1.nofigure=EXP.nofigure; end
cc_suffix= sprintf('n%db%0.2f-%0.2f',EXP.num_pcs, EXP.bpf);
[~,name1,~]= myfileparts(fname_epi_ua);
fname_out= [p1,'/',name1,'_',cc_suffix,'_eigenvec.txt'];
if isfield(EXP,'cov_idx'), exp1.cov_idx=EXP.cov_idx; end
if isfield(EXP,'dir_fig'), exp1.dir_fig=EXP.dir_fig; end
if isfield(EXP,'out_prefix'), exp1.out_prefix=EXP.out_prefix; end
if isfield(EXP,'dir_fs'), exp1.dir_fs=EXP.dir_fs; end
if isfield(EXP,'subjid'), exp1.subjid=EXP.subjid; end
if ~exist(fname_out,'file') || EXP.overwrite || isfield(EXP,'cov_idx')
  disp(['[7] Creating ',num2str(EXP.num_pcs),' anatomical CompCor regressors..']);
  exp1 = myspm_denoise(exp1);
  if isfield(exp1,'fname_cov')
    fname_cov = exp1.fname_cov;
    EXP.fname_cov = fname_cov;
    fname_cc = exp1.fname_out;
  end
  isdone(fname_out,7);
end


%% (8). Smoothing
if isfield(EXP,'fwhm_mm')
  myspm_smooth(struct('fname',[p1,'/rua',f1,e1],'fwhm_mm',EXP.fwhm_mm));
  if exist('fname_cov','var')
    [p3,f3,e3]=myfileparts(fname_cov);
    myspm_smooth(struct('fname',[p3,'/rua',f3,e3],'fwhm_mm',EXP.fwhm_mm));
  end
end

%   catch ME
%     warning(ME.identifier, '##### ERROR DURING PROCESSING: %s\n',EXP.fname_epi{j});
%     warning(ME.identifier, '%s', ME.message);
%     for iStack=1:numel(ME.stack)
%       fprintf('#####  In %s (line %i):\n', ...
%         ME.stack(iStack).file, ME.stack(iStack).line);
%     end
%     EXP.ME = ME;
%     EXP.error_with = EXP.fname_epi{j};
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
