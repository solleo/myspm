function JOB = myspm_fmriprep12_dartel (JOB)
% JOB = myspm_fmriprep12_dartel (JOB)
%
% For preprocessing of EPI and T1w images. This script is based on a batch
% script included in SPM12.
% It does slice-timing-correction before realignment/unwarping.
%
% This does:
%  (1) fname_epi -> slice timing correction -> unified unwarp+realign
%  (2) fname_t1w -> unified segmentation+normalization -> brain-masking (skull-stripping) -> normalize into mni
%  (3) fname_epi -> coregister into t1w (only header) -> normalize into mni (resampling) -> smoothing (default: fwhm=2.5 vox)
%
% JOB requires:
% [data]
%  .fname_epi  : a string for a single session; a cell for multiple sessions
%  .fname_t1w  : a string for anatomical image (full path)
% (.fname_dcm) : a string for a dicom file of any volume of the session
%                to read acquisition parameters TE, TR, readout time, ...
%
% [temporal process]
% (.TR_sec)    : it's better to manually enter..
% (.fname_vdm) : if a fieldmap is available (or .fname_mag, .fname_pha, .TEs_fmap)
%
% [Slice timing correction]
% (.slice_order)    : can be 1-based indices or actual time from the pulse (msec)
% (.ref_slice) : reference slice (in the same unit as .slice_order)
% (.fname_dcm) :
% (.noSTC)     : for a sparse sampling (e.g. TR>6) STC would make no effect
%                or could introduce a strange interpolation between distant time points
% [Misc.]
% (.dir_fig)   : a directory to save head motion plots of multiple subjects together
% (.vox_mm)    : resampling resolution in MNI152
%
% ["Denoising" or motion artifacts suppression]
% (.bpf)       : (default=[1/128 inf]) or ([0.009 0.08] from Satterthwaite.2013)
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
spm('Defaults','fmri');
spm_jobman('initcfg');

%% ---T1w processing---
%% 1. Unified segmentation
% outputs:
% [1]  ${t1w}_seg8.mat : deformation field in native space & meta
% [2]  y_${t1w}.nii    : cosine functions in template space (wtf?)
% [3]  c?${t1w}.nii    : tissue prob maps in native space
% [4]  mwc?${t1w}.nii  : JD modulated TPM in template space
% [5]  m${t1w}.nii     : bias-corrected t1w
% [6]  bm${t1w}.nii    : skull-stripped brain image
[p2,f2,e2]=myfileparts(JOB.fname_t1w);
JOB.fname_t1w=[p2,'/',f2,e2];
dir_tpm=[spm('dir'),'/tpm'];
matlabbatch={};
preproc=[];
preproc.channel.vols = {[JOB.fname_t1w,',1']};
preproc.channel.biasreg = 0.001;
preproc.channel.biasfwhm = 60;
preproc.channel.write = [0 1];
preproc.tissue(1).tpm = {[dir_tpm,'/TPM.nii,1']};
preproc.tissue(1).ngaus = 1;
preproc.tissue(1).native = [1 1]; % [native dartel]
preproc.tissue(1).warped = [0 1];
preproc.tissue(2).tpm = {[dir_tpm,'/TPM.nii,2']};
preproc.tissue(2).ngaus = 1;
preproc.tissue(2).native = [1 1]; % [native dartel]
preproc.tissue(2).warped = [0 1];
preproc.tissue(3).tpm = {[dir_tpm,'/TPM.nii,3']};
preproc.tissue(3).ngaus = 2;
preproc.tissue(3).native = [1 0];
preproc.tissue(3).warped = [0 1];
preproc.tissue(4).tpm = {[dir_tpm,'/TPM.nii,4']};
preproc.tissue(4).ngaus = 3;
preproc.tissue(4).native = [0 0];
preproc.tissue(4).warped = [0 0];
preproc.tissue(5).tpm = {[dir_tpm,'/TPM.nii,5']};
preproc.tissue(5).ngaus = 4;
preproc.tissue(5).native = [0 0];
preproc.tissue(5).warped = [0 0];
preproc.tissue(6).tpm = {[dir_tpm,'/TPM.nii,6']};
preproc.tissue(6).ngaus = 2;
preproc.tissue(6).native = [0 0];
preproc.tissue(6).warped = [0 0];
preproc.warp.mrf = 1;
preproc.warp.cleanup = 1;
preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
preproc.warp.affreg = 'mni';
if isfield(JOB,'isestern') && JOB.isestern
 preproc.warp.affreg = 'eastern';
end
preproc.warp.fwhm = 0;
preproc.warp.samp = 3;
preproc.warp.write = [0 1];
matlabbatch{1}.spm.spatial.preproc = preproc;
fname_in = JOB.fname_t1w;
ls(fname_in);
fname_out = [p2,'/rc1',f2,e2];
if ~exist(fname_out,'file') || JOB.overwrite
 disp('[1] Unified segmentation of T1w..');
 fname_matlabbatch=[p2,'/spm12_fmriprep1.mat'];
 save(fname_matlabbatch,'matlabbatch');
 spm_jobman('run', matlabbatch);
 nii=cell(1,3);
 for c=1:3
  nii{c} = load_untouch_nii([p2,'/c',num2str(c),f2,e2]);
 end
 manat = load_untouch_nii([p2,'/m',f2,e2]);
 manat.img(~(nii{1}.img | nii{2}.img | nii{3}.img))=0;
 save_untouch_nii(manat, fname_out);
 clear manat nii
 isdone(fname_out,1);
end

%% 2. register the native brain to MNI152 using DARTEL
% output:
% [1]  wb${t1w}.nii : transformed skull-stripped brain in template space
fname_in = [p2,'/bm',f2,e2];
ls(fname_in)
fname_def = [p2,'/y_',f2,e2];
ls(fname_def)
fname_out = [p2,'/wbm',f2,e2];
matlabbatch={};
write1=[];
write1.subj.resample = {fname_in};
write1.subj.def = {fname_def};
write1.woptions.bb = [-78 -112 -70; 78 76 85];
write1.woptions.vox = [1 1 1];
write1.woptions.interp = 4;
matlabbatch{1}.spm.spatial.normalise.write = write1;
if ~exist(fname_out,'file') || JOB.overwrite
 disp('[2] Registering T1w to MNI using DARTEL..');
 fname_matlabbatch=[p2,'/spm12_fmriprep2.mat'];
 save(fname_matlabbatch,'matlabbatch');
 spm_jobman('run', matlabbatch);
 myspm_denan(fname_out,1);
 isdone(fname_out,2);
end
return

%% ---EPI processing---
% Find filenames
if ~iscell(JOB.fname_epi)
 JOB.fname_epi = {JOB.fname_epi};
end
for j = 1:numel(JOB.fname_epi)
 % Find parameters of the EPI scan
 fname_epi = JOB.fname_epi{j};
 [path1,name1,ext1] = myfileparts(fname_epi);
 fname_epi = [path1,'/',name1,ext1];
 text1=['[Session ',num2str(j),']:'];
 disp([text1,repmat('-',1,60-numel(text1))]);
 ls(fname_epi)
 hdr = load_nii_hdr(fname_epi);
 NumFrames = hdr.dime.dim(5);
 disp(['> Number of frames = ',num2str(NumFrames)]);
 if ~isfield(JOB,'TR_sec')
  JOB.TR_sec=hdr.dime.pixdim(5);
  if ~JOB.TR_sec
   error(['TR=0 in the file header; Enter JOB.TR_sec !']);
  end
 end
 disp(['> TR = ',num2str(JOB.TR_sec),' sec']);
 if JOB.TR_sec >= 6 && ~isfield(JOB,'noSTC') && ~isfield(JOB,'fname_dcm')
  warning(['TR is long (>= 6 s) and found no dicom file (.fname_dcm) to read actual slice time. Thus setting .noSTC=1']);
  disp(['[!] It may possible to use vectors in .slice_order and .ref_slice']);
  JOB.noSTC=1;
 end
 % find slice timing in msec and repetition time in sec from a example DICOM
 if isfield(JOB,'fname_dcm')
  hdr2 = spm_dicom_headers(JOB.fname_dcm);
  if isfield(hdr2{1},'Private_0019_1029') % recent Siemens scanners
   slice_order = hdr2{1}.Private_0019_1029; % this could be also time
  else
   if ~isfield(JOB,'slice_order')
    error('Slice timing information cannot be found in the DICOM header!')
   end
  end
  JOB.TR_sec  = hdr2{1}.RepetitionTime/1000;
  ref_slice = JOB.TR_sec/2*1000;  % in msec
  NumFrames = hdr.dime.dim(5);
 elseif isfield(JOB,'slice_order') && isfield(JOB,'ref_slice')
  disp(['Found .slice_order and .ref_slice; applying those values']);
  ref_slice = JOB.ref_slice;
  slice_order    = JOB.slice_order;
  JOB.noSTC = 0;
 else
  warning('NO slice timing information is given. Creating a link instead of STC...');
  JOB.noSTC=1;
 end
 
 %% 3. slice timing correction
 % output:
 % [1]  a${epi}.nii   : slice-timing corrected EPI
 ls(fname_epi);
 [p1,f1,e1] = myfileparts(fname_epi);
 fname_output = [p1,'/a',f1,e1];
 if ~exist(fname_output,'file') || JOB.overwrite
  disp('[3] slice timing correction..');
  if isfield(JOB,'noSTC') && JOB.noSTC
   disp(['Create a link for ',fname_epi,' as ',fname_output,'..']);
   unix(['ln -sf ',fname_epi,' ',fname_output])
  else
   st1.scans={};
   for t=1:NumFrames
    st1.scans{1}{t,1} = [fname_epi,',',num2str(t)];
   end
   st1.nslices = numel(slice_order);
   st1.tr = JOB.TR_sec;
   if min(slice_order)<1
    st1.ta = 0; %OR timing = [0 TR] when previous inputs are specified in milliseconds
   else
    st1.ta = JOB.TR_sec-(JOB.TR_sec/st1.nslices); 
   end
   st1.so = slice_order;
   st1.refslice = ref_slice;
   st1.prefix = 'a';
   matlabbatch={};
   matlabbatch{1}.spm.temporal.st = st1;
   fname_matlabbatch=[p1,'/spm12_fmriprep3.mat'];
   save(fname_matlabbatch,'matlabbatch');
   spm_jobman('run', matlabbatch);
  end
  isdone(fname_out,3);
 end
 
 %% 4a. preparing VDM
 if isfield(JOB,'fname_mag') && isfield(JOB,'fname_pha') ...
   && isfield(JOB,'TEs_fmap') && isfield(JOB,'totalreadout_msec')
  if ~isfield(JOB,'iPAT'), JOB.iPAT=1; end
  [~,f_ph,e_ph]=myfileparts(JOB.fname_pha);
  JOB.fname_vdm=[p1,'/vdm5_sc',f_ph,e_ph];
  if ~exist(JOB.fname_vdm,'file')
   myspm_prepare_vdm(JOB.fname_mag, JOB.fname_pha, JOB.TEs_fmap, fname_epi, ...
    JOB.totalreadout_msec, JOB.iPAT, JOB.fname_t1w);
  end
  ls(JOB.fname_vdm)
 end
 
 %% 4. unwarp+realign to MEAN IMAGE
 % outputs:
 % [1]  rp_a${epi}.txt   : six rigid-body motion parameters
 % [2]  a${epi}.mat      : [4x4xT] realign transform
 % [3]  a${epi}_uw.mat   : unwarping meta data
 % [4]  ua${epi}.nii     : unwarped/realigned image
 % [5]  meanua${epi}.nii : mean image of [4]
 realignunwarp1=[];
 for t=1:NumFrames
  realignunwarp1.data.scans{t,1} = [p1,'/a',f1,e1,',',num2str(t)];
 end
 if isfield(JOB,'fname_vdm')
  realignunwarp1.data.pmscan = {[JOB.fname_vdm,',1']};
 else
  realignunwarp1.data.pmscan = {''};
 end
 realignunwarp1.eoptions.quality = 0.9;
 realignunwarp1.eoptions.sep = 4;
 realignunwarp1.eoptions.fwhm = 5;
 realignunwarp1.eoptions.rtm = 1; % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
 realignunwarp1.eoptions.einterp = 4;
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
  fname_matlabbatch=[p1,'/spm12_fmriprep4.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  myspm_denan(fname_out,1);
  isdone(fname_out,4);
 end
  
 %% 5. Intensity-bias correction of EPI (for better coregistration)
 fname_epi_unbiased = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
 if ~exist(fname_epi_unbiased,'file')
  unix(['mri_nu_correct.mni --i ',p1,'/meanua',f1,e1,' --o ',fname_epi_unbiased]);
 end
 
 %% 6. Coregistration of EPI to native T1w
 % outputs: <modifying transform matrices in headers>
 % [1]  mmeanua${epi}.nii
 % [2]  ua${epi}.nii        : "other" images, header modified
 % [3]  ua${epi}.mat        : "other" images, header modified
 estimate1=[];
 fname_t1w_brain = [p2,'/bm',f2,e2]; % bias-corrected T1w
 estimate1.ref{1}    = fname_t1w_brain;
 ls(fname_t1w_brain)
 fname_epi_unbiased = [p1,'/mmeanua',f1,e1]; % bias-corrected EPI
 estimate1.source{1} = fname_epi_unbiased ;
 ls(fname_epi_unbiased )
 % ------------------here you specified OTHER imgaes-----------------
 fname_epi_ua = [p1,'/ua',f1,e1];
 ls(fname_epi_ua);
 estimate1.other{1} = fname_epi_ua; % NIFTI-format have one xfm for all volumes
 % ------------------------------------------------------------------
 estimate1.eoptions.cost_fun = 'nmi';
 estimate1.eoptions.sep = [4 2];
 estimate1.eoptions.tol = ...
  [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
 estimate1.eoptions.fwhm = [7 7];
 matlabbatch={};
 matlabbatch{1}.spm.spatial.coreg.estimate = estimate1;
 fname_out=[p1,'/ua',f1,'.mat'];
 if ~exist(fname_out,'file') || JOB.overwrite
  disp('[6] Coregistration of EPI to native T1w..');
  fname_matlabbatch=[p1,'/spm12_fmriprep6.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  isdone(fname_out,6);
 end
 %% resample the first volume (to check coregistration quality)
 [p4,f4,e4]=myfileparts(fname_epi_unbiased);
 if ~exist([p4,'/r',f4,e4],'file')  || JOB.overwrite
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
  unix(['slices ',fname_epi_in_t1w,' ',fname_t1w_brain, ...
   ' -o ',p1,'/rmmeanua',f1,'_in_',f2,'.gif'])
 end
 %% 7. Compcor
 if isfield(JOB,'NOCOMPCOR') && JOB.NOCOMPCOR
  disp('[7] Skipping CompCor.');
 else
  if ~isfield(JOB,'num_pcs'), JOB.num_pcs=3; end
  if isfield(JOB,'restbpf'), JOB.bpf=[0.009 0.08]; end
  if ~isfield(JOB,'bpf'), JOB.bpf=[0 inf]; end
  job1=struct('path1',p1, 't1w_suffix',f2, 'bpf', JOB.bpf, ...
   'TR_sec', JOB.TR_sec, 'num_pcs', JOB.num_pcs, ...
   'name_epi', ['ua',f1,e1], 'name_t1w', ['bm',f2,e2], ...
   'dir_data', p1, 'name_rp',['rp_a',f1,'.txt']);
  cc_suffix= sprintf('n%db%0.2f-%0.2f',JOB.num_pcs, JOB.bpf);
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
    JOB.fname_cov = fname_cov ;
    fname_cc = job1.fname_out;
   end
   isdone(fname_out,7);
  end
 end
 %% 8. Apply forward deformation on EPI for registration into MNI152
 matlabbatch={};
 write1=[];
 write1.subj.def{1} = fname_def;
 write1.subj.resample{1} = fname_epi_ua;
 write1.woptions.bb = [-78 -112 -70; 78 76 85];
 if ~isfield(JOB,'vox_mm')
  JOB.vox_mm = hdr.dime.pixdim(2:4);
  disp(['Resampling at [',num2str(JOB.vox_mm),'] mm in MNI152 space..'])
 end
 if numel(JOB.vox_mm) == 1
  JOB.vox_mm=[1 1 1]*JOB.vox_mm;
 end
 write1.woptions.vox = JOB.vox_mm;
 write1.woptions.interp = 4;
 matlabbatch{1}.spm.spatial.normalise.write = write1;
 fname_out=[p1,'/wua',f1,e1];
 if ~exist(fname_out,'file') || JOB.overwrite
  disp('[8] Registration of EPI to MNI152..');
  fname_matlabbatch=[p1,'/spm12_fmriprep8.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  myspm_denan(fname_out,1);
  isdone(fname_out,8);
 end
 fname_gif=[p1,'/wua',f1,'_in_mni152.gif'];
 fname_epi_in_mni1=[p1,'/wua',f1,'1.nii.gz'];
 setenv('FSLOUTPUTTYPE','NIFTI_GZ');
 unix(['fslroi ',p1,'/wua',f1,e1,' ',fname_epi_in_mni1,' 0 1']);
 myspm_denan(fname_epi_in_mni1);
 ls(fname_epi_in_mni1);
 fname_mni=[getenv('FSLDIR'),'/data/standard/MNI152_T1_2mm.nii.gz'];
 fname_mni1=[p1,'/mni152.nii.gz'];
 unix(['mri_convert --like ',fname_epi_in_mni1,' ',fname_mni,' ',fname_mni1]);
 ls(fname_mni1)
 unix(['slices ',fname_epi_in_mni1,' ',fname_mni1,' -o ',fname_gif]);
 ls(fname_gif)
 
 %% (9). Transforming residual image:
 if exist('fname_cov','var')
  matlabbatch={};
  write1=[];
  write1.subj.def{1} = fname_def;
  write1.subj.resample{1} = fname_cov;
  write1.woptions.bb = [-78 -112 -70; 78 76 85];
  write1.woptions.vox = JOB.vox_mm;
  write1.woptions.interp = 4;
  matlabbatch{1}.spm.spatial.normalise.write = write1;
  [p3,f3,e3]=myfileparts(fname_cov);
  fname_out=[p3,'/w',f3,e3];
  if ~exist(fname_out,'file') || JOB.overwrite
   disp('[9] Registration of denoised EPI to MNI152..');
   fname_matlabbatch=[p1,'/spm12_fmriprep9.mat'];
   save(fname_matlabbatch,'matlabbatch');
   spm_jobman('run', matlabbatch);
   nii0=load_untouch_nii(fname_cov);
   bgval=nii0.img(1);
   nii=load_untouch_nii(fname_out);
   nii.img(isnan(nii.img))=bgval;
   save_untouch_nii(nii,fname_out);
   myspm_denan(fname_out,1);
   isdone(fname_out,9);
  end
 end
 
end % for each session
end

%% Subroutines:
function isdone(fname_out,proci)
if exist(fname_out,'file')
 disp(['[',num2str(proci),'] done: ',fname_out])
else
 error(['[',num2str(proci),'] failed: ',fname_out,' was not created'])
end
end
