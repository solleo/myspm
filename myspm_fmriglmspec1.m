function EXP=myspm_fmriglmspec1 (EXP)
% EXP=myspm_fmriglmspec1 (EXP)
%
% This script helps you to run GLMs for the 1st-level fMRI analysis using fMRI model
%   specification.
% The important difference
%   from the simple GLM is that it deals with temporal aspects of fMRI data including
%   slow hedonamics, pre-whitening for auto-correlation, as well as high-pass
%   filtering for detrending.
%
% For one subject, one run
%  EXP.subjid
%  EXP.ridx
%
% Inputs:
% EXP requires fields for myspm_fmriglm.m:
% -input files
%  .files_query  'Nx1' for query to find image filenames, "${subjid}" will
%                be replaced by given subjIDs
% or
%  .filenames    {Nimgx1} (instead of files_query)
%
% -output directory
%  .dir_glm      'Nx1' directory to save SPM results
% or
%  .dir_base     'Nx1' directory for a subdirectory that has SPM.mat
% (.dir_prefix)
%
% - filenames (*.ons | *.param) for fMRI design matrix
%  .COND(k).name
%  .COND(k).onset
% (.COND(k).param)
% or
%  .name_cond {kx1}   only fname1 in [path1,'/',fname1,ext1]
%  .name_param {jx1}
%
%  .StimDurSec (kx1)
%  .TR
%  .reg
%  .fname_rp
% (.hpfcutoff)   : max([EXP.TR*EXP.NumFrames(j)/10 128])
%
% (cc) 2015. sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com/


%% 0. parsing & check inputs
global overwrite; if isempty(overwrite), overwrite=0; end
spm('Defaults','fmri');
spm_jobman('initcfg');

if isfield(EXP,'model')
  model_desc = fsss_model_desc(EXP.model);
else
  model_desc = EXP.model_desc;
end
if ~isfield(EXP,'dir_prefix'), EXP.dir_prefix=''; end

if ~isfield(EXP,'dir_glm')
  if ~isfield(EXP,'dir_base')
    EXP.dir_glm = [pwd,'/glm_',EXP.dir_prefix,model_desc];
  else
    EXP.dir_glm = [EXP.dir_base,'/glm_',EXP.dir_prefix,model_desc];
  end
end
[~,~]=mkdir(EXP.dir_glm);

% find epi fnames...
if isfield(EXP,'files_query')
  
  idx = strfind(EXP.files_query,'${subj}');
  prefix = EXP.files_query(1:idx-1);
  suffix = EXP.files_query(idx+7:end);
  EXP.files_query = [prefix,EXP.subjid,suffix];
  
  idx = strfind(EXP.files_query,'${ridx}');
  prefix = EXP.files_query(1:idx-1);
  suffix = EXP.files_query(idx+7:end);
  fnames{1} = [prefix,num2str(EXP.ridx),suffix];
else
  error('You need to specify inputs in EXP.files_query or EXP.filenames');
end
ls(fnames{1});
[path1,~,~] = fileparts_gz(fnames{1});
hdr = load_untouch_header_only(fnames{1});
EXP.NumFrames = hdr.dime.dim(5);
if ~isfield(EXP,'COND')
  if isfield(EXP,'name_cond')
    EXP.NumCond = numel(EXP.name_cond);
    readOnsetFile=1;
  else
    EXP.NumCond = 0;
  end
else
  EXP.NumCond = numel(EXP.COND);
  readOnsetFile=0;
end

% find log files...
NumCond = numel(EXP.name_cond);
for k=1:NumCond
  fnamelog = [path1,'/',EXP.name_cond{k}];
  idx = strfind(fnamelog,'${ridx}');
  prefix = fnamelog(1:idx-1);
  suffix = fnamelog(idx+7:end);
  fnamelog = [prefix,num2str(EXP.ridx),suffix];
  EXP.prefix_cond{k} = [prefix,num2str(EXP.ridx),suffix];
  [~,EXP.namestr_cond{k}] = fileparts(EXP.prefix_cond{k});
  
  if exist([fnamelog,'.ons'],'file')
    src = [fnamelog,'.ons'];
    [~,f1,e1] = fileparts(src);
    trg = [EXP.dir_base,'/',f1,e1];
    copyfile(src, trg);
  end
  if exist([fnamelog,'.par'],'file')
    src = [fnamelog,'.par'];
    [~,f1,e1] = fileparts(src);
    trg = [EXP.dir_base,'/',f1,e1];
    copyfile(src, trg);
  end
end

%% 1. model specification
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};     % for spm12
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; % for spm12
matlabbatch{1}.spm.stats.fmri_spec.dir = {EXP.dir_glm};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % middle point

sess.scans=cell(EXP.NumFrames,1);
for t=1:EXP.NumFrames
  sess.scans{t,1} = [fnames{1},',',num2str(t)]; % this must be <Tx1>
end

for k=1:EXP.NumCond
  if readOnsetFile  % to create EXP.COND
    EXP.COND(k).name  = EXP.namestr_cond{k};
    EXP.COND(k).onset = dlmread([EXP.prefix_cond{k},'.ons']);
  end
  sess.cond(k).name     = EXP.COND(k).name;
  sess.cond(k).onset    = EXP.COND(k).onset;
  if numel(EXP.StimDurSec)>1
    sess.cond(k).duration = EXP.StimDurSec(k);
  else
    sess.cond(k).duration = EXP.StimDurSec;
  end
  sess.cond(k).tmod = 0; % time modulation
  
  % now for parametric design
  isParamFile = isfield(EXP,'param_name');
  isParamCOND = isfield(EXP.COND(k),'param');
  if isParamFile || isParamCOND
    if isParamFile
      for i=1:size(EXP.param_name,1);
        if isfield(EXP,'param_name')
          EXP.COND(k).pmod(i).name  = EXP.param_name{i,k};
        else
          EXP.COND(k).pmod(i).name  = 'Modulation parameter';
        end
        EXP.COND(k).pmod(i).param = dlmread([dir0,'/',EXP.param_name{i,k},'.par']);
      end
      for i=1:size(EXP.param_name)
        sess.cond(k).pmod(i).name  = EXP.COND(k).pmod(i).name;
        sess.cond(k).pmod(i).param = EXP.COND(k).pmod(i).param;
        if ~isfield(EXP.COND(k).pmod,'poly')
          sess.cond(k).pmod(i).poly  = 1;
        else
          sess.cond(k).pmod(i).poly  = EXP.COND(k).pmod(i).poly;
        end
      end
    else
      sess.cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
    end
  end
end % loop for conditions

if ~isfield(EXP,'hpfcutoff');
  sess.hpf = max([128 EXP.TR*EXP.NumFrames/10]); % one 10th of total length
else
  sess.hpf = EXP.hpfcutoff;
end
EXP.hpfcutoff = sess.hpf;

matlabbatch{1}.spm.stats.fmri_spec.sess = sess;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
% this calls spm_fMRI_design, spm_fmri_spm_ui.

%% 2. review the design matrix and save it.
if ~exist([EXP.dir_glm,'/SPM.mat'],'file') || overwrite
  save([EXP.dir_glm,'/glm_design.mat'], 'matlabbatch');
  spm_jobman('run', matlabbatch)
end
fname_resultps=[EXP.dir_glm,'/spm_',datestr(now,'yyyymmmdd'),'.ps'];
if exist(fname_resultps,'file'), delete(fname_resultps); end
spm_print(fname_resultps);
end