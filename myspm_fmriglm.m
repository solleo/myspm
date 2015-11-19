function EXP=myspm_fmriglm (EXP)
% EXP=myspm_fmriglm (EXP)
%
% This script helps you to run GLMs for the 1st-level fMRI analysis using fMRI model
%   specification & specification. A result report of the 1st-level GLM is created
%   using myspm_result.m and myspm_graph.m. The important difference
%   from the simple GLM is that it deals with temporal aspects of fMRI data including
%   slow hedonamics, pre-whitening for auto-correlation, as well as high-pass
%   filtering for detrending.
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
% - for fMRI design matrix
%  .TR
%  .COND(k).name
%  .COND(k).onset
% (.COND(k).param)
%  .reg
%  .fname_rp
% (.hpfcutoff)   : EXP.TR*EXP.NumFrames(j)/10
% (.masking)     'Nx1' filename for an explicit (inclusive) mask
% (.fwhm)        [1x1|1x3] 3-D smoothing kernel size in mm
%
% - for PPI,
%  .fname_ppi
%
% - optionally for myspm_result.m:
% (.thresh.desc)    'Nx1'  'FWE','none', or 'cluster'(default)
% (.thresh.alpha)   [1x1]  alpha level (default=0.05)
% (.thresh.extent)  [1x1]  extent threshold of clusters in voxels (default=0)
% (.thresh.clusterInitAlpha)   <1x1> cluster forming height threshold (default=0.001)
% (.thresh.clusterInitExtent)  <1x1> cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct)   'Nx1' fullpath filename for background anatomical image for orthogonal slices
%                         (defulat='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)       {1xNcont} Title text for SPM result report (default={'positive','negative'})
% (.dir_sum)        'Nx1' a summary directory into where you want to copy significant results
% (.append)         [1x1] whether to append results into an existing report (default=0)
% (.print)          [1x1] whether to generate results (default=1)
% (.mygraph.x_name) 'Nx1' x-axis label in a scatterplot (default='x')
% (.mygraph.y_name) 'Nx1' y-axis label in a scatterplot (default='y')
% (.atlas)          'Nx1' atlas to find anatomical names: 'fsl' (default) or 'spm12'
% (.fname_rp)
% (.fname_cc)
%
% Example:
%
% Results:
%
%
% (cc) 2015. sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com/

%% -1. PPI?
if isfield(EXP,'fname_ppi')
  load(EXP.fname_ppi, 'PPI');
  EXP.reg=[];
  EXP.reg(1).name = 'ppi';
  EXP.reg(1).val  = PPI.ppi;
  EXP.reg(2).name = 'psy';
  EXP.reg(2).val  = PPI.P;
  EXP.reg(3).name = 'phy';
  EXP.reg(3).val  = PPI.Y;
  EXP.cntrstMtx = [
    +1  0  0
    -1  0  0];
  if isfield(EXP,'phy1')
    load(EXP.fname_ppi, 'PPI');
    EXP.reg=[];
    EXP.reg(1).name = 'phy';
    EXP.reg(1).val  = zscore(PPI.Y);
    EXP.cntrstMtx = [1; -1];
  end
end

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

if isfield(EXP,'files_query')
  fnames={}; % this MUST be a column vector <Sx1>
  files = dir(EXP.files_query);  % this will sort the filename!
  [mypath,~,~] = fileparts(EXP.files_query);
  for n=1:numel(files)
    fnames{n,1} = [mypath,'/',files(n).name,',1'];
  end
elseif isfield(EXP,'filenames')
  for n=1:numel(EXP.filenames)
    fnames{n,1} = [EXP.filenames{n},',1'];
  end
else
  error('You need to specify inputs in EXP.files_query or EXP.filenames');
end
ls(fnames{n,1}(1:end-2));
EXP.NumSess = numel(fnames);
for j=1:EXP.NumSess
  hdr   = load_untouch_header_only(fnames{j}(1:end-2));
  NF(j) = hdr.dime.dim(5);
end
EXP.NumFrames=NF;
if ~isfield(EXP,'COND')
  if isfield(EXP,'fname_cond')
    EXP.NumCond = numel(EXP.fname_cond);
    readOnsetFile=1;
  else
    EXP.NumCond = 0;
  end
else
  EXP.NumCond = numel(EXP.COND);
  readOnsetFile=0;
end
EXP.NumCntrst = size(EXP.cntrstMtx,1);

%% isotropic smoothing
if isfield(EXP,'fwhm')
  EXP.fnames = fnames;
  EXP = myspm_smooth(EXP);
  fnames = EXP.fnames;
end

%% 1. model specification
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''}; matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; % for spm12
matlabbatch{1}.spm.stats.fmri_spec.dir = {EXP.dir_glm};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % middle point
sess=[];
for j=1:EXP.NumSess
  for t=1:EXP.NumFrames
    sess(j).scans{t,1} = [EXP.filenames{j},',',num2str(t)]; % this must be <Tx1>
  end
  if EXP.NumCond
    for k=1:EXP.NumCond
      if readOnsetFile
        EXP.COND(k).name = EXP.fname_cond{k};
        [dir0,~,~] = fileparts_gz(fnames{j,1});
        EXP.COND(k).onset = dlmread([dir0,'/',EXP.fname_cond{k},'.ons']);
      end
      sess(j).cond(k).name     = EXP.COND(k).name;
      sess(j).cond(k).onset    = EXP.COND(k).onset;
      if numel(EXP.StimDurSec)>1
        sess(j).cond(k).duration = EXP.StimDurSec(k);
      else
        sess(j).cond(k).duration = EXP.StimDurSec;
      end
      sess(j).cond(k).tmod = 0; % time modulation
      
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
            sess(j).cond(k).pmod(i).name  = EXP.COND(k).pmod(i).name;
            sess(j).cond(k).pmod(i).param = EXP.COND(k).pmod(i).param;
            if ~isfield(EXP.COND(k).pmod,'poly')
              sess(j).cond(k).pmod(i).poly  = 1;
            else
              sess(j).cond(k).pmod(i).poly  = EXP.COND(k).pmod(i).poly;
            end
          end
        else
          sess(j).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
        end
      end
    end % loop for conditions
  else
    sess(j).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
  end
  
  if isfield(EXP,'reg')
    NumReg = numel(EXP.reg);
    sess(j).regress = EXP.reg;
  else
    NumReg=0;
  end
  
  % rigid-motion parameters & compcor
  if isfield(EXP,'fname_rp')
    rp=dlmread(EXP.fname_rp);
    rpname={'dx','dy','dz','rx','ry','rz'};
    for k=[1:6]+NumReg
      sess(j).regress(k).name = rpname{k-NumReg};
      sess(j).regress(k).val = rp(:,k-NumReg);
    end
    % and l2norm of trans. and rot.
    sess(j).regress(NumReg+6+1).val = [0;l2norm(diff(rp(:,1:3)))];
    sess(j).regress(NumReg+6+1).name = '|dd/dt|';
    sess(j).regress(NumReg+6+2).val = [0;l2norm(diff(rp(:,4:6)))];
    sess(j).regress(NumReg+6+2).name = '|dr/dt|';
  end
  if isfield(EXP,'fname_cc');
    pci=1;
    pc=dlmread(EXP.fname_cc);
    numpc=size(pc,2);
    for k=[1:numpc]+NumReg+6+2
      sess(j).regress(k).name = ['CC#',num2str(pci)];
      sess(j).regress(k).val  = pc(:,k-NumReg-6-2);
      pci=pci+1;
    end
  end
  
  if ~isfield(EXP,'hpfcutoff');
    sess(j).hpf = EXP.TR*EXP.NumFrames(j)/10; % one 10th of total length
  else
    sess(j).hpf = EXP.hpfcutoff;
  end
  
  matlabbatch{1}.spm.stats.fmri_spec.sess = sess;
  matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
  matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
  matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
  matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
  if isfield(EXP,'masking')
    if isnumeric(EXP.masking)
      fname_func1 = fnames{1,1}(1:end-2);
      fname_mask = [EXP.dir_glm,'/masking_',num2str(EXP.masking),'.nii'];
      spmversion=spm('version');
      if strcmp(spmversion(4:5),'8 ')
        fname_moving = [spm('Dir'),'/tpm/grey.nii'];
      else
        fname_tpm = [spm('Dir'),'/tpm/TPM.nii'];
        fname_moving = '/tmp/grey.nii';
        unix(['fslroi ',fname_tpm,' ',fname_moving,' 0 1']);
      end
      unix(['mri_convert --like ',fname_func1,' ',fname_moving,' ',fname_mask]);
      nii1 = load_untouch_nii(fname_func1,1);
      nii2 = load_untouch_nii(fname_mask);
      figure; imageorth(threshold_prctile(nii1.img,90))
      figure; imageorth(threshold_prctile(nii2.img,80))
      matlabbatch{1}.spm.stats.fmri_spec.mask = {fname_mask};
    elseif ischar(EXP.masking)
      matlabbatch{1}.spm.stats.fmri_spec.mask = {EXP.masking};
    else
      warning('Cannot determine an explicit mask. Using implicit mask instead..');
    end
  else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  end
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
  
  %% 3. model estimation (takes SO LONG!)
  matlabbatch={};
  fmri_est.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  fmri_est.method.Classical = 1; % ReML
  matlabbatch{1}.spm.stats.fmri_est = fmri_est;
  % Estimation if not done before... let's check mask.img as SPM does
  if ~exist([EXP.dir_glm,'/mask.img'],'file') || overwrite
    save ([EXP.dir_base,'/glm_estimation.mat'],'matlabbatch');
    spm_jobman('run', matlabbatch)
  end
  
  %% and set different contrast and printout all results
  if ~isfield(EXP,'thresh')
  EXP.thresh.desc  = 'FWE';
  EXP.thresh.alpha = 0.05;
  end
  myspm_cntrst (EXP);
  
end % for session

end
