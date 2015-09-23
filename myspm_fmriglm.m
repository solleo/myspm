function EXP=myspm_fmriglm (EXP)
% EXP=myspm_fmriglm (EXP)
%
% This script helps you to run GLMs for the 1st-level fMRI analysis
%   using fMRI model specification & specification.
%   Followingly 1st-level GLM results are created using
%   myspm_result.m and myspm_graph.m to create cluster visualization,
%   and plots using SPM
%
% Inputs:
% required fields for myspm_glm.m:
%   EXP.dir_glm     <string> directory to save SPM results
%   EXP.files_query  <string> for query to find image filenames,
%     or
%   EXP.filenames    <cell:1xN> for N images
%
%   EXP.model_desc   <str>
%
% optional fields for myspm_glm.m:
%   EXP.fwhm         <num:1x1> a scalr for of an isotropic smoothing
%   EXP.masking      <string> filename for an explicit (inclusive) mask
%
% optional fields for myspm_result.m:
%   EXP.thresh.desc    <string>  'FWE','none', or 'cluster'(default)
%   EXP.thresh.alpha   <1x1:num> alpha level (default= 0.05)
%   EXP.thresh.extent  <1x1:num> a number of voxels for extent threshold
%   EXP.thresh.clusterInitAlpha   <1x1 num> 0.001 (default) for a cluster-threshold
%   EXP.thresh.clusterInitExtent  <1x1 num> 10 (voxels; default)
%   EXP.fname_struct   <string>  '$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz' (default)
%   EXP.titlestr       <1xNumCont:cell> {'positive','negative'} (default)
%   EXP.dir_sum        <string>  '' (a summary directory where you want to copy 'significant' results into)
%   EXP.append         <1x1 num> = 0 (default)
%   EXP.print          <1x1 num> = 1 (default)
%   EXP.mygraph.x_name <string>  'x'
%   EXP.mygraph.y_name <string>  'y' for the scatterplots and summary tables
%   EXP.atlas          <string>  'fsl' (default) or 'spm12'
%   EXP.fname_spm_fig  <string>
%
% Example:
%
% Results:
%
%
% (cc) 2015. sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com/

%% -1. PPI?
if isfield(EXP,'fname_ppi');
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
    -1  0  0 ];
end

%% 0. parsing & check inputs
spm('Defaults','fmri');
spm_jobman('initcfg');

if isfield(EXP,'model')
  model_desc = fsss_model_desc(EXP.model);
else
  model_desc = EXP.model_desc;
end
if ~isfield(EXP,'dir_base')
  EXP.dir_glm = [pwd,'/glm_',model_desc];
else
  EXP.dir_glm = [EXP.dir_base,'/glm_',model_desc];
end
[~,~]=mkdir(EXP.dir_glm);

if isfield(EXP,'files_query')
  fnames={};
  files = dir(EXP.files_query);
  [mypath,~,~] = fileparts(EXP.files_query);
  for n=1:numel(files)
    fnames{n,1} = [mypath,'/',files(n).name,',1'];
  end
elseif isfield(EXP,'filenames')
  for n=1:numel(EXP.filenames)
    fnames{n,1} = [EXP.filenames{n},',1'];
  end
else
  error('You need to specify inputs in EXP.files_queryend or EXP.filenames');
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

%% 0.5. isotropic smoothing

if isfield(EXP,'fwhm')
  matlabbatchs={};
  if numel(EXP.fwhm) == 1
    fwhm = [EXP.fwhm EXP.fwhm EXP.fwhm];
  else
    fwhm = EXP.fwhm;
  end
  matlabbatchs{1}.spm.spatial.smooth.data = fnames;
  matlabbatchs{1}.spm.spatial.smooth.fwhm = fwhm;
  matlabbatchs{1}.spm.spatial.smooth.dtype = 0;
  matlabbatchs{1}.spm.spatial.smooth.im = 0;
  prefix=['s' num2str(round(mean(fwhm))) '_'];
  matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;
  for j=1:size(fnames,2)
    for n=1:size(fnames,1)
      [a,b,c]=fileparts(fnames{n,j});
      fnames{n,j} = [a '/' prefix b c];
    end
  end
  if ~exist(fnames{end,1}(1:end-2),'file')
    spm_jobman('run', matlabbatchs)
  end
end

%% 1. model specification
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.dir = {EXP.dir_glm};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % middle point

sess=[];
for j=1:EXP.NumSess
  for t=1:EXP.NumFrames
    sess(j).scans{t} = [EXP.filenames{j},',',num2str(t)];
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
  
  % rigid-motion parameters
  if isfield(EXP,'fname_rp')
    rp=dlmread(EXP.fname_rp);
    rpname={'dx','dy','dz','rx','ry','rz'};
    for k=[1:6]+NumReg
      sess(j).regress(k).name = rpname{k-NumReg};
      sess(j).regress(k).val = rp(:,k-NumReg);
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
  if isfield(EXP,'GMmask')
    if isnumeric(EXP.GMmask)
      fname_func1 = fnames{1,1}(1:end-2);
      fname_mask = [EXP.dir_glm,'/GMmask_',num2str(EXP.GMmask),'.nii'];
      fname_moving = [spm('Dir'),'/tpm/grey.nii'];
      unix(['mri_convert --like ',fname_func1,' ',fname_moving,' ',fname_mask]);
      nii1 = load_untouch_nii(fname_func1,1);
      nii2 = load_untouch_nii(fname_mask);
      figure; imageorth(threshold_prctile(nii1.img,90))
      figure; imageorth(threshold_prctile(nii2.img,80))
      matlabbatch{1}.spm.stats.fmri_spec.mask = {fname_mask};
    elseif ischar(EXP.GMmask)
      matlabbatch{1}.spm.stats.fmri_spec.mask = {EXP.GMmask};
    else
      warning('Cannot determine an explicit mask. Implicit mask used instead');
    end
  else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  end
  matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  % this calls spm_fMRI_design, spm_fmri_spm_ui.
  
  %% 2. review the design matrix and save it.
  if ~exist([EXP.dir_glm,'/SPM.mat'],'file')
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
  if ~exist([EXP.dir_glm,'/mask.img'],'file')
    save ([EXP.dir_base,'/glm_estimation.mat'],'matlabbatch');
    spm_jobman('run', matlabbatch)
  end
  
  %% and set different contrast and printout all results
  EXP.thresh.desc  = 'FWE';
  EXP.thresh.alpha = 0.05;
  EXP = myspm_cntrst (EXP);
  
end % for session

end
