function EXP=myspm_fmriglm_par (EXP)
% EXP=myspm_fmriglm_par (EXP)
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


%% 0. parsing & check inputs
[~,hostname]=whereami;
if ~strcmp(hostname,'hayd')
  error('Run MATLAB parallel at Hayd!!');
end
delete(gcp);
spm('Defaults','fmri');
spm_jobman('initcfg');
matlabpool('open',12);

if isfield(EXP,'model')
  model_desc = fsss_model_desc(EXP.model);
else
  model_desc = EXP.model_desc;
end

mrtID = fsss_subjID(EXP.subjID);
NumSubj = numel(mrtID);
dir0  = EXP.dir0;

for i=1:NumSubj
  % find the filenames with comma-frameIndex
  FNAMES_comma={};
  EXP.dir_base  = [dir0,mrtID{i}];
  EXP.dir_glm = [EXP.dir_base,'/glm_',model_desc];
  [~,~]=mkdir(EXP.dir_glm);
  
  EXP.NumSess   = numel(EXP.name_data);
  for j=1:EXP.NumSess % for each session
    EXP.filenames{j} = [dir0,mrtID{i},'/',EXP.name_data{j}];
    ls(EXP.filenames{j});
    hdr   = load_untouch_header_only(EXP.filenames{j});
    EXP.NumFrames(j) = hdr.dime.dim(5);
    for t=1:EXP.NumFrames(j)  % for each frame
      FNAMES_comma{j}{t} = [EXP.filenames{j},',',num2str(t)];
    end
  end
  
%   % find conditions
%   if ~isfield(EXP,'COND')
%     EXP.NumCond = numel(EXP.name_cond);
%   else
%     EXP.NumCond = numel(EXP.COND);
%   end
  EXP.NumCond = numel(EXP.name_cond);
  
  %% 0. isotropic smoothing
  if isfield(EXP,'fwhm')
    for j=1:EXP.NumSess
      matlabbatchs={};
      if numel(EXP.fwhm) == 1
        fwhm = [EXP.fwhm EXP.fwhm EXP.fwhm];
      else
        fwhm = EXP.fwhm;
      end
      matlabbatchs{1}.spm.spatial.smooth.data = FNAMES_comma{j};
      matlabbatchs{1}.spm.spatial.smooth.fwhm = fwhm;
      matlabbatchs{1}.spm.spatial.smooth.dtype = 0;
      matlabbatchs{1}.spm.spatial.smooth.im = 0;
      prefix=['s' num2str(round(mean(fwhm))) '_'];
      matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;
      for t=1:EXP.NumFrames(j)
        [a,b,c]=fileparts(FNAMES_comma{j}{t});
        FNAMES_comma{j}{t} = [a '/' prefix b c];
      end
      if ~exist(FNAMES_comma{j}{end}(1:end-2),'file')
        spm_jobman('run', matlabbatchs)
      end
    end
  end
end

%% 1. model specification
anything2run=0;
mrtID2run={};
for i=1:NumSubj
  if ~exist([dir0,'/',mrtID{i},'/glm_',model_desc,'/SPM.mat'],'file')
    anything2run=1;
    mrtID2run = [mrtID2run mrtID{i}];
  end
end

%% condition setting in EXP
for k=1:EXP.NumCond
  EXP.COND(k).name = EXP.name_cond{k};
  EXP.COND(k).onset = dlmread([EXP.dir_base,'/',EXP.name_cond{k},'.ons']);
  if isfield(EXP,'name_param')
    if exist([EXP.dir_base,'/',EXP.name_param{k},'.par'],'file')
      EXP.isParamFile = 1;
      EXP.COND(k).pmod.name  = EXP.name_param{k};
      EXP.COND(k).pmod.param = dlmread([EXP.dir_base,'/',EXP.name_cond{k},'.par']);
    else
      EXP.isParamFile=0;
    end
  else
    EXP.isParamFile=0;
  end
end

%%
NumSubj2run = numel(mrtID2run);
EXPS = cell(1,NumSubj2run);
for i=1:NumSubj2run, EXPS{i}=EXP; end
MATLABBATCH = cell(1,NumSubj2run);

% NumCores=min(12,numel(mrtID2run));
% matlabpool('open',NumCores);

parfor i=1:NumSubj2run
  MATLABBATCH{i}={};
  MATLABBATCH{i}{1}.spm.stats.fmri_spec.dir = {EXPS{i}.dir_glm};
  MATLABBATCH{i}{1}.spm.stats.fmri_spec.timing.units = 'secs';
  MATLABBATCH{i}{1}.spm.stats.fmri_spec.timing.RT = EXPS{i}.TR;
  MATLABBATCH{i}{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
  MATLABBATCH{i}{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % middle point
  
  for j=1:EXPS{i}.NumSess
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).scans = FNAMES_comma{j};
    for k=1:EXPS{i}.NumCond
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).name     = EXPS{i}.COND(k).name;
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).onset    = EXPS{i}.COND(k).onset;
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).duration = EXPS{i}.StimDurSec;
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).tmod = 0; % time modulation
      
      % now for parametric design
      if EXPS{i}.isParamFile
        MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod.name  = EXPS{i}.COND(k).pmod.name;
        MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod.param = EXPS{i}.COND(k).pmod.param;
        if ~isfield(EXPS{i}.COND(k).pmod,'poly')
          MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod.poly  = 1;
        else
          MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod.poly  = EXPS{i}.COND(k).pmod.poly;
        end
      else
        MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
      end
    end
    
    % rigid-motion parameters
    if isfield(EXPS{i},'name_rp')
      rp=dlmread([EXPS{i}.dir_base,'/',EXPS{i}.name_rp]);
      rpname={'dx','dy','dz','rx','ry','rz'};
      for k=1:6
        MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).regress(k).name = rpname{k};
        MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).regress(k).val = rp(:,k);
      end
    end
    
    if ~isfield(EXPS{i},'hpfcutoff');
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).hpf = EXPS{i}.TR*EXPS{i}.NumFrames(j)/10; % one 10th of total length
    else
      MATLABBATCH{i}{1}.spm.stats.fmri_spec.sess(j).hpf = EXPS{i}.hpfcutoff;
    end
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.volt = 1;
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.global = 'None';
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.mask = {''};
    MATLABBATCH{i}{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  end
  spm_jobman('initcfg');
  spm_jobman('run', MATLABBATCH{i})
end

%% 2. review the design matrix and save it.
%   fname_resultps=[EXPS{i}.dir_glm,'/spm_',datestr(now,'yyyymmmdd'),'.ps'];
%   if exist(fname_resultps,'file'), delete(fname_resultps); end
%   spm_print(fname_resultps);



%% 3. model estimation (takes SO LONG!)
anything2run=0;
mrtID2run={};
for i=1:NumSubj
  if ~exist([dir0,'/',mrtID{i},'/glm_',model_desc,'/mask.img'],'file')
    anything2run=1;
    mrtID2run = [mrtID2run mrtID{i}];
  end
end

if anything2run
%   NumCores=min(12,numel(mrtID2run));
%   matlabpool('open',NumCores);
  MATLABBATCH=cell(1,NumCores);
  NumSubj2run = numel(mrtID2run);
  EXPS=cell(1,NumSubj2run);
  for i=1:NumSubj2run, EXPS{i}=EXP; end
  parfor i=1:NumSubj2run
    EXPS{i}.dir_base = [dir0,'/',mrtID2run{i}];
    EXPS{i}.dir_glm  = [EXPS{i}.dir_base,'/glm_',model_desc];
    MATLABBATCH{i}{1}=[];
    MATLABBATCH{i}{1}.spm.stats.fmri_est.spmmat = {[EXPS{i}.dir_glm,'/SPM.mat']};
    MATLABBATCH{i}{1}.spm.stats.fmri_est.method.Classical = 1; % ReML
    %% Estimation if not done before... let's check mask.img as SPM does
    spm_jobman('initcfg');
    spm_jobman('run', MATLABBATCH{i});
  end
%   matlabpool('close');
end

%% 4. If all done, set different contrast and printout all results
EXP.NOREPORT=1;
parfor i=1:NumSubj
  %EXP.dir_glm = [dir0,'/',mrtID{i},'/glm_',model_desc];
  myspm_cntrst (EXP{i});
end
matlabpool('close');
end
