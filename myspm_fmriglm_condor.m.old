function EXP=myspm_fmriglm_condor (EXP)
% EXP=myspm_fmriglm_condor (EXP)
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
spm('Defaults','fmri');
spm_jobman('initcfg');

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
  
  % find conditions
  if ~isfield(EXP,'COND')
    EXP.NumCond = numel(EXP.name_cond);
  else
    EXP.NumCond = numel(EXP.COND);
  end
  
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
  
  %% 1. model specification
  matlabbatch={};
  matlabbatch{1}.spm.stats.fmri_spec.dir = {EXP.dir_glm};
  matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
  matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8; % middle point
  
  sess=[];
  if ~isfield(EXP,'COND')
    readOnsetFile=1;
  end
  for j=1:EXP.NumSess
    sess(j).scans = FNAMES_comma{j};
    for k=1:EXP.NumCond
      if exist('readOnsetFile','var')
        EXP.COND(k).name = EXP.name_cond{k};
        EXP.COND(k).onset = dlmread([EXP.dir_base,'/',EXP.name_cond{k},'.ons']);
      end
      sess(j).cond(k).name     = EXP.COND(k).name;
      sess(j).cond(k).onset    = EXP.COND(k).onset;
      sess(j).cond(k).duration = EXP.StimDurSec;
      sess(j).cond(k).tmod = 0; % time modulation

      % now for parametric design
      isParamFile = exist([EXP.dir_base,'/',EXP.name_cond{k},'.par'],'file');
      isParamCOND = isfield(EXP.COND(k),'param');
      if isParamFile || isParamCOND
        if isParamFile
          if isfield(EXP,'name_param')
            EXP.COND(k).pmod.name  = EXP.name_param{k};
          else
            EXP.COND(k).pmod.name  = 'Modulation parameter';
          end
          EXP.COND(k).pmod.param = dlmread([EXP.dir_base,'/',EXP.name_cond{k},'.par']);
        end
        sess(j).cond(k).pmod.name  = EXP.COND(k).pmod.name;
        sess(j).cond(k).pmod.param = EXP.COND(k).pmod.param;
        if ~isfield(EXP.COND(k).pmod,'poly')
          sess(j).cond(k).pmod.poly  = 1;
        else
          sess(j).cond(k).pmod.poly  = EXP.COND(k).pmod.poly;
        end
      else
        sess(j).cond(k).pmod = struct('name', {}, 'param', {}, 'poly', {});
      end
    end
    
    % rigid-motion parameters
    if isfield(EXP,'name_rp')
      rp=dlmread([EXP.dir_base,'/',EXP.name_rp]);
      rpname={'dx','dy','dz','rx','ry','rz'};
      for k=1:6
        sess(j).regress(k).name = rpname{k};
        sess(j).regress(k).val = rp(:,k);
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
%     if isfield(EXP,'GMmask')
%       if isnumeric(EXP.GMmask)
%         fname_func1 = FNAMES_comma{1,1}(1:end-2);
%         fname_mask = [EXP.dir_glm,'/GMmask_',num2str(EXP.GMmask),'.nii'];
%         fname_moving = [spm('Dir'),'/tpm/grey.nii'];
%         unix(['mri_convert --like ',fname_func1,' ',fname_moving,' ',fname_mask]);
%         nii1 = load_untouch_nii(fname_func1,1);
%         nii2 = load_untouch_nii(fname_mask);
%         figure; imageorth(threshold_prctile(nii1.img,90))
%         figure; imageorth(threshold_prctile(nii2.img,80))
%         matlabbatch{1}.spm.stats.fmri_spec.mask = {fname_mask};
%       elseif ischar(EXP.GMmask)
%         matlabbatch{1}.spm.stats.fmri_spec.mask = {EXP.GMmask};
%       else
%         warning('Cannot determine an explicit mask. Implicit mask used instead');
%       end
%     else
%       matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
%     end
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  end
  % this calls spm_fMRI_design, spm_fmri_spm_ui.
  
  %% 2. review the design matrix and save it.
  if ~exist([EXP.dir_glm,'/SPM.mat'],'file')
    save([EXP.dir_glm,'/glm_design.mat'], 'matlabbatch');
    spm_jobman('run', matlabbatch)
  end
  fname_resultps=[EXP.dir_glm,'/spm_',datestr(now,'yyyymmmdd'),'.ps'];
  if exist(fname_resultps,'file'), delete(fname_resultps); end
  spm_print(fname_resultps);
end

%% 3. model estimation (takes SO LONG!)
anything2run=0;
[~,~]=mkdir([dir0,'/log_job/']);
fname_cmd = [dir0,'/log_job/SPM_est.cmd'];
fid=fopen(fname_cmd,'w');
for i=1:NumSubj
  EXP.dir_base = [dir0,'/',mrtID{i}];
  EXP.dir_glm  = [EXP.dir_base,'/glm_',model_desc];
  matlabbatch={};
  fmri_est.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  fmri_est.method.Classical = 1; % ReML
  matlabbatch{1}.spm.stats.fmri_est = fmri_est;
  %% Estimation if not done before... let's check mask.img as SPM does
  if ~exist([EXP.dir_glm,'/mask.img'],'file')
    save ([EXP.dir_base,'/glm_estimation.mat'],'matlabbatch');
    anything2run=1;
    
    fname_m = [EXP.dir_base,'/glm_estimation_run.m'];
    fid_m = fopen(fname_m,'w');
    fprintf(fid_m,'spm(''Defaults'',''fmri''); spm_jobman(''initcfg'');\n');
    fprintf(fid_m,'load %s;\n', [EXP.dir_base,'/glm_estimation.mat']);
    fprintf(fid_m,'spm_jobman(''run'', matlabbatchs)\n');
    fclose(fid_m);

    fprintf(fid,'MATLAB --version 8.2 matlab -nosplash -r "run(''%s'');exit"\n', fname_m);
  end
end
fclose(fid);

if anything2run
  % 1. submit condor jobs
  [~,qid] = unix(['fsl_sub -t ',fname_cmd]);
  % 2. wait for it
  unix(['waitForCONDORJobs.sh 100 ',qid(end-4:end-1)]);
  
end
%% 4. If all done, set different contrast and printout all results
for i=1:NumSubj
  EXP.dir_glm = [];
  EXP = myspm_cntrst (EXP);
end

end
