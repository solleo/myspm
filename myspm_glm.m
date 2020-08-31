function JOB = myspm_glm(JOB)
% JOB = myspm_glm(JOB)
%
% This script helps you to set and run GLMs. A result report of
% the 1st-level GLM will be created using myspm_result.m and myspm_graph.m
%
% JOB requires for myspm_glm.m:
% -output directory
%  .dir_base     '1xN' directory for a subdirectory that has SPM.mat
% (.dir_glm)     '1xN' directory to save SPM results
% (.dir_prefix)  '1xN' prefix for GLM directory name
%
% -input files
%  .subjID       [NxM] or {Nx1}
%  .files_query  'Nx1' a query to find image filenames with "${subj}",
%                 which will be replaced by .subjID
% or
%  .filenames    {Nsubjx1} (instead of files_query)
%
% (.masking)     '1xN' filename for an explicit (inclusive) mask
%  .design       '1xN' GLM design: either multiple regression ('mreg')
%                 or one-sample t-test ('t1') or paired t-test ('pt')
%
% -model specification for JOB.design='mreg':
%  .vi.val       [Nsubjx1] a vector of interest
%  .vi.name      'string' a name of interest
%  .vn(c).val    [Nsubjx1] a vector of c-th nuissance variable
%  .vn(c).name   'string' a name of c-th nuissance variable
% or
%  .model        <term> SurfStat term structure that describes a GLM
%  .cidx         [1x1] 1-based index for the contrast of interest
%
% -for myspm_cntrst.m:
% (.cntrstMtx)
% (.titlestr)
% (.effectOfInterest)
% (.FcntrstMtx)
% (.Ftitlestr)
%
% optionally for myspm_result.m:
% (.thres.desc)    '1xN'  'FWE','none', or 'cluster'(default)
% (.thres.alpha)   [1x1]  alpha level (default=0.05)
% (.thres.extent)  [1x1]  extent threshold of clusters in voxels (default=0)
% (.thres.clusterInitAlpha) [1x1] cluster forming height threshold (default=0.001)
% (.thres.clusterInitExtent)  <1x1> cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct)   '1xN' fullpath filename for background anatomical image for orthogonal slices
%                         (defulat='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)       {1xNcont} Title text for SPM result report (default={'positive','negative'})
% (.dir_sum)        '1xN' a summary directory into where you want to copy significant results
% (.append)         [1x1] whether to append results into an existing report (default=0)
% (.print)          [1x1] whether to generate results (default=1)
% (.mygraph.x_name) '1xN' x-axis label in a scatterplot (default='x')
% (.mygraph.y_name) '1xN' y-axis label in a scatterplot (default='y')
% (.atlas)          '1xN' atlas to find anatomical names: 'fsl' (default) or 'spm12'
% (.keepResidual)   [1x1] Keeping standardized residuals
%
% Example: one-sample t-test with a covariate
%
% JOB=[];
% JOB.dir_base = '/where/I/want/to/create/dir_glm/';
% JOB.subjID   = {'subj1','subj2'};
% JOB.files_query = '/where/I/have/data/${subj}/fmridata/preproced_epi.nii';
% JOB.design  = 't1';
% JOB.vn.val  = age;
% JOB.vn.name = 'age';
% myspm_glm(JOB)
%
% Results:
%
% (cc) 2015. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if nargin == 0, help(mfilename); return; end
if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end
overwrite=JOB.overwrite;
if ~isfield(JOB,'design'),design='mreg'; else design = JOB.design; end
spm('Defaults','fmri')

matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {};
% check data
if isfield(JOB,'fnames')
  fnames = JOB.fnames;
elseif isfield(JOB,'files_query')
  JOB.subjID=fsss_subjID(JOB.subjID);
  fnames=cell(numel(JOB.subjID),1); % this MUST be a column vector <Nx1>
  idx = strfind(JOB.files_query,'${subj}');
  prefix = JOB.files_query(1:idx-1);
  suffix = JOB.files_query(idx+7:end);
  for n=1:numel(JOB.subjID)
    [~,res] = mydir([prefix,JOB.subjID{n},suffix]);
    if isempty(res)
      error(['File not found: ',prefix,JOB.subjID{n},suffix]);
    end
    fnames{n,1}=[res,',1'];
  end
elseif isfield(JOB,'filenames')
  for n=1:size(JOB.filenames,1)
    fnames{n,1} = [JOB.filenames{n,1},',1'];
  end
else
  error('You need to specify inputs in JOB.files_query or JOB.filenames');
end
Nsubj=numel(fnames);
for n=1:Nsubj
  ls(fnames{n,1}(1:end-2));
end

% check second scans for paired t-test
if strcmpi(design,'pt')
  if isfield(JOB,'files_query2')
    idx = strfind(JOB.files_query2,'${subj}');
    prefix = JOB.files_query2(1:idx-1);
    suffix = JOB.files_query2(idx+7:end);
    for n=1:numel(JOB.subjID)
      [~,res] = mydir([prefix,JOB.subjID{n},suffix]);
      if isempty(res)
        error(['File not found: ',prefix,JOB.subjID{n},suffix]);
      end
      fnames{n,2} = [res,',1'];
    end
  elseif size(JOB.filenames,2)
    for n=1:size(JOB.filenames,1)
      fnames{n,2} = [JOB.filenames{n,2},',1'];
      
    end
  else
    error('You need to specify inputs in JOB.files_queryend or JOB.filenames');
  end
  for n=1:Nsubj
    ls(fnames{n,2}(1:end-2));
  end
end

%% design specification
switch design
  case 'mreg' % multiple regression
    JOB = myspm_strcNterm(JOB);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames; % THIS HAS TO BE <Nx1>, not <1xN>!!!
    % variable of interest
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = double(JOB.vi.val);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = JOB.vi.name;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1; % centering
    % including intercept
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    if isfield(JOB,'nointercept')
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.multi_cov ...
      = struct('files', {}, 'iCFI', {}, 'iCC', {}); % for spm12
    
  case 't1' % one-sample t-test
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fnames;
    JOB.vi.name = '1';
    JOB.cntrstMtx = [1; -1];
    
  case 'pt' % paired t-test
    for n=1:size(fnames,1)
      matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(n).scans ...
        = {fnames{n,1};fnames{n,2}};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    JOB.vi.name='1';
end
if isfield(JOB,'vn') % of variable of nuissance
  for c=1:numel(JOB.vn)
    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = double(JOB.vn(c).val);
    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = JOB.vn(c).name;
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interaction none
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing overall mean
  end
else
  matlabbatch{1}.spm.stats.factorial_design.cov ...
    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
% if ~isfield(JOB,'masking')
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
% else
%   matlabbatch{1}.spm.stats.factorial_design.masking.em = {[JOB.masking,',1']};
% end
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
if ~isfield(JOB,'model_desc')
  JOB.model_desc = design;
end
if ~isfield(JOB,'dir_glm')%&&isfield(JOB,'dir_base');
  if ~isfield(JOB,'dir_prefix'), JOB.dir_prefix = ''; end
  JOB.dir_glm = fullfile(JOB.dir_base,[JOB.dir_prefix,JOB.model_desc]);
end
[~,~] = mkdir(JOB.dir_glm);
matlabbatch{1}.spm.stats.factorial_design.dir = {JOB.dir_glm};
save([JOB.dir_glm,'/glm_design.mat'], 'matlabbatch');
need2est = 1;
if overwrite
  unix(['rm -f ',JOB.dir_glm,'/SPM.mat']);
else
  if exist([JOB.dir_glm,'/SPM.mat'],'file')
    need2est = 0;
  end
end
if need2est
  spm_jobman('initcfg')
  spm_jobman('run', matlabbatch)
end
cd(JOB.dir_glm)
spm_print;

%% Coefficeint estimation
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[JOB.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
save([JOB.dir_glm,'/glm_estimation.mat'], 'matlabbatch');
if overwrite || ~exist([JOB.dir_glm,'/beta_0001.nii'],'file')
  spm_jobman('run', matlabbatch)
end


%% when required, create a NIFTI file before deleting ResI.*
pwd0 = pwd;
cd (JOB.dir_glm)
if isfield(JOB,'keepResidual') && JOB.keepResidual
  % create NIFTI
  setenv('FSLOUTPUTTYPE','NIFTI_GZ')
  myunix('fslmerge -t ResI ResI_????.nii');
  ls('ResI.nii.gz')
  setenv('FSLOUTPUTTYPE','NIFTI')
end
% delete
myunix('rm -f ResI_????.nii');
cd (pwd0);


%% 4. Design review & Orthogonality check
% set figure filename
if ~isfield(JOB,'fname_spm_fig')
  today = datestr(now,'yyyymmmdd');
  JOB.fname_spm_fig = fullfile(JOB.dir_glm,[JOB.model_desc,'_spm_',today,'.ps']);
  JOB.fname_spm_fig = strrep(JOB.fname_spm_fig,'>','-gt-');
  JOB.fname_spm_fig = strrep(JOB.fname_spm_fig,'<','-lt-');
end
% delete previous figure files today (rewriting)
if exist(JOB.fname_spm_fig,'file')
  delete(JOB.fname_spm_fig);
end

% review the design matrix and save it:
toDisplay = {'matrix','orth','covariance'};
spm_figure('clear')
for j=1:numel(toDisplay)
  matlabbatch = {};
  matlabbatch{1}.spm.stats.review.spmmat = {[JOB.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.review.display.(toDisplay{j})= 1;
  matlabbatch{1}.spm.stats.review.print = false;
  spm_jobman('run', matlabbatch);
  spm_print(JOB.fname_spm_fig);
end


%% Now create result reports
if ~isfield(JOB,'thres')
  JOB.thres.desc  = 'cluster';
  JOB.thres.alpha = 0.05;
end
if ~(isfield(JOB,'NOCNTRST') && JOB.NOCNTRST)
  myspm_cntrst (JOB);
end

end
