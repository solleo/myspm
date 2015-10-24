function EXP=myspm_glm (EXP)
% EXP=myspm_glm (EXP)
%
% This script helps you to set and run GLMs. A result report of the 1st-level GLM 
%   will be created using myspm_result.m and myspm_graph.m
%
%
% EXP requires for myspm_glm.m:
% -output directory
%  .dir_glm      'Nx1' directory to save SPM results
% or
%  .dir_base     'Nx1' directory for a subdirectory that has SPM.mat
% (.dir_prefix)
%
% -input files
%  .subjID       [NxM] or {Nx1}
%  .files_query  'Nx1' a query to find image filenames. Wildcard (*) can be used
%                except "${subj}", which will be replaced by given .subjID
% or
%  .filenames    {Nsubjx1} (instead of files_query)
%
% (.fwhm)        [1x1|1x3] 3-D smoothing kernel size in mm
% (.masking)     'Nx1' filename for an explicit (inclusive) mask
%  .design       'Nx1' type of GLM design: either multiple regression ('mreg'),
%                       or one-sample t-test ('t1') or paired t-test ('pt')
%
% -for model specification (i.e. when EXP.design='mreg'):
%  .vi.val       [Nsubjx1] a vector of interest
%  .vi.name      'string' a name of interest
%  .vn(c).val    [Nsubjx1] a vector of c-th nuissance variable
%  .vn(c).name   'string' a name of c-th nuissance variable
% or 
%  .model        <term> SurfStat term structure that describes a GLM
%  .cidx         [1x1] 1-based index for the contrast of interest
%
% optionally for myspm_result.m:
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
%
% Example:
%
% Results:
%
%
% (cc) 2015. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if nargin<1, help myspm_glm; return; end

if ~isfield(EXP,'design')
  design='mreg';
else
  design = EXP.design;
end
spm('Defaults','fmri')

matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {};
% check data
if isfield(EXP,'fnames')
  fnames = EXP.fnames;
elseif isfield(EXP,'files_query')
  EXP.subjID=fsss_subjID(EXP.subjID);
  fnames=cell(numel(EXP.subjID),1); % this MUST be <Nx1>
  idx = strfind(EXP.files_query,'${subj}');
  prefix = EXP.files_query(1:idx-1);
  suffix = EXP.files_query(idx+7:end);
  for n=1:numel(EXP.subjID)
    [~,res] = mydir([prefix,EXP.subjID{n},suffix]);
    fnames{n,1}=[res{1},',1'];
  end
elseif isfield(EXP,'filenames')
  for n=1:size(EXP.filenames,1)
    fnames{n,1} = [EXP.filenames{n,1},',1'];
  end
else
  error('You need to specify inputs in EXP.files_query or EXP.filenames');
end

Nsubj=numel(fnames);

% check second scans for paired t-test
if strcmpi(design,'pt')
  if isfield(EXP,'files_query2')
    files = dir(EXP.files_query2);
    [mypath,~,~] = fileparts(EXP.files_query2);
    for n=1:numel(files)
      fnames{n,2} = [mypath,'/',files(n).name,',1'];
    end
  elseif size(EXP.filenames,2)
    for n=1:size(EXP.filenames,1)
      fnames{n,2} = [EXP.filenames{n,2},',1'];
    end
  else
    error('You need to specify inputs in EXP.files_queryend or EXP.filenames');
  end
  for n=1:Nsubj
    ls(fnames{n,2}(1:end-2));
  end
end

%% isotropic smoothing
if isfield(EXP,'fwhm')
  EXP.fnames = fnames;
  EXP = myspm_smooth(EXP);
  fnames = EXP.fnames;
end

%% design specification

switch design
  case 'mreg' % multiple regression
    EXP = myspm_strcNterm(EXP);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames; % THIS HAS TO BE <Nx1>, not <1xN>!!!
    % variable of interest
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = double(EXP.vi.val);
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = EXP.vi.name;
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1; % centering
    % including intercept
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
    if isfield(EXP,'nointercept')
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
    end
    matlabbatch{1}.spm.stats.factorial_design.multi_cov ...
      = struct('files', {}, 'iCFI', {}, 'iCC', {}); % for spm12
  case 't1' % one-sample t-test
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fnames;
    EXP.vi.name='1';
  case 'pt' % paired t-test
    for n=1:size(fnames,1)
      matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(n).scans ...
        = {fnames{n,1};fnames{n,2}};
    end
    matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
    EXP.vi.name='1';
end
if isfield(EXP,'vn') % of variable of nuissance
  for c=1:numel(EXP.vn)
    matlabbatch{1}.spm.stats.factorial_design.cov(c).c = double(EXP.vn(c).val);
    matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = EXP.vn(c).name;
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interation
    matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing
  end
else
  matlabbatch{1}.spm.stats.factorial_design.cov ...
    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
end

matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
if ~isfield(EXP,'masking')
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
else
  if isnumeric(EXP.masking)
    %EXP.masking ='/scr/vatikan1/skim/matlab/spm12/tpm/mask_ICV.nii';
     mask_mni = [spm('dir'),'/tpm/gray',num2str(EXP.masking),'.nii'];
%     %if ~exist(mask_mni,'file')
      nii = load_untouch_nii([spm('dir'),'/tpm/TPM.nii']);
      nii.img = double(nii.img(:,:,:,1)>EXP.masking);
      nii.hdr.dime.datatype=2;
      nii.hdr.dime.dim(5)=1;
      save_untouch_nii(nii,mask_mni);
%     %end
%     [path1,~,~] = fileparts(fnames{1,1}(1:end-2));
%     mask_indi = [path1,'/ogray_',num2str(EXP.masking),'.nii'];
%     unix(['mri_convert --like ',fnames{1,1}(1:end-2),' ',mask_mni,' ',mask_indi]);
     EXP.masking = mask_mni;
  end
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {EXP.masking};
end

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

if ~isfield(EXP,'model_desc')
  EXP.model_desc=design;
end

if ~isfield(EXP,'dir_glm')&&isfield(EXP,'dir_base');
  if ~isfield(EXP,'dir_prefix'), EXP.dir_prefix=''; end
  EXP.dir_glm=fullfile(EXP.dir_base,[EXP.dir_prefix,EXP.model_desc]);
end
[~,~]=mkdir(EXP.dir_glm);
matlabbatch{1}.spm.stats.factorial_design.dir = {EXP.dir_glm};


%% Run glm
save([EXP.dir_glm,'/glm_design.mat'], 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch)
cd(EXP.dir_glm)
spm_print;

%% Coefficeint estimation
matlabbatch={};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[EXP.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

%% Contrast setting
matlabbatch{2}.spm.stats.con.spmmat = {[EXP.dir_glm,'/SPM.mat']};
matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = ['+',EXP.vi.name];
matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = ['-',EXP.vi.name];
matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.delete = 1;
if isfield(EXP,'vn')
  X_cov=zeros(1,numel(EXP.vn));
else
  X_cov=[];
end

if isfield(EXP,'vi')
  if isfield(EXP.vi,'val')
    if isempty( EXP.vi.val )
      X_int=[];
    elseif EXP.vi.val == 1
      X_int=[];
    else
      X_int=0;
    end
  end
else
  X_int=0;
end
switch design
  case 'mreg'
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec = [X_int X_cov +1];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec = [X_int X_cov -1];
  case 't1'
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec = [+1 X_cov];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec = [-1 X_cov];
  case 'pt'
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec = [1 -1 X_cov zeros(1,size(fnames,1))];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec = [-1 1 X_cov zeros(1,size(fnames,1))];
end
EXP.NumCntrst=2;

%% Run glm
save([EXP.dir_glm,'/glm_estimation.mat'], 'matlabbatch');
spm_jobman('run', matlabbatch)

%% Now create result reports
if ~isfield(EXP,'mygraph')
  EXP.mygraph.y_name='y';
  EXP.mygraph.x_name='x';
end
if ~isfield(EXP,'thresh')
  EXP.thresh.desc='cluster';
  EXP.thresh.alpha=0.05;
end
if isfield(EXP,'noresult')&&EXP.noresult
  return
else
  myspm_result(EXP)
end

end
