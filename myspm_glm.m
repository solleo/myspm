function EXP=myspm_glm (EXP)
% EXP=myspm_glm (EXP)
%
% This script helps you to run GLMs, and forward the ouput
%   to myspm_result.m and myspm_graph.m to create cluster visualization,
%   and plots using SPM
%
% Inputs:
% required fields for myspm_glm.m:
%   EXP.dir_glm     <string> directory to save SPM results
%   EXP.files_query  <string> for query to find image filenames,
%                           or just NumSubj filenames in a <1xNumSubj> cell
%   EXP.filenames    <cell> (instead of files_query)
%   EXP.design       <string> type of GLM design: either multiple regression ('mreg'),
%                           or one-sample t-test ('t1') or paired t-test ('pt')
%
% required fields for myspm_glm.m when EXP.design='mreg'
%   EXP.vi.val       <NumSubjx1:num> a vector of interest
%   EXP.vi.name      <string> a name of interest
%   EXP.vn(c).val    <NumSubjx1:num> a vector of c-th nuissance variable
%   EXP.vn(c).name   <string> a name of c-th nuissance variable
% or SurfStat Term structure that describes a GLM
%   EXP.model        <term>
%   EXP.cidx         <scalar> 1-based index for the contrast of interest
%
% optional fields for myspm_glm.m:
%   EXP.fwhm         <1x1:num> a scalr for of an isotropic smoothing
%                 or <1x3:num> a vector for anisotropic smoothing
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

if ~isfield(EXP,'design')
  design='mreg';
else
  design = EXP.design;
end

spm('Defaults','fmri')
[~,~]=mkdir(EXP.dir_glm);

matlabbatch{1}.spm.stats.factorial_design.dir = {EXP.dir_glm};
if isfield(EXP,'files_query') % this could be bad
  fnames={};
  if isfield(EXP,'subjID')
    idx = strfind(EXP.files_query,'*');
    prefix = EXP.files_query(1:idx-1);
    suffix = EXP.files_query(idx+1:end);
    for n=1:numel(EXP.subjID)
      fnames{n,1}=[prefix,EXP.subjID{n},suffix,',1'];
    end
  else
    [~,fnames] = mydir(EXP.files_query);
    if ~isempty(fnames)
      for n=1:numel(fnames)
        fnames{n,1} = [fnames{n},',1'];
      end
    else
      txt = ls(EXP.files_query);
      idx = strfind(txt,sprintf('\n'));
      idx = [0 idx];
      for n=1:numel(idx)-1
        fnames{n,1} = [txt(idx(n)+1:idx(n+1)-1),',1'];
      end
    end
  end
elseif isfield(EXP,'filenames')
  for n=1:numel(EXP.filenames)
    fnames{n,1} = [EXP.filenames{n},',1'];
  end
else
  error('You need to specify inputs in EXP.files_query or EXP.filenames');
end

Nsubj=numel(fnames);
for n=1:Nsubj
  ls(fnames{n,1}(1:end-2));
end

if strcmpi(design,'pt')
  if isfield(EXP,'files_query2')
    files = dir(EXP.files_query2);
    [mypath,~,~] = fileparts(EXP.files_query2);
    for n=1:numel(files)
      fnames{n,2} = [mypath,'/',files(n).name,',1'];
    end
  elseif isfield(EXP,'filenames2')
    for n=1:numel(EXP.filenames2)
      fnames{n,2} = [EXP.filenames{n},',1'];
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
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatchs)
  end
end

%% design specification
if ~isfield(EXP,'model')
  switch design
    case 'mreg' % multiple regression
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames;
      % variable of interest
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = EXP.vi.val;
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = EXP.vi.name;
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1; % centering
      % including intercept
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
      if isfield(EXP,'nointercept')
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
      end
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
      matlabbatch{1}.spm.stats.factorial_design.cov(c).c = EXP.vn(c).val;
      matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = EXP.vn(c).name;
      matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interation
      matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing
    end
  else
    matlabbatch{1}.spm.stats.factorial_design.cov ...
      = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
  end
else % for model term
  error('Now CODE model term!');
  switch design
    case 'mreg' % multiple regression
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames;
      % variable of interest
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = EXP.vi.val;
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = EXP.vi.name;
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1; % centering
      % including intercept
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
      if isfield(EXP,'nointercept')
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 0;
      end
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
      matlabbatch{1}.spm.stats.factorial_design.cov(c).c = EXP.vn(c).val;
      matlabbatch{1}.spm.stats.factorial_design.cov(c).cname = EXP.vn(c).name;
      matlabbatch{1}.spm.stats.factorial_design.cov(c).iCFI = 1; % interation
      matlabbatch{1}.spm.stats.factorial_design.cov(c).iCC = 1;  % centing
      matlabbatch{1}.spm.stats.factorial_design.cov ...
        = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    end
  end
end
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
if ~isfield(EXP,'masking')
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
else
  matlabbatch{1}.spm.stats.factorial_design.masking.em = {EXP.masking};
end

matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

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
myspm_result(EXP)

end
