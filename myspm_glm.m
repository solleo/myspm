function EXP=myspm_glm (EXP)
% EXP=myspm_glm (EXP)
%
% This script helps you to run GLMs, and forward the ouput
%   to myspm_result.m and myspm_graph.m to create cluster visualization,
%   and plots using SPM
%
% Inputs:
% <for myspm_glm.m>
% EXP.dir_name <string> directory to save SPM results
% EXP.files_query <string> for query to find image filenames,
%   or just N filenames in <1xN> cell
% EXP.fwhm <1x1> a scalr for isotropic smoothing
%   or <1x3> a vector for anisotropic smoothing
% EXP.design <string> design type either multiple regression ('mreg')
%   or one-sample t-test ('t1') or paired t-test ('pt')
% EXP.vn(c).val <Nx1> a vector of c-th nuissance variable
% EXP.vn(c).name <string> a name of c-th nuissance variable
% EXP.vi.val <Nx1> a vector of interest
% EXP.vi.name <string> a name of interest
% EXP.masking <string> filename for an explicit mask
%   or method of normalization as "Dartel_vbm8" or "Dartel_vbm12"
%
% <for myspm_result.m>
% EXP.thresh.desc = either 'FWE','none', or 'cluster' (initial alpha=0.001, k=0)
% EXP.thresh.alpha = 0.05 (default)
% EXP.thresh.extent (# fo voxels, optional)
% EXP.thresh.clusterInitAlpha = 0.001 (default) for a cluster-threshold
% EXP.thresh.clusterInitExtent = 10 (voxels; default)
% EXP.titlestr <1 x #ofContrast>
% EXP.fname_struct = '$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz' (default)
% EXP.titlestr = {'positive','negative'} (default)
% EXP.dir_sum = '' (a summary directory where you want to copy 'significant' results into)
% EXP.append = 0 (default)
% EXP.print  = 1 (default)
% EXP.mygraph.y_name = 'x'
% EXP.mygraph.x_name = 'y' for the scatterplots and summary tables
% EXP.atlas = 'fsl' (default) or 'spm12'
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

mkdir(EXP.dir_name)
matlabbatch{1}.spm.stats.factorial_design.dir = {EXP.dir_name};

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
end

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

switch design
  case 'mreg' % multiple regression
    matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = fnames;
    if isfield(EXP,'vi') % variable of interest
      for c=1:numel(EXP.vi)
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(c).c = EXP.vi(c).val;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(c).cname = EXP.vi(c).name;
        matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(c).iCC = 1; % centering
      end
    else
      matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov ...
        = struct('c', {}, 'cname', {}, 'iCC', {});
    end
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

if isfield(EXP,'vn'); %isempty(EXP.vi.val) % of variable of nuissance
  matlabbatch{1}.spm.stats.factorial_design.cov.c = EXP.vn.val;
  matlabbatch{1}.spm.stats.factorial_design.cov.cname = EXP.vn.name;
  matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1; % interation
  matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;  % centing
else
  matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
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

%% Coefficeint estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).sname = 'Factorial design specification: SPM.mat File';
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Design specification
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep;
matlabbatch{3}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(1).value = 'mat';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.spm.stats.con.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.spm.stats.con.spmmat(1).sname = 'Model estimation: SPM.mat File';
matlabbatch{3}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ['+',EXP.vi.name];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['-',EXP.vi.name];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

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
    end
  end
else
  X_int=0;
end

% contrast setting
switch design
  case 'mreg'
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [X_int +1 X_cov];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [X_int -1 X_cov];
  case 't1'
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [+1 X_cov];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 X_cov];
  case 'pt'
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 -1 X_cov zeros(1,size(fnames,1))];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [-1 1 X_cov zeros(1,size(fnames,1))];
end


%% Run glm
save([EXP.dir_name,'/glm.mat'], 'matlabbatch');
spm_jobman('initcfg')
spm_jobman('run', matlabbatch)

%% Now create result reports
if ~isfield(EXP,'mygraph')
  EXP.mygraph.y_name='y';
  EXP.mygraph.x_name='x';
end
myspm_result(EXP)

end
