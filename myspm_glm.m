function EXP=myspm_glm (EXP)
% EXP=myspm_glm (EXP)
%
% examples of the fields of EXP:
% EXP.dir_name=['/scr/vatikan3/Roberta/TPIAMO/dparsf/GLMs/meanFD_s8_',ROI_NAMES{i},'_',BHV_NAMES{j}];
%
% EXP.files_query=['/scr/vatikan3/Roberta/TPIAMO/dparsf/Results/FC_Roberta/s8zROI',num2str(i),'FCMap_*.nii'];
% or EXP.filenames=<1xN> cell
%
% EXP.fwhm = a scalr for isotropic smoothing or a vector for aisotropic smoothing
% EXP.design = either 'mreg' 't1'
% EXP.vn(c).val, EXP.vn(c).name (variable of nuissance)
% EXP.vi.val, EXP.vi.name (variable of interest)
% EXP.masking = full path for an explicit mask or method of normalization
% "Dartel_vbm8" "Dartel_vbm12"
% EXP.thresh.desc = either 'FWE','none', or 'cluster' (initial alpha=0.001, k=0)
% EXP.thresh.alpha = 0.05 (default)
% EXP.thresh.extent (# fo voxels, optional)
% EXP.titlestr <1 x #ofContrast>
% EXP.fname_struct = '$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz' (by default)
%
% (cc) 2015. sgKIM. solleo@gmail.com https://ggooo.wordpress.com/

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
  for n=1:numel(fnames)
    [a,b,c]=fileparts(fnames{n,1});
    fnames{n,1} = [a '/' prefix b c];
  end
  if ~exist(fnames{end,1}(1:end-2),'file')
    spm_jobman('initcfg')
    spm_jobman('run', matlabbatchs)
  end
  
end


if ~isfield(EXP,'design')
  design='mreg';
else
  design = EXP.design;
end

%% design specification

switch design
  case 'mreg'
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
  case 't1'
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = fnames;
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

switch design
  case 'mreg'
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [X_int +1 X_cov];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [X_int -1 X_cov];
  case 't1'
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = +1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = -1;
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
