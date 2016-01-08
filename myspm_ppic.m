function EXP = myspm_ppic (EXP)
% EXP = myspm_ppic (EXP)
%
% EXP requires:
%  .dir_base   = [dir0,mrtID{i}];
%  .dir_voi
% (.dir_psy)
%  .filenames  = {[dir0,mrtID{i},'/swudata.nii']};
%  .voi
%   .name
%   .coord   (mm)
%   .radius  (mm)
%  .ppi
%   .NAMES
%   .CNTRSTMTX  [i j w];
%       first column: i-th condition in SPM.Sess.U(i)
%       second column: j-th parameter? in SPM.Sess.U(i).name{j} ('SOUND', 'SOUNDxLIKE')
%       third column: weight of w
% or
%   .cntrstVec   only contrast weight
%
% for myspm_fmriglm, EXP requires:
%  .TR         (sec)
%  .StimDurSec (sec)
%
%
% shape of the VOI, name, and other options
%
% (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

%% 0. parsing & check inputs

%global overwrite; if isempty(overwrite), overwrite=0; end
if isfield(EXP,'overwrite')&&EXP.overwrite, overwrite=1; else overwrite=0; end
if ~isfield(EXP,'dir_voi')&&isfield(EXP,'dir_psy')
  EXP.dir_voi = EXP.dir_psy;
else
  EXP.dir_psy = EXP.dir_voi;
end
for n=1:numel(EXP.filenames)
  if ~strcmp(EXP.filenames{n}(1),'/') % relative path?
    EXP.filenames{n} = fullfile(EXP.dir_base, EXP.filenames{n});
  end
  fnames{n,1} = [EXP.filenames{n},',1'];
end
try ls(fnames{n,1}(1:end-2));
catch ME
  error(['file ',fnames{n,1}(1:end-2),' not found']);
end
EXP.NumSess = numel(fnames);
for j=1:EXP.NumSess
  hdr   = load_untouch_header_only(fnames{j}(1:end-2));
  NF(j) = hdr.dime.dim(5);
end
EXP.NumFrames=NF;

%% isotropic smoothing
if isfield(EXP,'fwhm')
  EXP.fnames = fnames;
  EXP = myspm_smooth(EXP);
  fnames = EXP.fnames;
end

%% 1. create VOI (physiological factor)

matlabbatch={};
matlabbatch{1}.spm.util.voi.spmmat  = {[EXP.dir_voi,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust  = 0;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name    = [EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm'];
if isfield(EXP.voi,'local') % shift to nearest local maxima in each individual
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat   = {[EXP.dir_voi,'/SPM.mat']};
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = EXP.voi.cntrstVec;
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
  %   matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
  %   matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = EXP.voi.coord;
  %   matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = EXP.voi.radius;
  %   matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm  = 1;  % this is index of .util.voi.roi
  %   matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.mask ='';  % and what's this?
  %   matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
  %   % the updated coordinate is stored in VOI_*_1.mat, xY.xyz
else % fixed coordinate
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = EXP.voi.coord;
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = EXP.voi.radius;
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
  matlabbatch{1}.spm.util.voi.expression = 'i1';
end

EXP.fname_voi=[EXP.dir_voi,'/VOI_',EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm_1.mat'];
if ~exist(EXP.fname_voi,'file') || overwrite
  spm_jobman('initcfg');
  save([EXP.dir_voi,'/mb_voi_pca.mat']);
  spm_jobman('run',matlabbatch);
end

%% 2. compute psychological factor (contrast),
%  & product (PPI) of the neural response of physiological factor (deconvoluted BOLD) and psy

num_reg=numel(EXP.ppi.NAMES);
fname_ppi={};
for k=1:num_reg
  % first column:  i in SPM.Sess.U(i)
  % second column: j in SPM.Sess.U(i).name{j}
  % third column: weight
  if ~isfield(EXP,'dir_ppi')
    dir_ppi = [EXP.dir_base,'/ppi_',EXP.ppi.dirname];
  end
  [~,~] = mkdir(dir_ppi);
  matlabbatch={};
  matlabbatch{1}.spm.stats.ppi.spmmat = {[EXP.dir_psy,'/SPM.mat']};
  matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {EXP.fname_voi};
  EXP.ppi.CNTRSTMTX{k}
  matlabbatch{1}.spm.stats.ppi.type.ppi.u   = EXP.ppi.CNTRSTMTX{k};
  matlabbatch{1}.spm.stats.ppi.name = EXP.ppi.NAMES{k}; %'STG-LxBCFD';
  matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
  fname_ppi{k}=[dir_ppi,'/PPI_',EXP.ppi.NAMES{k},'.mat'];
  if ~exist(fname_ppi{k},'file') || overwrite
    spm_jobman('initcfg')
    save([dir_ppi,'/mb_ppi_deconv.mat'],'matlabbatch');
    spm_jobman('run',matlabbatch);
    fname_output=[EXP.dir_psy,'/PPI_',EXP.ppi.NAMES{k},'.mat'];
    copyfile(fname_output,fname_ppi{k})
  end
end
EXP.fname_ppi = fname_ppi;

%% 3. combine regressors [phy, psy1*phy, psy2*phy, psy3*phy, ... ]

load(EXP.fname_ppi{1}, 'PPI');
EXP.reg(1).name = 'PHY';
EXP.reg(1).val  = PPI.Y;
EXP.cntrstMtx = 0;
num_reg  = numel(EXP.fname_ppi);
PSYNAMES = EXP.ppi.NAMES;
i=2;
for k=1:num_reg
  load(EXP.fname_ppi{k}, 'PPI');
  EXP.reg(i).name = PSYNAMES{k};
  EXP.reg(i).val  = PPI.P;
  EXP.reg(i+1).name = [PSYNAMES{k},'xPHY'];
  EXP.reg(i+1).val  = PPI.ppi;
  i=i+2;
  EXP.cntrstMtx=[EXP.cntrstMtx 0 EXP.PPICntrstVec(k)];
end
EXP.cntrstMtx=[EXP.cntrstMtx 0]; % for interceptor
EXP.cntrstMtx=[EXP.cntrstMtx; -EXP.cntrstMtx];

%% 4. run fmriglm!

EXP.dir_glm = dir_ppi;
EXP = myspm_fmriglm(EXP);
