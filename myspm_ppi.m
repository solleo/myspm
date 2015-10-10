function EXP = myspm_ppi (EXP)
% EXP = myspm_ppi (EXP)
%
% requires:
% EXP.dir_glm
%
% EXP
%  .voi
%   .name
%   .coord
%   .radius
%   .fixed
%  .ppi
%   .name
%   .cntrstVec
%
% EXP.TR = 1;
% EXP.StimDurSec = 30;
% EXP.dir_base   = [dir0,mrtID{i}];
% EXP.filenames  = {[dir0,mrtID{i},'/swudata.nii']};
% EXP.fname_rp   = [dir0,mrtID{i},'/rp_data.txt'];
% EXP.model_desc = 'FC+FD+BC+BD_RP';
%
% shape of the VOI, name, and other options
%
% (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

%% 0. parsing & check inputs
spm('Defaults','fmri');
spm_jobman('initcfg');

model_desc = EXP.model_desc;
if ~isfield(EXP,'dir_base')
  EXP.dir_base=pwd;
end
EXP.dir_glm = [EXP.dir_base,'/glm_',model_desc];
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

try ls(fnames{n,1}(1:end-2));
catch ME
  fnames
  error('file not found');
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

%% 1. compute PPI
spm('Defaults','fmri');
spm_jobman('initcfg');

% 1-1. compute Physioloigcal vector
matlabbatch={};
matlabbatch{1}.spm.util.voi.spmmat = {[EXP.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust = 0;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name = [EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm']; %'STG-left';
if isfield(EXP.voi,'local')
  matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = EXP.voi.local.spmindex;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
  matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
  
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = EXP.voi.coord; % [54 -2 -10]
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = EXP.voi.radius;
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm  = 1;  % this could be the contrast index...
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.mask ='';  % and what's this?
  matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
  % the updated coordinate is stored in VOI_*_1.mat, xY.xyz
else
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = EXP.voi.coord; % [54 -2 -10]
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = EXP.voi.radius;
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
  matlabbatch{1}.spm.util.voi.expression = 'i1';
end

EXP.fname_voi=[EXP.dir_glm,'/VOI_',EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm_1.mat'];
if ~exist(EXP.fname_voi,'file')
  save([EXP.dir_glm,'/mb_voi_pca.mat']);
  spm_jobman('run',matlabbatch);
end

%% 2. compute Psychological vector,
%  neural response of Physiological vector (deconvoluted BOLD),
%  and a product of them
NumCond = numel(EXP.ppi.cntrstVec);
EXP.ppi.cntrstMtx=[1:NumCond; ones(1,NumCond); EXP.ppi.cntrstVec]';
[~,~] = mkdir([EXP.dir_base,'/glm_PPI_',EXP.ppi.name]);
EXP.fname_ppi = [EXP.dir_base,'/glm_PPI_',EXP.ppi.name,'/PPI.mat'];

if EXP.TR<6 && ~isfield(EXP,'noDeconv') % 6 seconds of time-lag in BOLD
  matlabbatch={};
  matlabbatch{1}.spm.stats.ppi.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {EXP.fname_voi};
  matlabbatch{1}.spm.stats.ppi.type.ppi.u = EXP.ppi.cntrstMtx;
  matlabbatch{1}.spm.stats.ppi.name = EXP.ppi.name; %'STG-LxBCFD';
  matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
  spm_jobman('initcfg')
  if ~exist(EXP.fname_ppi,'file')
    save([EXP.dir_glm,'/mb_ppi_deconv.mat']);
    spm_jobman('run',matlabbatch);
    src=[EXP.dir_glm,'/PPI_',EXP.ppi.name,'.mat'];
    movefile(src,EXP.fname_ppi);
  end
else % just read psy (the design & contrast); read phy, and psyphy
  load ([EXP.dir_glm,'/SPM.mat']);
  load (EXP.fname_voi);
  PPI=[];
  PPI.P = SPM.xX.X(:,1:numel(EXP.ppi.cntrstVec)) * EXP.ppi.cntrstVec'; % if you used FIR for sparse sampling, H*X is also fine.
  PPI.Y = Y;
  PPI.ppi = PPI.P .* PPI.Y;
  PPI.xY = xY;
  save (EXP.fname_ppi, 'PPI');
  figure;
  subplot(311); plot(PPI.Y); ylabel('Phy');
  subplot(312); plot(PPI.P); ylabel('Psy');
  subplot(313); plot(PPI.ppi); ylabel('PsyPhy');
  screen2png([EXP.dir_glm,'/PPI.png'],120);
  close(gcf);
  clear SPM
end

%% 3. run fmri-glm with the PPI
EXP.model_desc = ['PPI_',EXP.ppi.name];
EXP = myspm_fmriglm(EXP);

end