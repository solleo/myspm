function EXP = myspm_ppi (EXP)
% EXP = myspm_ppi (EXP)
%
% requires:
% EXP
% .dir_glm
% .dir_voi
% .dir_psy
%
% EXP
%  .voi
%   .name
%   .coord
%   .radius
% or
%   .local
%   .cntrstVec
%  .ppi
%   .name
%   .cntrstVec     % for parametric design: [0 0 2] means 3nd condition and 2nd parameter
%  .redovoi
%
% EXP.TR = 1;
% EXP.StimDurSec = 30;
% EXP.dir_base   = [dir0,mrtID{i}];
% EXP.filenames  = {[dir0,mrtID{i},'/swudata.nii']};
% EXP.fname_rp   = [dir0,mrtID{i},'/rp_data.txt'];
%
% shape of the VOI, name, and other options
%
% (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

%% 0. parsing & check inputs

global overwrite; if isempty(overwrite), overwrite=0; end
spm('Defaults','fmri');
spm_jobman('initcfg');
for n=1:numel(EXP.filenames)
  if ~strcmp(EXP.filenames{n}(1),'/') % relative path?
    EXP.filenames{n} = fullfile(EXP.dir_base, EXP.filenames{n});
  end
  fnames{n,1} = [EXP.filenames{n},',1'];
end
try ls(fnames{n,1}(1:end-2));
catch ME
  fnames
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

%% 1. create VOI (for physiological vector)
spm('Defaults','fmri');
spm_jobman('initcfg');

matlabbatch={};
matlabbatch{1}.spm.util.voi.spmmat  = {[EXP.dir_voi,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust  = 0;
matlabbatch{1}.spm.util.voi.session = 1;
matlabbatch{1}.spm.util.voi.name    = [EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm']; %'STG-left';
if isfield(EXP.voi,'local')
  matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat   = {[EXP.dir_voi,'/SPM.mat']};
  matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = EXP.voi.cntrstVec;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none';
  matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.001;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0;
  matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});
  
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = EXP.voi.coord; % [54 -2 -10]
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = EXP.voi.radius;
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.spm  = 1;  % this is index of .util.voi.roi
  matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.local.mask ='';  % and what's this?
  matlabbatch{1}.spm.util.voi.expression = 'i1&i2';
  % the updated coordinate is stored in VOI_*_1.mat, xY.xyz
else
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = EXP.voi.coord; % [54 -2 -10]
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = EXP.voi.radius;
  matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
  matlabbatch{1}.spm.util.voi.expression = 'i1';
end

EXP.fname_voi=[EXP.dir_voi,'/VOI_',EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm_1.mat'];
if ~exist(EXP.fname_voi,'file') || isfield(EXP,'redovoi') || overwrite
  save([EXP.dir_voi,'/mb_voi_pca.mat']);
  spm_jobman('run',matlabbatch);
end

%% 2. compute Psychological vector,
%  neural response of Physiological vector (deconvoluted BOLD),
%  and a product of them (ppi)

% first column: i in SPM.Sess.U(i)
% second column: j in SPM.Sess.U(i).name{j}
% third column: weight
if ~isfield(EXP.ppi,'cntrstMtx')
  NumCond = numel(EXP.ppi.cntrstVec);
  EXP.ppi.cntrstMtx = [1:NumCond; ones(1,NumCond); EXP.ppi.cntrstVec;]';
end

if ~isfield(EXP,'dir_ppi')
  dir_ppi = [EXP.dir_base,'/ppi_',EXP.ppi.name];
end
[~,~] = mkdir(dir_ppi);
EXP.fname_ppi = [dir_ppi,'/PPI.mat'];

%if EXP.TR<6 && ~isfield(EXP,'noDeconv') % less than 6 seconds of time-lag in BOLD
  matlabbatch={};
  matlabbatch{1}.spm.stats.ppi.spmmat = {[EXP.dir_psy,'/SPM.mat']};
  matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {EXP.fname_voi};
  matlabbatch{1}.spm.stats.ppi.type.ppi.u = EXP.ppi.cntrstMtx;
  matlabbatch{1}.spm.stats.ppi.name = EXP.ppi.name; %'STG-LxBCFD';
  matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
  spm_jobman('initcfg')
  if ~exist(EXP.fname_ppi,'file') || overwrite
    save([dir_ppi,'/mb_ppi_deconv.mat'],'matlabbatch');
    spm_jobman('run',matlabbatch);
    src=[EXP.dir_psy,'/PPI_',EXP.ppi.name,'.mat'];
    movefile(src,EXP.fname_ppi);
  end
% else % just read psy (the design & contrast); read phy, and psyphy
%   load ([EXP.dir_voi,'/SPM.mat']);
%   load (EXP.fname_voi);
%   PPI=[];
%   PPI.P = SPM.xX.X(:,1:numel(EXP.ppi.cntrstVec)) * EXP.ppi.cntrstVec'; % if you use FIR for sparse sampling, H*X is also fine.
%   PPI.Y = Y;
%   PPI.ppi = PPI.P .* PPI.Y;
%   PPI.xY = xY;
%   save (EXP.fname_ppi, 'PPI');
%   figure;
%   subplot(311); plot(PPI.Y); ylabel('Phy');
%   subplot(312); plot(PPI.P); ylabel('Psy');
%   subplot(313); plot(PPI.ppi); ylabel('PsyPhy');
%   screen2png([EXP.dir_voi,'/PPI.png'],120);
%   close(gcf);
%   clear SPM
% end

%% 3. run fmri-glm with the PPI
EXP.model_desc = ['PPI_',EXP.ppi.name];
EXP.dir_glm = dir_ppi;
EXP.dir_ppi = dir_ppi;
EXP = myspm_fmriglm(EXP);

end