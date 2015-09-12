function EXP = myspm_ppi (EXP)
% EXP = myspm_ppi (EXP)
%
% requires:
% EXP.dir_glm
% EXP.voi.name
% EXP.voi.coord
% EXP.ppi.name
% EXP.ppi.cntrstVec
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
  EXP.dir_glm = [pwd,'/glm_',model_desc];
else
  EXP.dir_glm = [EXP.dir_base,'/glm_',model_desc];
end
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

ls(fnames{n,1}(1:end-2));

EXP.NumSess = numel(fnames);
for j=1:EXP.NumSess
  hdr   = load_untouch_header_only(fnames{j}(1:end-2));
  NF(j) = hdr.dime.dim(5);
end
EXP.NumFrames=NF;

%% 0.5. isotropic smoothing

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
    spm_jobman('run', matlabbatchs)
  end
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
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = EXP.voi.coord; % [54 -2 -10]
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = EXP.voi.radius;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';
EXP.fname_voi=[EXP.dir_glm,'/VOI_',EXP.voi.name,'_r',num2str(EXP.voi.radius),'mm_1.mat'];
if ~exist(EXP.fname_voi,'file')
  spm_jobman('run',matlabbatch);
end

%% 2. compute Psychological vector and neural response of Physiological
% vector, and product them
NumCond = numel(EXP.ppi.cntrstVec);
EXP.ppi.cntrstMtx=[1:NumCond; ones(1,NumCond); EXP.ppi.cntrstVec]';
matlabbatch={};
matlabbatch{1}.spm.stats.ppi.spmmat = {[EXP.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {EXP.fname_voi};
matlabbatch{1}.spm.stats.ppi.type.ppi.u = EXP.ppi.cntrstMtx;
matlabbatch{1}.spm.stats.ppi.name = EXP.ppi.name; %'STG-LxBCFD';
matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
EXP.fname_ppi  = [EXP.dir_glm,'/PPI_',EXP.ppi.name,'.mat'];
spm_jobman('initcfg')
if ~exist(EXP.fname_ppi,'file')
  spm_jobman('run',matlabbatch);
end

%% 3. run fmri-glm with the PPI
EXP.model_desc = ['PPI_',EXP.ppi.name];
EXP = myspm_fmriglm(EXP);

end