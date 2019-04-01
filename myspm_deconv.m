function EXP = myspm_deconv (EXP)
%{
EXP = myspm_deconv (EXP)

requires:
EXP
.dir_base
.dir_phi  (default: dir_base)

.phi
.name
.fname_voi | .fname_mask | .cluster | .sphere (either one of them)
.fname_voi (for already computed VOI .mat file)
or
.fname_mask (for a mask image .nii or .img file)
or
.cluster (for a cluster ID maps (sigclus*.nii) created from myspm_fmriglm)
.index
.sign (either '+' or '-')
or
.sphere (spherical VOI)
.centre (mm)
.radius (mm)

.psy
.name
.cntrstVec (psychological contrast as in GLM)

(.ppi.name) default: [phi.name,'x',psy.name]

OUTPUT:

Time of the first sample is SPM.xBF.T0-th mirco-time sample
If dt is 0.0625 (or 1/16) sec as usual, the T(0) is 0.5 sec.

See:
RT=SPM.xY.RT
dt=SPM.xBF.dt
fMRI_T0=SPM.xBF.T0

(cc) 2015, 2016, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com
%}

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

%if ~isfield(EXP,'dir_phi'), EXP.dir_phi = EXP.dir_base; end
% if ~isfield(EXP,'dir_psy'), EXP.dir_psy = EXP.dir_base; end
[p1,~,~]=fileparts(EXP.dir_glm);
EXP.dir_base=p1;
EXP.dir_phi=EXP.dir_glm;

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
matlabbatch{1}.spm.util.voi.spmmat  = {[EXP.dir_phi,'/SPM.mat']};
matlabbatch{1}.spm.util.voi.adjust  = 0;
matlabbatch{1}.spm.util.voi.session = 1;
if ~isfield(EXP,'fname_voi')
if isfield(EXP.phi,'fname_mask') % just a file
matlabbatch{1}.spm.util.voi.name = EXP.phi.name;
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[EXP.phi.fname_mask,',1']};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
matlabbatch{1}.spm.util.voi.expression = 'i1';
EXP.fname_voi=[EXP.dir_phi,'/VOI_',EXP.phi.name,'_1.mat'];
elseif isfield(EXP.phi,'cluster') % create a mask file for i-th contrast, j-th cluster
matlabbatch{1}.spm.util.voi.name = ['clus_',EXP.phi.name];
if ~isfield(EXP.phi.cluster,'dir'); EXP.phi.cluster.dir=EXP.dir_phi; end
fname_sig = [EXP.phi.cluster.dir,'/sigclus_',EXP.phi.cluster.sign,'1.nii'];
nii = load_uns_nii(fname_sig);
cluster_indices=sort(unique(nii.img(:)));
cluster_indices(1)=[];
nii.img = uint8(nii.img == cluster_indices(EXP.phi.cluster.index));
nii.hdr.dime.datatype=2;
fname_cluster = [EXP.dir_phi,'/clusmask_',EXP.phi.name,'.nii'];
save_untouch_nii(nii, fname_cluster);
clear nii
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = {[fname_cluster,',1']};
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0;
matlabbatch{1}.spm.util.voi.expression = 'i1';
EXP.fname_voi = [EXP.dir_phi,'/', ...
'VOI_',matlabbatch{1}.spm.util.voi.name,'_1.mat'];
elseif isfield(EXP.phi,'sphere') % sphere
matlabbatch{1}.spm.util.voi.name = [EXP.phi.name,'_r',num2str(EXP.phi.sphere.radius),'mm'];
matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = EXP.phi.sphere.centre;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = EXP.phi.sphere.radius;
matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
matlabbatch{1}.spm.util.voi.expression = 'i1';
EXP.fname_voi=[EXP.dir_phi,'/', ...
'VOI_',EXP.phi.name,'_r',num2str(EXP.phi.sphere.radius),'mm_1.mat'];
end
% % in case of Y is not residual
%  if isfield(EXP,'voiadjust') && EXP.phiadjust
%   matlabbatch{1}.spm.util.voi.adjust = NaN; % This seems like taking eig1 from residuals... is it a right way to do?
%  end
end
if ~exist(EXP.fname_voi,'file') || overwrite
save([EXP.dir_phi,'/myspm_ppi_1_create_voi.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);
end

%% 2. do simple deconvolution (sd)
matlabbatch={};
matlabbatch{1}.spm.stats.ppi.spmmat = {[EXP.dir_phi,'/SPM.mat']};
matlabbatch{1}.spm.stats.ppi.type.sd.voi = {EXP.fname_voi};
matlabbatch{1}.spm.stats.ppi.name=EXP.phi.name;
matlabbatch{1}.spm.stats.ppi.disp = 0; % not very informative
EXP.fname_ppi=[EXP.dir_phi,'/PPI_',EXP.phi.name,'.mat'];

spm_jobman('initcfg')
if ~exist(EXP.fname_ppi,'file') || overwrite
save([EXP.dir_phi,'/myspm_ppi_2_construct_ppi.mat'],'matlabbatch');
spm_jobman('run',matlabbatch);
%copyfile([EXP.dir_phi,'/PPI_',EXP.ppi.name,'.mat'],dir_ppi);
end

end
