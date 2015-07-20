% NOT DONE YET!!!

function EXP = myspm_long12 (EXP)
% performes longitudinal vbm using spm12
%
% input:
% EXP.files_query: {1xTP} fresh raw NIFTIs from all subjects for TP timepoints
% EXP.tdiff: [Nx1] interval between two scans, or a constant [1x1]
%
% workflow:
% [intrasubject]
% 1. Segmentation of images (save: native c1/c2/c3 and bias-corrected images)
% 2. Combine GM/WM/CSF (c1/c2/c3) and mask the bias corrected image to create skullstripped
%     brain from (1)
% 3. Pair-wise registration of skull-stripped/bias-corrected images from (2)
%    resulting in intrasubject template, Jacobian rate, divergence rate.
%
% [intersubject]
% 4. Segmentation of avg images from (3) (preproc for Dartel; so save only Dartel import in native space)
% 5. Dartel template creation from (4)
% 6. Warping spc_mwc1 (change of GMV/WMV) into the group template from (6) and (7)
% 7. Warping jd_ (rate of JD of whole brain) into the group template from (3) and (7)
%
% [stats]
% 8. Run statistical tests on (6) and (7)
% 9. Warping results into MNI-space for a publication from (8)
%
%
% (cc) 2015, sgKIM. solleo@gmail.com


%% 0a. Set spm12 for pairwise intrasubject registration
% ver=spm('version');
% if str2double(ver(4:5)) < 12
%   %error('We need spm12 for pairwise intrasubject registration!')
%   warning('setting spm12 up..');
%   myspm_set_spm12
% end
% spm('Defaults','fmri');

%% 0b. get filenames
Ntp = numel(EXP.files_query);
scans=cell(1,Ntp);
subjID={};
for t=1:Ntp
  files = dir(EXP.files_query{t});
  [mypath,~,~] = fileparts(EXP.files_query{t});
  for n=1:numel(files)
    scans{t}{n} = [mypath,'/',files(n).name,',1'];
    [~,b,c] = fileparts(files(n).name);
    subjID{n,t} = b;
  end
end
EXP.scans = scans;
EXP.subjID = subjID;
Nsubj = numel(files);

%% 0c. setting default parameters
if ~isfield(EXP,'isAsian')
  EXP.isAsian=0;
end
EXP.spmpath=spm('dir');
if ~isfield(EXP,'tdiff')
  EXP.tdiff = ones(1, numel(scans{1}));
elseif numel(EXP.tdiff) == 1
  EXP.tdiff = EXP.tdiff * ones(1, numel(scans{1}));
end

%% [intrasubject registration]=======================================

%% 1. Segmentation using vbm8
myspm_vbm8ss([EXP.scans{:}]);

%% 2. Skullstripping of all (prefix 'B') for better registration
for t=1:Ntp
  for n=1:Nsubj
    subjid = EXP.subjID{n,t};
    system(['FSLOUTPUTTYPE=NIFTI; fslmaths ',mypath,'/c1',subjid,'.nii -add ', ...
      mypath,'/c2',subjid,'.nii -add ',mypath,'/c3',subjid,'.nii -bin -mul ', ...
      mypath,'/m',subjid,'.nii ',mypath,'/Bm',subjid,'.nii']);
  end
end

%% 3. Pair-wise registration (for longitudianl or within-subject design)
for t=1:Ntp
  files = dir(['b',EXP.files_query{t}]);
  [mypath,~,~] = fileparts(EXP.files_query{t});
  for n=1:numel(files)
    scans{t}{n} = [mypath,'/',files(n).name,',1'];
  end
end
matlabbatch={};
matlabbatch{1}.spm.tools.longit{1}.pairwise.vols1 = EXP.scans{1};
matlabbatch{1}.spm.tools.longit{1}.pairwise.vols2 = EXP.scans{2};
matlabbatch{1}.spm.tools.longit{1}.pairwise.tdif = EXP.tdiff;
matlabbatch{1}.spm.tools.longit{1}.pairwise.noise = NaN;
matlabbatch{1}.spm.tools.longit{1}.pairwise.wparam = [0 0 100 25 100];
matlabbatch{1}.spm.tools.longit{1}.pairwise.bparam = 1000000;
matlabbatch{1}.spm.tools.longit{1}.pairwise.write_avg = 1;
matlabbatch{1}.spm.tools.longit{1}.pairwise.write_jac = 1;
matlabbatch{1}.spm.tools.longit{1}.pairwise.write_div = 1;
if ~isfield(EXP,'vbm8')
  write_def=0;
else
  write_def=1;
end
matlabbatch{1}.spm.tools.longit{1}.pairwise.write_def = write_def;

%% 1. Segmentation spm12 (but the segmentation really sucks...)
matlabbatch={};
matlabbatch{1}.spm.spatial.preproc.channel.vols = [EXP.scans{:}];
matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
NGaus=[1 1 2 3 4 2];
tissue=[];
for c=1:6
  tissue(c).tpm = {[fullfile(EXP.spmpath,'tpm','TPM.nii'),',',num2str(c)]};
  tissue(c).ngaus = NGaus(c);
  tissue(c).native = [1 0];
  tissue(c).warped = [0 0];
end
tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue = tissue;
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
if EXP.isAsian
  matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'eastern';
else
  matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
end
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

end