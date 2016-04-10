function EXP = myspm_norm1(EXP)
% EXP = myspm_norm1(EXP)
%
% EXP requires:
%
% -for single subject:
%  .fname_moving     fullpath
%  .isAsian
%
% (cc) 2015, sgKIM

spm('Defaults','fmri');
% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12'),  error(['Run ',mfilename,' on SPM12!']);  end
if ~isfield(EXP,'interp'), EXP.interp = 1; end % linear as default

if ~strcmp(EXP.fname_moving(1),'/')
  EXP.fname_moving=[pwd,'/',EXP.fname_moving];
end

matlabbatch={};
matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {[EXP.fname_moving,',1']};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii'};
matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
if isfield(EXP,'eastern') && EXP.eastern
  matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'eastern';
end
matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

end
