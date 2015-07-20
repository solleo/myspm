function EXP=myspm_seg12(EXP)
% EXP.mrf (Markov Random Field strengh
% EXP.cleanup (0=no, 1=light, 2=thorough

for c=1:numel(EXP.fname_scans)
  matlabbatch{1}.spm.spatial.preproc.channel(c).vols = {[EXP.fname_scans{c},',1']};
  matlabbatch{1}.spm.spatial.preproc.channel(c).biasreg = 0.001;
  matlabbatch{1}.spm.spatial.preproc.channel(c).biasfwhm = 60;
  matlabbatch{1}.spm.spatial.preproc.channel(c).write = [0 0];
end

matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/afs/cbs.mpg.de/software/spm/12.6225/8.2/2.15/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = EXP.mrf;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = EXP.cleanup;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni'; %european ICBM
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
end