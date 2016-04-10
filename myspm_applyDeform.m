function myspm_applyDeform(EXP)
% function myspm_applyDeform(EXP)

% '/scr/vatikan3/Myelin2/spm12/y_t1w_brain_2001.nii'
% '/scr/vatikan3/Myelin2/spm12/qR1_2001.nii'

spm('Defaults','fmri');

matlabbatch={};
matlabbatch{1}.spm.util.defs.comp{1}.def = {EXP.deform};
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {EXP.moving};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];

spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end