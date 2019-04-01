function myspm_apply_deform(EXP)
% EXP = myspm_apply_deform(EXP)
%
% EXP requires:
% .fname_deform
% .isInv
% .fname_moving
% .fname_fixed
% .intorder
%
% (cc) 2017, sgKIM.

if ~isfield(EXP,'intorder'), intorder=1; else, intorder=EXP.intorder; end
if ~isfield(EXP,'isInv'), isInv=0; else, isInv=EXP.isInv; end

matlabbatch={};
if isInv
 matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {EXP.fname_deform};
 matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {EXP.fname_fixed};
else
 matlabbatch{1}.spm.util.defs.comp{1}.def = {EXP.fname_deform};
end
matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {EXP.fname_moving};
matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.savepwd = 1;
matlabbatch{1}.spm.util.defs.out{1}.pull.interp = intorder;
matlabbatch{1}.spm.util.defs.out{1}.pull.mask   = 0;
matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm   = [0 0 0];
matlabbatch{1}.spm.util.defs.out{1}.prefix = 'd';

spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
end
