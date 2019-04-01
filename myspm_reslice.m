function EXP = myspm_reslice(EXP)
% EXP requires
%  .fname_ref
%  .fname_src
% (.interp)
% (.prefix)
%
% (cc) 2017, sgKIM
if ~nargin, help(mfilename); return; end
if isfield(EXP,'interp')
 interp=EXP.interp;
else
 interp=4;
end
matlabbatch={};
matlabbatch{1}.spm.spatial.coreg.write.ref = {[EXP.fname_ref,',1']};
matlabbatch{1}.spm.spatial.coreg.write.source = cellstr(spm_select('expand',EXP.fname_src));
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = interp;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
if isfield(EXP,'prefix')
 matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = EXP.prefix;
end

spm_jobman('initcfg')
spm_jobman('run', matlabbatch)

end