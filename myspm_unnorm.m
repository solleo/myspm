function myspm_unnorm(EXP)
% myspm_unnorm(EXP)
%
% EXP requires:
%  .fname_y        : T1w deformation field (e.g. y_t1w.nii)
%  .fname_ref      : EPI coregistered in native-T1w space (e.g., uaEPI.nii)
%  .fname_input    : Image in MNI152 space
% (.prefix)        : (default=x)
% (.interp)        : (default=4; 4-th order B-spline)

defs=[];
ls(EXP.fname_y)
defs.comp{1}.inv.comp{1}.def = {EXP.fname_y};
ls(EXP.fname_ref)
defs.comp{1}.inv.space = {EXP.fname_ref};
ls(EXP.fname_input)
hdr=load_untouch_header_only(EXP.fname_input);
% T=hdr.dime.dim(5);
% for t=1:T
% defs.out{1}.pull.fnames{t,1} = [EXP.fname_input,',',num2str(t)];
% end
defs.out{1}.pull.fnames{1}=EXP.fname_input;
defs.out{1}.pull.savedir.savesrc = 1;
if ~isfield(EXP,'interp'), EXP.interp=4; end
defs.out{1}.pull.interp = EXP.interp;
defs.out{1}.pull.mask = 1;
defs.out{1}.pull.fwhm = [0 0 0];
if ~isfield(EXP,'prefix'), EXP.prefix='x'; end
defs.out{1}.pull.prefix = EXP.prefix;
matlabbatch={};
matlabbatch{1}.spm.util.defs = defs;
spm_jobman('run',matlabbatch)
[p1,f1,e1]=fileparts(EXP.fname_input);
fname_out=[p1,'/',EXP.prefix,f1,e1];
ls(fname_out)
myspm_denan(fname_out);
end