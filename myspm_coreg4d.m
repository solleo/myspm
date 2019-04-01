function EXP = myspm_coreg4d(EXP)
% EXP = myspm_coreg(EXP)
%
% EXP requires:
%  .prefix
%  .interp
%  .fname_moving
%  .fname_fixed
%
% (cc) 2015, sgKIM.
if ~nargin,  help myspm_coreg;  return; end
if ~isfield(EXP,'interp'), EXP.interp=4; end
if ~isfield(EXP,'overwrite'); EXP.overwrite=0; end
if ~isfield(EXP,'prefix'), EXP.prefix='o'; end
fnames_moving{1,1} = [EXP.fname_moving,',1'];
ls(EXP.fname_moving);
ls(EXP.fname_fixed);
[p1,f1,e1]=fileparts(EXP.fname_moving);
if isempty(p1), p1=pwd; end
cd(p1);
% hdr = load_untouch_header_only(EXP.fname_moving);
% NT = hdr.dime.dim(5);
% fnames_others = cell(NT,1);
% for t=1:NT
%  fnames_others{t,1}=[EXP.fname_moving,',',num2str(t)];
% end
fnames_others=cellstr(spm_select('expand',EXP.fname_moving));
fname_out=fullfile(p1,[EXP.prefix,f1,e1]);
if ~exist(fname_out,'file') || EXP.overwrite
 matlabbatch={};
 matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = {[EXP.fname_fixed,',1']};
 matlabbatch{1}.spm.spatial.coreg.estwrite.source = fnames_moving;
 matlabbatch{1}.spm.spatial.coreg.estwrite.other  = fnames_others;
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol ...
  = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
 matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm   = [7 7];
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = EXP.interp;
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap   = [0 0 0];
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask   = 0;
 matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = EXP.prefix;
 save([p1,'/coreg_',datestr(now,'yyyymmdd'),'.mat'], 'matlabbatch');
 spm('Defaults','fmri')
 spm_jobman('initcfg');
 spm_jobman('run', matlabbatch);
end
end
