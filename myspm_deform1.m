function myspm_deform1(EXP)
% EXP = myspm_deform1(EXP)
%
% EXP requires:
%  .fname_moving
%  .fname_def
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

fnames={};
if ~iscell(EXP.fname_moving)
  EXP.fname_moving={EXP.fname_moving};
end
j=1;
for c=1:numel(EXP.fname_moving) % for files
  hdr = load_untouch_header_only(EXP.fname_moving{c});
  for t=1:hdr.dime.dim(5) % for timepoints
    fnames{j} = [EXP.fname_moving{c},',',num2str(t)];
    j=j+1;
  end
end

normalise=[];
normalise.write.subj.def = {EXP.fname_def};
normalise.write.subj.resample = fnames;
normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
if ~isfield(EXP,'vox_mm'), EXP.vox_mm=1; end
normalise.write.woptions.vox = [1 1 1]*EXP.vox_mm;
normalise.write.woptions.interp = EXP.interp; 

matlabbatch={};
matlabbatch{1}.spm.spatial.normalise = normalise;

% save([path1,'/deform1.mat'], 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end