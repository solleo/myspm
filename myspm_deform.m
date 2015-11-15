function EXP = myspm_deform(EXP)
% EXP = myspm_deform(EXP)
%
% EXP requires:
%
% -for single subject:
%  .fname_moving     fullpath
%  .fname_def
% (.prefix
% (.interp) [1x1] default=1 (linear)
% (.vox_mm) [1x1] default=2
%
% -for multiple subjects:
%  .subjID
%  .dir_base
%  .name_moving
%  .name_def
%
% -optionally:
% (.prefix)  output prefix (default:'w')
%
%
% (cc) 2015, sgKIM

path0=pwd;
if ~isfield(EXP,'subjID')
  myspm_deform1(EXP);
else
  subjID = fsss_subjID(EXP.subjID);
  for n=1:numel(subjID)
    subjid=subjID{n};
    path1=fullfile(EXP.dir_base,subjid);
    cd(path1)
    
    exp1=EXP;
    exp1.name_moving=[];
    error('now code batch process!');
    myspm_deform1(exp1);    
  end
  cd(path0);
end
end

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

fnames={};
if iscell(EXP.fname_moving)
  [path1,name1,ext1]= fileparts(EXP.fname_moving{1});
  j=1;
  for c=1:numel(EXP.fname_moving)
    hdr = load_untouch_header_only(EXP.fname_moving{c});
    for t=1:hdr.dime.dim(5)
      fnames{j} = [EXP.fname_moving{c},',',num2str(t)];
      j=j+1;
    end
  end
else
  [path1,name1,ext1]= fileparts(EXP.fname_moving);
  j=1;
  hdr = load_untouch_header_only(EXP.fname_moving);
  for t=1:hdr.dime.dim(5)
    fnames{j} = [EXP.fname_moving,',',num2str(t)];
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

save([path1,'/deform1.mat'], 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end