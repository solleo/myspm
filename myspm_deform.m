function EXP = myspm_deform(EXP)
% EXP = myspm_deform(EXP)
%
% EXP requires:
%
% -for single subject:
%  .fname_moving     fullpath
%  .fname_def
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
    
      
%   exp1=[];
%   exp1.name_moving=[dir_base,subjid,'.probtrkX/',trkname];
%   exp1.name_def=[dir_rest,subjid,'/y_t1w.nii'];
%   exp1.prefix='w';
%   exp1.vox_mm=1.5;
    
    
  end
  cd(path0);
end
end

function myspm_deform1(EXP)
% EXP = myspm_deform1(EXP)
%
% EXP requires:
%  .name_moving
%  .name_def
%  .name_others
%
% (cc) 2015, sgKIM

spm('Defaults','fmri');
% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12'),  error(['Run ',mfilename,' on SPM12!']);  end
if ~isfield(EXP,'prefix'), EXP.prefix='w'; end
if ~isfield(EXP,'interp'), EXP.interp=1; end % linear as default

[path1,name1,ext1]= fileparts(EXP.fname_moving);

if isfield(EXP,'fname_moving')
  fnames={};
  if iscell(EXP.fname_moving)
    j=1;
    for c=1:numel(EXP.fname_moving)
      hdr = load_untouch_header_only(EXP.fname_others{c});
      for t=1:hdr.dime.dim(5)
        fnames{j} = [EXP.fname_moving{c},',',num2str(t)];
        j=j+1;
      end
    end
  else
    j=1;
    hdr = load_untouch_header_only(EXP.fname_moving);
    for t=1:hdr.dime.dim(5)
      fnames{j} = [EXP.fname_moving,',',num2str(t)];
      j=j+1;
    end
  end
end

% applying deformation to 'other' image?
normalise=[];
normalise.write.subj.def(1) = {EXP.fname_def};
normalise.write.subj.resample(1) = fnames;
normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
normalise.write.woptions.vox = [1 1 1]*EXP.vox_mm;
normalise.write.woptions.interp = EXP.interp; 

matlabbatch={};
matlabbatch{1}.spm.spatial.normalise = normalise;

save([path1,'/deform1.mat'], 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

% remove negative values from B-splines! (???)
% if EXP.interp > 2
%   [path1,name1,ext1]=fileparts(EXP.fname_moving);
%   nii = load_untouch_nii(EXP.fname_moving);
%   if ~sum(nii.img(:)<0) % no negative values in the original image
%     nii = load_untouch_nii([path1,'/',EXP.prefix,name1,ext1]);
%     nii.img(nii.img<0) = 0;
%     save_untouch_nii(nii,  [path1,'/',EXP.prefix,name1,ext1]);
%   end
% end

end