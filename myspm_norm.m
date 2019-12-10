function JOB = myspm_norm(JOB)
% JOB = myspm_norm(JOB)
%
% JOB requires:
%  .fname_moving
%  .fname_deform (y_*)
% (.vox_mm)
% (.interp) B-spline order from 0 to 6 (default=4)
%
% (cc) 2015, sgKIM, solleo@gmail.com
if nargin==0, help(mfilename); return; end
spm('Defaults','fmri');
% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12'),  error(['Run ',mfilename,' on SPM12!']);  end
if ~isfield(JOB,'interp'), JOB.interp = 4; end % 4-spline as default

fnames={};
if iscell(JOB.fname_moving)
 [path1,~,~]= fileparts(JOB.fname_moving{1});
 j=1;
 for c=1:numel(JOB.fname_moving)
  hdr = load_untouch_header_only(JOB.fname_moving{c});
  for t=1:hdr.dime.dim(5)
   fnames{j,1} = [JOB.fname_moving{c},',',num2str(t)];
   j=j+1;
  end
 end
else
 [path1,~,~]= fileparts(JOB.fname_moving);
 j=1;
 hdr = load_untouch_header_only(JOB.fname_moving);
 for t=1:hdr.dime.dim(5)
  fnames{j,1} = [JOB.fname_moving,',',num2str(t)];
  j=j+1;
 end
end

normalise=[];
normalise.write.subj.def = {JOB.fname_deform};
normalise.write.subj.resample = fnames;
normalise.write.woptions.bb = [-78 -112 -70; 78 76 85];
if ~isfield(JOB,'vox_mm')
 normalise.write.woptions.vox = hdr.dime.pixdim(2:4);
else
 normalise.write.woptions.vox = JOB.vox_mm;
end
normalise.write.woptions.interp = JOB.interp;

matlabbatch={};
matlabbatch{1}.spm.spatial.normalise = normalise;

save(fullfile(path1,['spm12_norm.mat']), 'matlabbatch');
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end
