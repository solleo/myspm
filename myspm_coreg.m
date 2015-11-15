function EXP = myspm_coreg(EXP)
% EXP = myspm_coreg(EXP)
%
% .prefix
% .interp
% .name_moving
% .name_fixed
% .name_others
%
% (cc) 2015, sgKIM.
if ~nargin,  help myspm_coreg;  return; end
[~,myname] = fileparts(mfilename('fullpath'));
disp(['### ',myname,': starting..']);

if ~isfield(EXP,'subjID')
  myspm_coreg1(EXP);
else
  path0=pwd;
  subjID = fsss_subjID(EXP.subjID);
  for n=1:numel(subjID)
    subjid=subjID{n};
    path1=fullfile(EXP.dir_base,subjid);
    cd(path1)
    myspm_coreg1(EXP);
  end
  cd(path0);
end
end

function myspm_coreg1(EXP)
if ~isfield(EXP,'interp'), EXP.interp=1; end
spm('Defaults','fmri')
spm_jobman('initcfg');
if ~isfield(EXP,'prefix'), EXP.prefix='o'; end
[path1,name1,ext1]= fileparts(EXP.name_moving);
if isempty(path1), path1=pwd; end
cd(path1);

fnames={''};
if isfield(EXP,'name_others')
  if iscell(EXP.name_others)
    j=1;
    for c=1:numel(EXP.name_others)
      hdr = load_untouch_header_only(EXP.name_others{c});
      for t=1:hdr.dime.dim(5)
        fnames{j} = [EXP.name_others{c},',',num2str(t)];
        j=j+1;
      end
    end
  else
    j=1;
    hdr = load_untouch_header_only(EXP.name_others);
    for t=1:hdr.dime.dim(5)
      fnames{j} = [EXP.name_others,',',num2str(t)];
      j=j+1;
    end
  end
end

matlabbatch={};
matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = {[EXP.name_fixed,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[EXP.name_moving,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other  = fnames;
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol ...
  = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = EXP.interp;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = EXP.prefix;

save(fullfile(path1,'coreg.mat'), 'matlabbatch');
spm_jobman('run', matlabbatch);

%sometimes it includes NaN (but for label iamge?)
if isfield(EXP,'name_others')
  for j=1:numel(fnames)
    [path1,name1,~] = fileparts(fnames{j});
    nii = load_untouch_nii(fullfile(path1,[EXP.prefix,name1,'.nii']));
    nii.img(isnan(nii.img)) = 0;
    save_untouch_nii(nii,  fullfile(path1,[EXP.prefix,name1,'.nii']));
  end
end



end