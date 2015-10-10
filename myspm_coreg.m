function EXP = myspm_coreg(EXP)
% EXP = myspm_coreg(EXP)
%
% (cc) 2015, sgKIM.

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
[path1,~,~]= fileparts(EXP.name_moving);
if isempty(path1), path1=pwd; end
cd(path1);

if isfield(EXP,'name_others')
  fnames={};
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
if ~exist(fullfile(path1,[EXP.prefix,EXP.name_moving]),'file')
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.estwrite.ref    = {[EXP.name_fixed,',1']};
  matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[EXP.name_moving,',1']};
  matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol ...
    = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
  matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = EXP.interp;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
  matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = EXP.prefix;
  
  save([path1,'/coreg.mat'], 'matlabbatch');
  spm_jobman('run', matlabbatch)
end
if isfield(EXP,'name_others')
  matlabbatch={};
  matlabbatch{1}.spm.spatial.coreg.write.ref =  {[EXP.name_fixed,',1']};
  matlabbatch{1}.spm.spatial.coreg.write.source = fnames;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = EXP.interp;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
  matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
  matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = EXP.prefix;
end
save([path1,'/coreg_others.mat'], 'matlabbatch');
spm_jobman('run', matlabbatch)
end