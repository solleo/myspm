function JOB = myspm_smooth(JOB)
% JOB = myspm_smooth(JOB)
%
% JOB requires:
%  .fnames
%  .fwhm_mm
%
% (cc) 2015, sgKIM.

if ~nargin, help(mfilename); return; end

fnames=JOB.fnames;
if ~iscell(fnames), fnames={fnames}; end
fnames=reshape(fnames,[],1);
if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end

matlabbatchs={};
matlabbatchs{1}.spm.spatial.smooth.dtype = 0;
matlabbatchs{1}.spm.spatial.smooth.im = 0;
matlabbatchs{1}.spm.spatial.smooth.data = fnames;
if numel(JOB.fwhm_mm)==3 % input is [1x3] or [3x1]
 if (JOB.fwhm_mm(1) == JOB.fwhm_mm(2)) ...
   && (JOB.fwhm_mm(2) == JOB.fwhm_mm(3)) % isotropic?
  f=JOB.fwhm_mm(3);
 else % anisotropic
  f=[num2str(JOB.fwhm_mm(1)),'_',num2str(JOB.fwhm_mm(2)), ...
   '_',num2str(JOB.fwhm_mm(3))];
 end
elseif numel(JOB.fwhm_mm)==1 % input is [1x1]
 f=JOB.fwhm_mm;
 JOB.fwhm_mm = [1 1 1].*JOB.fwhm_mm;
else
 error(['# JOB.fwhm_mm must be either a scalar or [1x3] vector!']);
end
matlabbatchs{1}.spm.spatial.smooth.fwhm = JOB.fwhm_mm;
prefix=['s',num2str(f)];
matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;

need2run=0;
for j=1:size(fnames,2)
 for n=1:size(fnames,1)
  [a,b,c]=fileparts(fnames{n,j});
  fnames{n,j} = [a '/' prefix b c];
  need2run=~exist(fnames{n,j},'file')||need2run;
 end
end
if need2run || JOB.overwrite
 spm_jobman('initcfg')
 spm_jobman('run', matlabbatchs)
end
JOB.fnames = fnames;
end
