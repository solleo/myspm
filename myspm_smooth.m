function EXP = myspm_smooth(EXP)
% EXP = myspm_smooth(EXP)
%
% EXP requires:
%  .fnames
%  .fhwm
%
% (cc) 2015, sgKIM. 

if ~nargin, help(mfilename); return; end
fnames=EXP.fnames;
if ~iscell(fnames), fnames={fnames}; end

matlabbatchs={};
matlabbatchs{1}.spm.spatial.smooth.data = fnames;
matlabbatchs{1}.spm.spatial.smooth.fwhm = repmat(EXP.fwhm,[1,3]);
matlabbatchs{1}.spm.spatial.smooth.dtype = 0;
matlabbatchs{1}.spm.spatial.smooth.im = 0;
if double(int64(EXP.fwhm)) == EXP.fwhm
  prefix=['s',sprintf('%d',EXP.fwhm)];
else
  prefix=['s',sprintf('%0.1f',EXP.fwhm)];
end
matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;

need2run=0;
for j=1:size(fnames,2)
  for n=1:size(fnames,1)
    [a,b,c]=fileparts(fnames{n,j});
    fnames{n,j} = [a '/' prefix b c];
    need2run=~exist(fnames{n,j}(1:end-2),'file')||need2run;
  end
end
if need2run
  spm_jobman('initcfg')
  spm_jobman('run', matlabbatchs)
end
EXP.fnames = fnames;
end