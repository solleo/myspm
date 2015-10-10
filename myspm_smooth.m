function EXP = myspm_smooth(EXP)
% EXP = myspm_smooth(EXP)
%
% EXP requires:
%  .fnames
%  .fhwm
%
% (cc) 2015, sgKIM. 

if ~nargin, help(mfilename); return; end

% if numel(EXP.fwhm) == 1
%   fwhm = [EXP.fwhm EXP.fwhm EXP.fwhm];
% else
%   fwhm = EXP.fwhm;
% end
% fnames={};
% if iscell(EXP.fnames)
%   j=1;
%   for c=1:numel(EXP.fnames)
%     hdr = load_untouch_header_only(EXP.fnames{c});
%     for t=1:hdr.dime.dim(5)
%       fnames{j} = [EXP.fnames{c},',',num2str(t)];
%       j=j+1;
%     end
%   end
% else % only string
%   j=1;
%   hdr = load_untouch_header_only(EXP.fnames);
%   for t=1:hdr.dime.dim(5)
%     fnames{j} = [EXP.fnames,',',num2str(t)];
%     j=j+1;
%   end
% end
fnames=EXP.fnames;
if ~iscell(fnames), fnames={fnames}; end

matlabbatchs={};
matlabbatchs{1}.spm.spatial.smooth.data = fnames;
matlabbatchs{1}.spm.spatial.smooth.fwhm = repmat(EXP.fwhm,[1,3]);
matlabbatchs{1}.spm.spatial.smooth.dtype = 0;
matlabbatchs{1}.spm.spatial.smooth.im = 0;
prefix=['s',sprintf('%0.1f',EXP.fwhm)];
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