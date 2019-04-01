function EXP = myspm_smooth(EXP)
% EXP = myspm_smooth(EXP)
%
% EXP requires:
%  .fname    '1xN' for a file name
%  .fwhm_mm  [1x1] or [1x3]
%
% (cc) 2015, sgKIM.

if ~nargin, help(mfilename); return; end
if ~isfield(EXP,'overwrite'), EXP.overwrite=0; end
%% for compatibility with the old code:
if isfield(EXP,'fnames') && iscell(EXP.fnames)
 for y=1:numel(EXP.fnames)
  EXP2=[];
  EXP2.fwhm_mm=EXP.fhwm_mm;
  EXP2.fname=EXP.fnames{y};
  myspm_smooth(EXP2);
 end
 return
end
%%
ls(EXP.fname);
matlabbatchs={};
matlabbatchs{1}.spm.spatial.smooth.dtype = 0;  % same precision
matlabbatchs{1}.spm.spatial.smooth.im = 0;     % implicit masking
%% Find # of volumes and set them as inputs
hdr=spm_vol_nifti(EXP.fname);
if numel(hdr.private.dat.dim) > 3
 T=hdr.private.dat.dim(4);
else
 T=1;
end
matlabbatchs{1}.spm.spatial.smooth.data = cell(T,1);
for t=1:T
 matlabbatchs{1}.spm.spatial.smooth.data{t,1} = [EXP.fname,',',num2str(t)];
end
%% Set a prefix for the kernel size
if numel(EXP.fwhm_mm)==3 % input is [1x3] or [3x1]
 if (EXP.fwhm_mm(1) == EXP.fwhm_mm(2)) ...
   && (EXP.fwhm_mm(2) == EXP.fwhm_mm(3)) % isotropic?
  f=num2str(EXP.fwhm_mm(3));
 else % anisotropic
  f=[num2str(EXP.fwhm_mm(1)),'x',num2str(EXP.fwhm_mm(2)), ...
   'x',num2str(EXP.fwhm_mm(3))];
 end
elseif numel(EXP.fwhm_mm)==1 % input is [1x1]
 f=num2str(EXP.fwhm_mm);
 EXP.fwhm_mm = [1 1 1].*EXP.fwhm_mm;
else
 error(['# EXP.fwhm_mm must be either a scalar or [1x3] vector!']);
end
matlabbatchs{1}.spm.spatial.smooth.fwhm = EXP.fwhm_mm;
prefix=['s',f];
matlabbatchs{1}.spm.spatial.smooth.prefix = prefix;
%% run if needed
[p1,f1,e1]=myfileparts(EXP.fname);
fname_out = fullfile(p1,['/' prefix f1 e1]);
if ~exist(fname_out,'file') || EXP.overwrite
 spm_jobman('initcfg')
 spm_jobman('run', matlabbatchs)
end
end
