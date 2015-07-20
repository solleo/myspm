function myspm_vbm8ss (filenames, cfg)
% myspm_vbm8ss (filenames, cfg)
% 
% runs skullstripping using vbm8
%
% input:
% filenames <1xN> is a CELL structure for T1w filename(s)
% cfg.isAsian = 1; when the subjects are asian
% 
% (cc) sgKIM, solleo@gmail.com

ver=spm('version');
if str2double(ver(4)) ~= 8
  error('vbm8 requires spm8');
  %myspm_set_spm8
end
spm('Defaults','pet')

if ~exist('cfg','var')
  cfg=[];
end

estwrite=[];
% data
for i=1:numel(filenames)
estwrite.data{i} = [filenames{i},',1'];
end

% Estimation options
%=======================================================================
estwrite.opts.tpm       = {fullfile(spm('dir'),'toolbox','Seg','TPM.nii')}; % TPM.nii
estwrite.opts.ngaus     = [2 2 2 3 4 2];  % Gaussians per class
if isfield(cfg,'isAsian')
  estwrite.opts.affreg    = 'eastern';    % Affine regularisation
else
  estwrite.opts.affreg    = 'mni';    % Affine regularisation
end
estwrite.opts.warpreg   = 4;      % Warping regularisation
estwrite.opts.biasreg   = 0.0001; % Bias regularisation
estwrite.opts.biasfwhm  = 60;     % Bias FWHM
estwrite.opts.samp      = 3;      % Sampling distance

% Writing options
%=======================================================================

% segmentations:
%   native    0/1   (none/yes)
%   warped    0/1   (none/yes)
%   modulated 0/1/2 (none/affine+nonlinear/nonlinear only)
%   dartel    0/1/2 (none/rigid/affine)

estwrite.output.bias.native  = 1;
estwrite.output.bias.warped  = 0;
estwrite.output.bias.affine  = 0;

estwrite.output.label.native = 0;
estwrite.output.label.warped = 0;
estwrite.output.label.dartel = 0;

estwrite.output.GM.native = 1;  % GM
estwrite.output.GM.warped = 0;  % GM
estwrite.output.GM.mod    = 0;  % GM
estwrite.output.GM.dartel = 0;  % GM

estwrite.output.WM.native = 1;  % WM
estwrite.output.WM.warped = 0;  % WM
estwrite.output.WM.mod    = 0;  % WM
estwrite.output.WM.dartel = 0;  % WM

estwrite.output.CSF.native = 1; % CSF
estwrite.output.CSF.warped = 0; % CSF
estwrite.output.CSF.mod    = 0; % CSF
estwrite.output.CSF.dartel = 0; % CSF

% jacobian determinant 0/1 (none/yes)
estwrite.output.jacobian.warped = 0;

% order is [forward inverse]
estwrite.output.warps = [0 0];

% Extended writing options
%=======================================================================
estwrite.extopts.dartelwarp.normlow = struct([]); % spm default
estwrite.extopts.print       = 0; % Display and print results
estwrite.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
estwrite.extopts.mrf       = 0.15; % MRF weighting
estwrite.extopts.ornlm     = 0.7;  % use ORNLM filter with weighting: 0 - no ORNLM

matlabbatch={};
matlabbatch{1}.spm.tools.vbm8.estwrite=estwrite;
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

for i=1:numel(filenames)
  [path,name,~]=fileparts(filenames{i});
  setenv('FSLOUTPUTTYPE','NIFTI')
  cd(path)
  system(['fslmaths p1',name,' -add p2',name,' -add p3',name,' -bin -mul ',...
    'm',name,' bm',name]);
end

end
