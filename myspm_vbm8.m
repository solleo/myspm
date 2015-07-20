% This function performs brain extraction using vbm8
% >> vbm8_brainExt (filenames)
% "filenames" is a CELL structure for T1w filename(s)
% 
% (cc) sgKIM, solleo@gmail.com

function myspm_vbm8 (filenames)
%spm pet
estwrite=[];

% data
for i=1:numel(filenames)
estwrite.data{i} = [filenames{i},',1'];
end

% Estimation options
%=======================================================================
estwrite.opts.tpm       = {fullfile(spm('dir'),'toolbox','Seg','TPM.nii')}; % TPM.nii
estwrite.opts.ngaus     = [2 2 2 3 4 2];  % Gaussians per class
estwrite.opts.affreg    = 'mni';    % Affine regularisation
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

% order is [native normalised modulated dartel]
estwrite.output.GM.native = 1;  % GM
estwrite.output.GM.warped = 0;  % GM
estwrite.output.GM.mod    = 1;  % GM
estwrite.output.GM.dartel = 0;  % GM

estwrite.output.WM.native = 1;  % WM
estwrite.output.WM.warped = 0;  % WM
estwrite.output.WM.mod    = 1;  % WM
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
[path,~,~]=fileparts(which('spm'));
estwrite.extopts.dartelwarp.normhigh.darteltpm  = {[path,'/toolbox/vbm8/Template_1_IXI550_MNI152.nii']}; % dartel normalization: 0 - spm default; 1 - yes
estwrite.extopts.print       = 1; % Display and print results
estwrite.extopts.cleanup     = 1;    % Cleanup: 0 - no; 1 - light; 2 -thorough
estwrite.extopts.mrf       = 0.15; % MRF weighting
estwrite.extopts.ornlm     = 0.7;  % use ORNLM filter with weighting: 0 - no ORNLM

matlabbatch={};
matlabbatch{1}.spm.tools.vbm8.estwrite=estwrite;
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);

end
