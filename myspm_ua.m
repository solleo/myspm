function EXP = myspm_ua (EXP)
% EXP = myspm_ua (EXP)
% unwarp, realignment using SPM8 "u"
% and slice-timing correction using SPM12 "a"
%
% needs EXP to have the following fields:
%  EXP.fname_epi  : [<#sessx1> cell OR <string>] of full path for the resting-state EPI file(s)
%  EXP.fname_dcm  : <string> full path of a dicom file to read slice timing
%  EXP.fmap       :
%  EXP.fmap.fname_query
% (EXP.dir_exp)   : <string> full path to save results (default: path of
%                   input image)
% (EXP.fname_t1w) : for VDM figure
%
% [TODO]: fieldmaps
% (cc) 2015, sgKIM. solleo@gmail.com

spm_jobman('initcfg');
spm_figure('GetWin','Graphics');

% unwarping
EXP = myspm_u(EXP);

% and slice-timing correction
EXP = myspm_a(EXP);
close all

if isfield(EXP,'dir_fig')
  cd(EXP.dir_exp)
  % and move the spm .ps file to 'fig' directory
  src = [EXP.dir_exp,'/spm_',datestr(date,'yyyymmmdd'),'.ps'];
  if exist(src,'file')
    movefile(src, ...
      [EXP.dir_fig,'/01_realign_unwarp_slicetiming_', ...
      datestr(date,'yyyymmmdd'),'.ps']);
  end
end

end


