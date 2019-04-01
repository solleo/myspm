function EXP = myspm_unwarp(EXP)
% EXP = myspm_unwarp(EXP)
%
% unified unwarping and realignment
%
% EXP requires:
%  .fname_epi 
% (.fname_mag, fname_pha) with json files!
%  .fname_t1w  for VDM preparation
%
% (cc) 2015, 2019, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

[p1,f1,e1] = myfileparts(EXP.fname_epi);

%% preparing VDM
if isfield(EXP,'fname_mag') && isfield(EXP,'fname_pha')
  if ~isfield(EXP,'TEs_fmap')
    [p5,f5,~] = myfileparts(EXP.fname_mag);
    try
      json1 = readjson([p5,'/',f5,'.json']);
      [p5,f5,~] = myfileparts(EXP.fname_pha);
      json2 = readjson([p5,'/',f5,'.json']);
      EXP.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      %         EXP.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  EXP.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  [~,f_ph,e_ph] = myfileparts(EXP.fname_pha);
  EXP.fname_vdm = [p1,'/vdm5_sc',f_ph,e_ph];
  myspm_prepare_vdm(EXP.fname_mag, EXP.fname_pha, EXP.TEs_fmap, fname_epi, ...
    EXP.totalreadout_msec, EXP.fname_t1w);
  ls(EXP.fname_vdm)
end
%% unwarp+realign to MEAN IMAGE
% outputs:
% [1]  rp_${epi}.txt   : six rigid-body motion parameters [mm & rad]
% [2]  ${epi}.mat      : [4x4xT] realign transform
% [3]  ${epi}_uw.mat   : unwarping meta data
% [4]  u${epi}.nii     : unwarped/realigned image
% [5]  meanu${epi}.nii : mean image of [4]
hdr = load_untouch_header_only([p1,'/',f1,e1]);
NumFrames = hdr.dime.dim(5);
realignunwarp1=[];
for t=1:NumFrames
  realignunwarp1.data.scans{t,1} = [p1,'/',f1,e1,',',num2str(t)];
end
if isfield(EXP,'fname_vdm')
  ls(EXP.fname_vdm)
  realignunwarp1.data.pmscan = {[EXP.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(EXP,'realign_sep')
  realignunwarp1.eoptions.sep = EXP.realign_sep;
end
realignunwarp1.eoptions.fwhm = 5;
realignunwarp1.eoptions.rtm = 1; % because MEAN image is used in coregistration. (although RP is still w.r.t. the 1st image.. could be confusing?)
realignunwarp1.eoptions.einterp = 2;
realignunwarp1.eoptions.ewrap = [0 0 0];
realignunwarp1.eoptions.weight = '';
realignunwarp1.uweoptions.basfcn = [12 12];
realignunwarp1.uweoptions.regorder = 1;
realignunwarp1.uweoptions.lambda = 100000;
realignunwarp1.uweoptions.jm = 0;
realignunwarp1.uweoptions.fot = [4 5];
realignunwarp1.uweoptions.sot = [];
realignunwarp1.uweoptions.uwfwhm = 4;
realignunwarp1.uweoptions.rem = 1;
realignunwarp1.uweoptions.noi = 5;
realignunwarp1.uweoptions.expround = 'Average';
realignunwarp1.uwroptions.uwwhich = [2 1];
realignunwarp1.uwroptions.rinterp = 4;
realignunwarp1.uwroptions.wrap = [0 0 0];
realignunwarp1.uwroptions.mask = 1;
realignunwarp1.uwroptions.prefix = 'u';
matlabbatch={};
matlabbatch{1}.spm.spatial.realignunwarp = realignunwarp1;
fname_in=[p1,'/',f1,e1];
ls(fname_in);
fname_output = [p1,'/u',f1,e1];
if ~exist(fname_output,'file') || EXP.overwrite
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  ls(fname_output)
  
  % visualize unwarping results:
  compare_unwarped([p1,'/',f1,e1], [p1,'/u',f1,e1], EXP.fname_t1w)
end
end



%-------------------------------------------------------------------------------
% spm('Defaults','fmri')
% if ~isfield(EXP,'overwrite'), EXP.overwrite =0; end
% pwd0=pwd;
% 
% disp('# Realigning and unwarping..');
% if iscell(EXP.fname_epi)
%  n_sess = numel(EXP.fname_epi);
% else
%  n_sess = 1;
%  EXP.fname_epi = {EXP.fname_epi};
% end
% EXP.n_sess = n_sess;
% if ~isfield(EXP,'dir_exp')
%  if strcmp(EXP.fname_epi{1}(1),'/')
%   [dir_exp,~,~] = fileparts(EXP.fname_epi{1});
%  else
%   dir_exp=pwd;
%  end
% else
%  dir_exp = EXP.dir_exp;
% end
% if ~isfield(EXP,'fname_t1w'),   EXP.fname_t1w='';  end
% cd (dir_exp);
% 
% for sess = 1:n_sess
%  disp(['# Session: ',num2str(sess),'/',num2str(n_sess)]);
%  %% 1. Realignment (coregister only) to the MEAN IMAGE
%  P   = EXP.fname_epi{sess};
%  hdr = load_nii_hdr(P);
%  FlagsC=[];
%  FlagsC.quality  = 1; % 1=max # of voxels used to compute rigid motion
%  FlagsC.fwhm     = 2*mean(hdr.dime.pixdim(2:4)); % fwhm = 2 x mean voxel size
%  FlagsC.sep      = 4; % down-sampling factor
%  
%  if isfield(EXP,'firstimage')&&EXP.firstimage
%   FlagsC.rtm      = 0;
%   disp(['# Realigning to the first image: ',P,'...']);
%  else
%   FlagsC.rtm      = 1; % register to mean (mean is computed after 1st alignment)
%   disp(['# Realigning to mean image: ',P,'...']);
%  end
%  FlagsC.interp   = 4; % 4th degree B-spline
%  FlagsC.PW       = ''; % weighting image (1/std)
%  FlagsC.graphics = 0;
%  FlagsC.lkp      = 1:6;
%  [dir1,name1,ext1] = fileparts(P);
%  if ~exist([dir_exp,'/rp_',name1,'.txt'],'file')
%   spm_realign(P,FlagsC);
%  end
%  
%  %% 2. Unwarp and resample
%  % 2-1. find field map files:
%  
%  if isfield(EXP,'fmap')&&~isempty(EXP.fmap)	% fieldmap scan available
%   if (~isfield(EXP.fmap,'shortmag') || ~isfield(EXP.fmap,'phasedif')) && ~isfield(EXP.fmap,'vdm')
%    [~,fullnames] = mydir(EXP.fmap.fname_query,1);
%    EXP.fmap.shortmag = fullnames{1};
%    EXP.fmap.phasedif = fullnames{3};
%    
%    % 2-1-2. find echo time from the filenames..
%    TEs=[0 0];
%    [~,fname2,~] = fileparts(EXP.fmap.shortmag);
%    idx1 = strfind(fname2,'_Te');
%    TEs(1) = str2double(fname2(idx1+3:end));
%    [~,fname2,~] = fileparts(EXP.fmap.phasedif);
%    idx1 = strfind(fname2,'_Te');
%    TEs(2) = str2double(fname2(idx1+3:end));
%    disp(['> TEs = ',num2str(TEs(1)),' / ',num2str(TEs(2))]);
%   end
%   % 2-1-3. check if VDM is already computed
%   if isfield(EXP.fmap,'vdm')
%    fmap = EXP.fmap.vdm;
%   else
%    [path1,a,b] = fileparts(EXP.fmap.phasedif);
%    fmap = fullfile(path1, ['vdm5_sc',a,b]);
%   end
%   
%   if exist(fmap,'file') && ~EXP.overwrite
%    disp('# Using pre-computed voxel displacement map (VDM):');
%    ls(fmap)
%   elseif isfield(EXP.fmap,'shortmag');
%    disp('# Computing voxel displacement map (VDM) from:');
%    fprintf('short-magnitude: ' ); ls(EXP.fmap.shortmag);
%    fprintf('phase-difference: '); ls(EXP.fmap.phasedif);
%    TRT=prepare_vdm ( EXP.fmap.shortmag, EXP.fmap.phasedif, TEs, ...
%     P, EXP.epitype, [], EXP.fname_t1w );
%    [~,a,b] = fileparts(P);
%    trg=[dir_exp,'/vdm5_sc',a,b];
%    if ~exist(trg,'file')
%     system(['ln -s ',fmap,' ',trg]);
%    end
%    EXP.fmap.fname_origvdm = fmap;
%    EXP.fmap.fname_vdm = trg;
%    fmap = trg;
%    
%    fid = fopen([dir_exp,'/fmap_params_for_',name1,ext1,'.txt'],'w');
%    if isfield(EXP,'epitype')
%     fprintf(fid, 'EPI type: %s\n', EXP.epitype);
%    end
%    fprintf(fid, 'Magnitude with a short TE: %s\n',EXP.fmap.shortmag);
%    fprintf(fid, 'Phase difference (deg): %s\n',EXP.fmap.phasedif);
%    fprintf(fid, 'VDM original filename: %s\n', EXP.fmap.fname_origvdm);
%    fprintf(fid, 'TE(msec)= %f/%f\n', TEs(1), TEs(2));
%    if ~exist('TRT','var')
%     TRT=nan;
%    end
%    fprintf(fid, 'TotalReadoutTime(msec)= %f\n', TRT);
%    fprintf(fid, '=====%s====\n\n',date);
%    fclose(fid);
%   end
%  else
%   disp(['# No fieldmap given, will estimate from the EPI file, ', ...
%    'which makes only limited correction']);
%   fmap = [];
%   EXP.fmap = [];
%  end
%  
%  % 2-2. Estimate unwarping parameters
%  fname_mat = fullfile(dir_exp, ['u',name1,'.mat']);
%  if ~exist(fname_mat,'file') || EXP.overwrite
%   uw_est_flags=[];
%   uw_est_flags.order = [12 12]; % order of basis functions for each dimension
%   uw_est_flags.sfP   = fmap; % static field file (VDM)
%   uw_est_flags.regorder = 1; % regularization on the N-th derivative of field
%   uw_est_flags.jm = 0; % jacobian modulation
%   uw_est_flags.fot = [4 5]; % first order of terms
%   uw_est_flags.sot = [];
%   uw_est_flags.fwhm = 4;
%   uw_est_flags.rem = 1; % re-estimation of RP
%   uw_est_flags.exp_round = 'Average';
%   uw_est_flags.noi = 5; % # of iterations
%   uw_est_flags.hold = [1 1 1 0 1 0];
%   ds = spm_uw_estimate(P,uw_est_flags);
%   save(fname_mat,'ds');
%  else
%   load(fname_mat,'ds');
%  end
%  
%  % 2-3. Resample unwarped images
%  fname_uwarped = fullfile(dir_exp,['u',name1,'.nii']);
%  if ~exist(fname_uwarped ,'file') || EXP.overwrite
%   uw_write_flags=[];
%   uw_write_flags.mask = 1;
%   uw_write_flags.mean = 1;
%   uw_write_flags.rinterp = 4;
%   uw_write_flags.wrap = [0 1 0];
%   uw_write_flags.which = [2 1];
%   uw_write_flags.udc = 1;
%   uw_write_flags.prefix ='u';
%   spm_uw_apply(ds, uw_write_flags);
%   disp(['> Unwarped and resampled: ']);
%   ls([dir_exp,'/u',name1,ext1]);
%  end
%  [p1,n1,e1] = fileparts(EXP.fname_epi{sess});
% %  if isempty(p1)
% %   EXP.fname_epi{sess} = [p1,'a',n1,e1];
% %  else
% %   EXP.fname_epi{sess} = [p1,'/a',n1,e1];
% %  end
% end
% 
% cd(pwd0);

