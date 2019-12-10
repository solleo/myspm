function JOB = myspm_unwarp(JOB)
% JOB = myspm_unwarp(JOB)
%
% unified unwarping and realignment
%
% JOB requires:
%  .fname_epi 
% (.fname_mag, fname_pha) with json files!
%  .fname_t1w  for VDM preparation
%
% (cc) 2015, 2019, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

[p1,f1,e1] = myfileparts(JOB.fname_epi);

%% preparing VDM
if isfield(JOB,'fname_mag') && isfield(JOB,'fname_pha')
  if ~isfield(JOB,'TEs_fmap')
    [p5,f5,~] = myfileparts(JOB.fname_mag);
    try
      json1 = readjson([p5,'/',f5,'.json']);
      [p5,f5,~] = myfileparts(JOB.fname_pha);
      json2 = readjson([p5,'/',f5,'.json']);
      JOB.TEs_fmap = [json1.EchoTime json2.EchoTime]*1000;
      %         JOB.totalreadout_msec = 1/json2.PixelBandwidth*1000;
    catch
      error(['Cannot find json files for fieldmap. ',...
        'Enter .TEs_fmap and .totalreadout_msec manually']);
    end
  end
  JOB.totalreadout_msec = json.TotalReadoutTime * 1000; % THIS IS of the EPI to unwarp!
  [~,f_ph,e_ph] = myfileparts(JOB.fname_pha);
  JOB.fname_vdm = [p1,'/vdm5_sc',f_ph,e_ph];
  myspm_prepare_vdm(JOB.fname_mag, JOB.fname_pha, JOB.TEs_fmap, fname_epi, ...
    JOB.totalreadout_msec, JOB.fname_t1w);
  ls(JOB.fname_vdm)
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
if isfield(JOB,'fname_vdm')
  ls(JOB.fname_vdm)
  realignunwarp1.data.pmscan = {[JOB.fname_vdm,',1']};
else
  realignunwarp1.data.pmscan = {''};
end
realignunwarp1.eoptions.quality = 1;
realignunwarp1.eoptions.sep = 4;
if isfield(JOB,'realign_sep')
  realignunwarp1.eoptions.sep = JOB.realign_sep;
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
if ~exist(fname_output,'file') || JOB.overwrite
  fname_matlabbatch=[p1,'/myspm_fmriprep4_',f1,'.mat'];
  save(fname_matlabbatch,'matlabbatch');
  spm_jobman('run', matlabbatch);
  ls(fname_output)
  
  % visualize unwarping results:
  compare_unwarped([p1,'/',f1,e1], [p1,'/u',f1,e1], JOB.fname_t1w)
end
end



%-------------------------------------------------------------------------------
% spm('Defaults','fmri')
% if ~isfield(JOB,'overwrite'), JOB.overwrite =0; end
% pwd0=pwd;
% 
% disp('# Realigning and unwarping..');
% if iscell(JOB.fname_epi)
%  n_sess = numel(JOB.fname_epi);
% else
%  n_sess = 1;
%  JOB.fname_epi = {JOB.fname_epi};
% end
% JOB.n_sess = n_sess;
% if ~isfield(JOB,'dir_exp')
%  if strcmp(JOB.fname_epi{1}(1),'/')
%   [dir_exp,~,~] = fileparts(JOB.fname_epi{1});
%  else
%   dir_exp=pwd;
%  end
% else
%  dir_exp = JOB.dir_exp;
% end
% if ~isfield(JOB,'fname_t1w'),   JOB.fname_t1w='';  end
% cd (dir_exp);
% 
% for sess = 1:n_sess
%  disp(['# Session: ',num2str(sess),'/',num2str(n_sess)]);
%  %% 1. Realignment (coregister only) to the MEAN IMAGE
%  P   = JOB.fname_epi{sess};
%  hdr = load_nii_hdr(P);
%  FlagsC=[];
%  FlagsC.quality  = 1; % 1=max # of voxels used to compute rigid motion
%  FlagsC.fwhm     = 2*mean(hdr.dime.pixdim(2:4)); % fwhm = 2 x mean voxel size
%  FlagsC.sep      = 4; % down-sampling factor
%  
%  if isfield(JOB,'firstimage')&&JOB.firstimage
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
%  if isfield(JOB,'fmap')&&~isempty(JOB.fmap)	% fieldmap scan available
%   if (~isfield(JOB.fmap,'shortmag') || ~isfield(JOB.fmap,'phasedif')) && ~isfield(JOB.fmap,'vdm')
%    [~,fullnames] = mydir(JOB.fmap.fname_query,1);
%    JOB.fmap.shortmag = fullnames{1};
%    JOB.fmap.phasedif = fullnames{3};
%    
%    % 2-1-2. find echo time from the filenames..
%    TEs=[0 0];
%    [~,fname2,~] = fileparts(JOB.fmap.shortmag);
%    idx1 = strfind(fname2,'_Te');
%    TEs(1) = str2double(fname2(idx1+3:end));
%    [~,fname2,~] = fileparts(JOB.fmap.phasedif);
%    idx1 = strfind(fname2,'_Te');
%    TEs(2) = str2double(fname2(idx1+3:end));
%    disp(['> TEs = ',num2str(TEs(1)),' / ',num2str(TEs(2))]);
%   end
%   % 2-1-3. check if VDM is already computed
%   if isfield(JOB.fmap,'vdm')
%    fmap = JOB.fmap.vdm;
%   else
%    [path1,a,b] = fileparts(JOB.fmap.phasedif);
%    fmap = fullfile(path1, ['vdm5_sc',a,b]);
%   end
%   
%   if exist(fmap,'file') && ~JOB.overwrite
%    disp('# Using pre-computed voxel displacement map (VDM):');
%    ls(fmap)
%   elseif isfield(JOB.fmap,'shortmag');
%    disp('# Computing voxel displacement map (VDM) from:');
%    fprintf('short-magnitude: ' ); ls(JOB.fmap.shortmag);
%    fprintf('phase-difference: '); ls(JOB.fmap.phasedif);
%    TRT=prepare_vdm ( JOB.fmap.shortmag, JOB.fmap.phasedif, TEs, ...
%     P, JOB.epitype, [], JOB.fname_t1w );
%    [~,a,b] = fileparts(P);
%    trg=[dir_exp,'/vdm5_sc',a,b];
%    if ~exist(trg,'file')
%     system(['ln -s ',fmap,' ',trg]);
%    end
%    JOB.fmap.fname_origvdm = fmap;
%    JOB.fmap.fname_vdm = trg;
%    fmap = trg;
%    
%    fid = fopen([dir_exp,'/fmap_params_for_',name1,ext1,'.txt'],'w');
%    if isfield(JOB,'epitype')
%     fprintf(fid, 'EPI type: %s\n', JOB.epitype);
%    end
%    fprintf(fid, 'Magnitude with a short TE: %s\n',JOB.fmap.shortmag);
%    fprintf(fid, 'Phase difference (deg): %s\n',JOB.fmap.phasedif);
%    fprintf(fid, 'VDM original filename: %s\n', JOB.fmap.fname_origvdm);
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
%   JOB.fmap = [];
%  end
%  
%  % 2-2. Estimate unwarping parameters
%  fname_mat = fullfile(dir_exp, ['u',name1,'.mat']);
%  if ~exist(fname_mat,'file') || JOB.overwrite
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
%  if ~exist(fname_uwarped ,'file') || JOB.overwrite
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
%  [p1,n1,e1] = fileparts(JOB.fname_epi{sess});
% %  if isempty(p1)
% %   JOB.fname_epi{sess} = [p1,'a',n1,e1];
% %  else
% %   JOB.fname_epi{sess} = [p1,'/a',n1,e1];
% %  end
% end
% 
% cd(pwd0);

