function EXP=myspm_fmriglm (EXP)
% EXP=myspm_fmriglm (EXP)
%
% This script helps you to run GLMs for the 1st-level fMRI analysis using 
%   fMRI model specification & specification. A result report of the 
%   1st-level GLM is created using myspm_result.m and myspm_graph.m. The 
%   important difference from the simple GLM is that it deals with temporal
%   aspects of fMRI data including slow hedonamics, pre-whitening for auto-
%   correlation, as well as high-pass filtering for detrending.
%
% [INPUTS]
% EXP requires fields for myspm_fmriglm.m:
% -input files
%  .files_query '1xN' for query to find image filenames using "dir" function
%                (e.g., "/path/to/subj/EPI_r*.nii")
% or
%  .filenames   {1xJ} EPI files for J sessions (instead of files_query)
%
% -output directory
%  .dir_glm     '1xN' directory to save SPM results
% or
%  .dir_base    '1xN' directory for a subdirectory that has SPM.mat
% (.dir_prefix)
%
% -HRF model: default is a canonical HRF
% (.microt0)        [1x1] where you STC [0,1] default = 0.5
% (.fir)            (1x1) finite-impulse-response model parameters
% (.fir.length_sec) [1x1] Post-stimulus window length (in seconds)
%                         (default: stimulus-duration)
% (.fir.order)      [1x1] the number of basis functions (or # of shifts/lags)
%                         (default: stimulus-duration/TR_sec)
%
% (.num_ISSS)       [1x1] interleaved silent steady state (Schwarzbauer, 2006)
%                         Number of images during the "scanning block"
%
% -for fMRI design matrix
%  .TR_sec
% (.stimdur_sec)       This could be a vector
% (.cond(j,k).name)    for j-th session, k-th condition
% (.cond(j,k).onset)
% (.cond(j,k).dur)
% (.cond(j,k).param)
%
% -autocorrelation modeling:
% (.hpfcutoff)     [1x1] high-pass cut-off period default=128 [s]
% (.autocorr)      '1xN'  "none" | "AR(1)" (default) | "FAST"
%
% -nuissance regressors
% (.suffix_rp)     {1xJ} rigid-body motion parameters for J sessions
%                        "rp_{suffix_rp{j}}.txt"
% (.num_rp)        [1x1] # of columns to include
% (.prefix_cc)     {1xJ} CompCor parameters for J sessions
%                        "{prefix_cc{j}}_eigenvec.txt"
% (.num_cc)        [1x1] # of columns to include
% (.scrb_thres)    [1x1] threshold to scrub on dShift/dt [mm]
%
% -custom regressor
% (.reg(j,k).name) '1xN' (additional regressors with/without FIR)
% (.reg(j,k).val)  [Tx1] for the k-th regressor in the j-th session
% (.reg(j,k).fir.order)   [1x1] number of equal-spaced shifts within a kernel
% (.reg(j,k).fir.length)  [1x1] length of a kernel in samples (TRs)
%
% -PPI,
% (.fname_ppi)
%
% -for myspm_result.m:
% (.noreport)       [1x1] creates no reports (stops after cntrst)
% (.masking)        '1xN' filename for an explicit (inclusive) mask
%                         * NECESSARY for images with negative values!
% (.thres.desc)     '1xN'  'FWE','none', or 'cluster'(default)
% (.thres.alpha)    [1x1]  alpha level (default=0.05)
% (.thres.extent)   [1x1]  extent threshold of clusters in voxels (default=0)
% (.thres.clusterInitAlpha)   [1x1] cluster forming height threshold (default=0.0001)
% (.thres.clusterInitExtent)  [1x1] cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct)   '1xN' fullpath filename for background anatomical image for orthogonal slices
%                         (default='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)       {1xK} Title texts for K contrasts for SPM result report (default={'positive','negative'})
% (.dir_sum)        '1xN' a summary directory into where you want to copy significant results
% (.append)         [1x1] whether to append results into an existing report (default=0)
% (.print)          [1x1] whether to generate results (default=1)
% (.mygraph.x_name) '1xN' x-axis label in a scatterplot (default='x')
% (.mygraph.y_name) '1xN' y-axis label in a scatterplot (default='y')
% (.atlas)          '1xN' atlas to find anatomical names: 'fsl' (default) or 'spm12'
%
% -for myspm_mycon.m:
% (.conname)       {1xK} names of K contrasts
% (.convec)        [KxP] contrast vectors for K contrasts, when the number
%                        of betas' is larger than P, zeros will be padded.
%
% -misc.
% (.keepResidual)  [1x1] to keep residual timeseries in RESI.nii
%
%
% [OUTPUT] 
%
% (cc) 2015-2019. sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com/


% Example:
%
% Results (from spm_spm.m):
% The following SPM.fields are set by spm_spm (unless specified)
%
%     xVi.V      - estimated non-sphericity trace(V) = rank(V)
%     xVi.h      - hyperparameters  xVi.V = xVi.h(1)*xVi.Vi{1} + ...
%     xVi.Cy     - spatially whitened <Y*Y'> (used by ReML to estimate h)
%
%                            ----------------
%
%     Vbeta     - struct array of beta image handles (relative)
%     VResMS    - file struct of ResMS image handle  (relative)
%     VM        - file struct of Mask  image handle  (relative)
%
%                            ----------------
%
%     xX.W      - if not specified W*W' = inv(x.Vi.V)
%     xX.V      - V matrix (K*W*Vi*W'*K') = correlations after K*W is applied
%     xX.xKXs   - space structure for K*W*X, the 'filtered and whitened'
%                 design matrix
%               - given as spm_sp('Set',xX.K*xX.W*xX.X) - see spm_sp
%     xX.pKX    - pseudoinverse of K*W*X, computed by spm_sp
%     xX.Bcov   - xX.pKX*xX.V*xX.pKX - variance-covariance matrix of
%                 parameter estimates
%                 (when multiplied by the voxel-specific hyperparameter ResMS
%                 of the parameter estimates (ResSS/xX.trRV = ResMS) )
%     xX.trRV   - trace of R*V
%     xX.trRVRV - trace of RVRV
%     xX.erdf   - effective residual degrees of freedom (trRV^2/trRVRV)
%     xX.nKX    - design matrix (xX.xKXs.X) scaled for display
%                 (see spm_DesMtx('sca',... for details)
%
%                            ----------------
%
%     xVol.M    - 4x4 voxel->mm transformation matrix
%     xVol.iM   - 4x4 mm->voxel transformation matrix
%     xVol.DIM  - image dimensions - column vector (in voxels)
%     xVol.XYZ  - 3 x S vector of in-mask voxel coordinates
%     xVol.S    - Lebesgue measure or volume       (in voxels)
%     xVol.R    - vector of resel counts           (in resels)
%     xVol.FWHM - Smoothness of components - FWHM, (in voxels)
%
%                            ----------------
%
%     xCon      - Contrast structure (created by spm_FcUtil.m)
%     xCon.name - Name of contrast
%     xCon.STAT - 'F', 'T' or 'P' - for F/T-contrast ('P' for PPMs)
%     xCon.c    - (F) Contrast weights
%     xCon.X0   - Reduced design matrix (spans design space under Ho)
%                 It is in the form of a matrix (spm99b) or the
%                 coordinates of this matrix in the orthogonal basis
%                 of xX.X defined in spm_sp.
%     xCon.iX0  - Indicates how contrast was specified:
%                 If by columns for reduced design matrix then iX0 contains
%                 the column indices. Otherwise, it's a string containing
%                 the spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
%                 (Usually this is the input argument F_iX0.)
%     xCon.X1o  - Remaining design space (orthogonal to X0).
%                 It is in the form of the coordinates of this matrix in
%                 the orthogonal basis of xX.X defined in spm_sp.
%     xCon.eidf - Effective interest degrees of freedom (numerator df)
%     xCon.Vcon - ...for handle of contrast/ESS image (empty at this stage)
%     xCon.Vspm - ...for handle of SPM image (empty at this stage)
%__________________________________________________________________________
%
% The following images are written to disk:
%
% mask.<ext>                                          - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
%
%                            ----------------
%
% beta_????.<ext>                                     - parameter images
% These are 32-bit (float32) images of the parameter estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.<ext>) are given
% value NaN.
%
%                            ----------------
%
% ResMS.<ext>                           - estimated residual variance image
% This is a 64-bit (float64) image of the residual variance estimate.
% Voxels outside the analysis mask are given value NaN.
%
%                            ----------------
%
% RPV.<ext>                              - estimated resels per voxel image
% This is a 64-bit (float64) image of the RESELs per voxel estimate.
% Voxels outside the analysis mask are given value 0.  These images
% reflect the nonstationary aspects the spatial autocorrelations.
%
%                            ----------------
%
% ResI_????.<ext>                - standardised residual (temporary) images
% These are 64-bit (float64) images of standardised residuals. At most
% maxres images will be saved and used by spm_est_smoothness, after which
% they will be deleted.
%__________________________________________________________________________
%
% References:
%
% Statistical Parametric Maps in Functional Imaging: A General Linear
% Approach. Friston KJ, Holmes AP, Worsley KJ, Poline JB, Frith CD,
% Frackowiak RSJ. (1995) Human Brain Mapping 2:189-210.
%
% Analysis of fMRI Time-Series Revisited - Again. Worsley KJ, Friston KJ.
% (1995) NeuroImage 2:173-181.


%% -2. prompt help
if nargin == 0, help(mfilename); return; end

%% 0. parsing & check inputs
if ~isfield(EXP,'overwrite'), EXP.overwrite=0; end
spm('Defaults','fmri');
spm_jobman('initcfg');

%% 0-1. find glm directory name and create one
if ~isfield(EXP,'dir_glm')
  if ~isfield(EXP,'dir_base')
    EXP.dir_base = pwd;
  end
  if ~isfield(EXP,'model_desc') && isfield(EXP,'model')
    EXP.model_desc = fsss_model_desc(EXP.model);
  end
  if ~isfield(EXP,'dir_prefix'), EXP.dir_prefix=''; end
  EXP.model_desc = [EXP.dir_prefix, EXP.model_desc]; % append prefix
  EXP.dir_glm = [EXP.dir_base, '/glm_', EXP.model_desc];
else
  [~,EXP.model_desc,~] = myfileparts(EXP.dir_glm);
end
[~,~] = mkdir(EXP.dir_glm);

%% 0-2. check files exist
disp('[0b] checking inputs..')
if isfield(EXP,'files_query')
  fnames={}; % this MUST be a column vector <Sx1>
  files = dir(EXP.files_query);  % this will sort the filename!
  [mypath,~,~] = myfileparts(EXP.files_query);
  for n=1:numel(files)
    fnames{n,1} = [mypath,'/',files(n).name,',1'];
  end
elseif isfield(EXP,'filenames')
  for n=1:numel(EXP.filenames)
    fnames{n,1} = [EXP.filenames{n},',1'];
  end
else
  error('You need to specify inputs in EXP.files_query or EXP.filenames');
end
ls(fnames{n,1}(1:end-2));

%% 0-3. set other params
EXP.NumSess = numel(fnames);
for j=1:EXP.NumSess
  hdr   = load_untouch_header_only(fnames{j}(1:end-2));
  NF(j) = hdr.dime.dim(5);
  disp(['[0c] # of frames: run(',num2str(j),')= ',num2str(NF(j))]);
end
EXP.NumFrames=NF;
% find conditions to contruct design matrix
if isfield(EXP,'cond')
  EXP.NumCond = size(EXP.cond,2);
else
  EXP.NumCond = 0; % no condition is given
end
if ~isfield(EXP,'cntrstMtx')
  EXP.NOCNTRST=1;
end

% % this could be a bit smarter...?
% if ~isfield(EXP,'suffix_rp')
%   EXP.suffix_rp=cell(1,EXP.NumSess);
%   for j=1:EXP.NumSess
%     [~,f1,~]=myfileparts(EXP.filenames{j});
%     % TODO: findout corresponding "rp_originalFilename.txt"
%     % from a given input "{unknownPrefix}originalFilename"
%     EXP.suffix_rp{j}=f1(idx_dot+1:end);
%   end
% end

%% 1. model specification
matlabbatch = {};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {''}; % for spm12
matlabbatch{1}.spm.stats.fmri_spec.dir = {EXP.dir_glm};

% Basis functions: hrf or fir
if isfield(EXP,'fir')
  if ~isfield(EXP.fir,'length_sec') && ~isfield(EXP.fir,'order')
    EXP.fir.length_sec = EXP.cond(1,1).dur;
    EXP.fir.order  = round(EXP.cond(1,1).dur/EXP.TR_sec);
  end
  matlabbatch{1}.spm.stats.fmri_spec.bases.fir.length = EXP.fir.length_sec;
  matlabbatch{1}.spm.stats.fmri_spec.bases.fir.order  = EXP.fir.order;
  matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR_sec;
  matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 1;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
else
  matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
  matlabbatch{1}.spm.stats.fmri_spec.timing.RT = EXP.TR_sec;
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16; % # of bins
  if ~isfield(EXP,'microt0')
    EXP.microt0 = 0.5;
  end
  matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = round(EXP.microt0 * 16); % middle point after STC
  matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
end
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';

sess=[];
for j=1:EXP.NumSess
  % find image files
  for t=1:EXP.NumFrames(j)
    sess(j).scans{t,1} = [EXP.filenames{j},',',num2str(t)]; % this must be <Tx1>
  end
  [dir_sess,~,~]=myfileparts(sess(j).scans{t,1});
  sess(j).cond = struct('name',{}, 'onset',{}, 'duration',{}, 'tmod',{});
  if EXP.NumCond
    for k=1:EXP.NumCond
      % onset as a vector
      sess(j).cond(k).onset = EXP.cond(j,k).onset;
      % condition name
      sess(j).cond(k).name  = EXP.cond(j,k).name;
      % condition duration
      sess(j).cond(k).duration = EXP.cond(j,k).dur;
      % For temporal modulation:
      sess(j).cond(k).tmod = 0;
      % For PARAMETRIC modulation (only one & 1st-order parameter)
      if isfield(EXP.cond(j,k),'paramname')
        sess(j).cond(j,k).pmod(1).name  = EXP.cond(j,k).paramname;
        sess(j).cond(k).pmod(1).poly  = 1;
        sess(j).cond(k).pmod(1).param = EXP.cond(j,k).param;
      end
    end % loop for conditions
  end
  
  l=1;
  % custom regressors (with or without FIR)
  if isfield(EXP,'reg')
    disp(['[1a] Covarying custom regressors']);
    NumReg = size(EXP.reg,2);
%     figure;
    for k=1:NumReg
      if ~isfield(EXP,'fir')
        if isfield(EXP.reg(j,k),'name')
          sess(j).regress(l).name = EXP.reg(j,k).name;
        else
          sess(j).regress(l).name = ['reg',num2str(k)];
        end
        sess(j).regress(l).val  = EXP.reg(j,k).val;
        l=l+1;
      else
        % 1. create a FIR kernel
        xBF = struct('name','Finite Impulse Response', ...
          'order',EXP.fir.order, ...
          'length',EXP.fir.length_sec, ...
          'dt',EXP.TR_sec,'UNITS','secs');
        xBF = spm_get_bf(xBF);   
        % 2. convolute a centered regressor & resample it for each order
        x = EXP.reg(j,k).val;
        x = [x; zeros(EXP.NumFrames(j)-size(x,1),1)]; % zero-padding
        xconv = zeros(size(x,1), xBF.order);
        for i = 1:xBF.order
          temp = conv(demean(x), xBF.bf(:,i));
          xconv(:,i) = temp(1:size(x,1));
        end
        % 3. orthogonalize them... or not?
        xconv = spm_orth(xconv);
        % 4. add regressors
        for i = 1:xBF.order
          sess(j).regress(l).name = sprintf('%s*fir(%i)',EXP.reg(j,k).name,i);
          sess(j).regress(l).val = xconv(:,i); 
          l=l+1;
        end
      end
    end
  end
  
  % rigid-motion parameters
  if ~isfield(EXP,'num_rp'), EXP.num_rp=0; end
  if EXP.num_rp
    if ~isfield(EXP,'suffix_rp')
      [~,fname] = myfileparts(EXP.filenames{j});
      [fwhm,isW,isU,isA,origname] = myspm_parse_filename(fname);
      prefixA={'','a'};
      fn_rp = [dir_sess,'/rp_',prefixA{isA+1},origname,'8.txt'];
    else
      fn_rp = [dir_sess,'/rp_',EXP.suffix_rp{j},'.txt'];
    end
    rpname={'dx','dy','dz','rx','ry','rz','dTrans/dt','dRot/dt'};
    disp(['[1b] Covarying ',num2str(EXP.num_rp),' motion parameters: ']);
    disp(rpname(1:EXP.num_rp))
    rp=load(fn_rp);
    rp=[rp(:,1:6) [0; l2norm(diff(rp(:,1:3)))] [0; l2norm(diff(rp(:,4:6)))]];
    for k=1:EXP.num_rp
      sess(j).regress(l).name = rpname{k};
      sess(j).regress(l).val  = rp(:,k);
      l=l+1;
    end
  end
  
  % CompCor parameters
  if ~isfield(EXP,'num_cc'), EXP.num_cc=0; end
  if EXP.num_cc
    disp(['[1c] Covarying ',num2str(EXP.num_cc),' CompCor parameters']);
    if ~isfield(EXP,'prefix_cc')
      [~,fname] = myfileparts(EXP.filenames{j});
      [fwhm,isW,isU,isA,origname] = myspm_parse_filename(fname);
      prefixA={'','a'}; prefixU={'','u'};
      fn_cc = [dir_sess,'/',prefixU{isU+1},prefixA{isA+1},origname,...
        '_n3b0.00-Inf_eigenvec.txt'];
    else
      fn_cc = [dir_sess,'/',EXP.prefix_cc{j},'_eigenvec.txt'];
    end
    cc=load(fn_cc);
    for k=1:EXP.num_cc
      sess(j).regress(l).name = ['cc',num2str(k)];
      sess(j).regress(l).val  = cc(:,k);
      l=l+1;
    end
  end
  
  % Scrubbing parameters
  if isfield(EXP,'scrb_thres')
    idx=find(rp(:,7)>EXP.scrb_thres);
    idx=unique([idx idx-1 idx+1]);
    idx(idx<1)=[];
    idx(idx>EXP.NumFrames(j))=[];
    disp(['[1d] Found ',num2str(numel(idx)),' frames exceeding a threshold of ',...
      num2str(EXP.scrb_thres),' mm in dTrans/dt']);
    if numel(idx)
      for k=1:numel(idx)
        sess(j).regress(l).name=['scrb',num2str(k)];
        sess(j).regress(l).val=zeros(EXP.NumFrames(j),1);
        sess(j).regress(l).val(idx(k))=1;
        l=l+1;
      end
    end
  end
  
  if isfield(EXP,'hpfcutoff')
    sess(j).hpf = EXP.hpfcutoff;
  else
    sess(j).hpf = 128;
  end
  disp(['[1e] high-pass cut-off: ',num2str(sess(j).hpf),' s']);
  
  matlabbatch{1}.spm.stats.fmri_spec.sess = sess;
  matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});  
  if isfield(EXP,'masking')
    disp(['[1f] A mask is given: ',EXP.masking]);
    matlabbatch{1}.spm.stats.fmri_spec.mask = {EXP.masking};
  else
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
  end
  
  % autocorrelation model
  if ~isfield(EXP,'autocorr')
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
  else
    matlabbatch{1}.spm.stats.fmri_spec.cvi = EXP.autocorr;
  end
end

if ~exist([EXP.dir_glm,'/SPM.mat'],'file') || EXP.overwrite
  save ([EXP.dir_glm,'/myspm_fmriglm_1spec.mat'],'matlabbatch');
  spm_jobman('initcfg');
  spm_jobman('run', matlabbatch);
end

%% 2. Model estimation
matlabbatch={};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[EXP.dir_glm,'/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1; % ReML
if isfield(EXP,'keepResidual') % since SPM12?
  matlabbatch{1}.spm.stats.fmri_est.write_residuals = EXP.keepResidual;
end

% running a batch
if (~exist([EXP.dir_glm,'/beta_0001.img'],'file') ...
    && ~exist([EXP.dir_glm,'/beta_0001.nii'],'file')) || EXP.overwrite
  save ([EXP.dir_glm,'/myspm_fmriglm_2est.mat'],'matlabbatch');
  spm_jobman('initcfg');
  spm_jobman('run', matlabbatch);
end

%% (optional: create a NIFTI file to keep ResI.*)
pwd0 = pwd;
cd (EXP.dir_glm)
if isfield(EXP,'keepResidual') && EXP.keepResidual
  % create NIFTI
  setenv('FSLOUTPUTTYPE','NIFTI_GZ')
  myunix('fslmerge -t RES Res_????.nii');
  ls('RES.nii.gz')
  setenv('FSLOUTPUTTYPE','NIFTI')
end
% delete
myunix('rm -f Res_????.nii');
myunix('rm -f ResI_????.nii');
cd (pwd0);

%% 4. Design review & Orthogonality check
% set figure filename
if ~isfield(EXP,'fname_spm_fig')
  today = datestr(now,'yyyymmmdd');
  EXP.fname_spm_fig = fullfile(EXP.dir_glm,[EXP.model_desc,'_spm_',today,'.ps']);
  EXP.fname_spm_fig = strrep(EXP.fname_spm_fig,'>','-gt-');
  EXP.fname_spm_fig = strrep(EXP.fname_spm_fig,'<','-lt-');
end
% delete previous figure files today (rewriting)
if exist(EXP.fname_spm_fig,'file')
  delete(EXP.fname_spm_fig);
end

% review the design matrix and save it:
toDisplay = {'matrix','orth','covariance'};
spm_figure('clear')
for j=1:numel(toDisplay)
  matlabbatch = {};
  matlabbatch{1}.spm.stats.review.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.review.display.(toDisplay{j})= 1;
  matlabbatch{1}.spm.stats.review.print = false;
  spm_jobman('run', matlabbatch);
  spm_print(EXP.fname_spm_fig);
end

%% 5. set contrasts and print out all results
if ~isfield(EXP,'thres')
  EXP.thres.desc  = 'cluster';
  EXP.thres.alpha = 0.05;
end
if (isfield(EXP,'NOCNTRST') && EXP.NOCNTRST) || (isfield(EXP,'NOREPORT') ...
    && EXP.NOREPORT)
else
  myspm_cntrst (EXP);
end

%% 6. Comutes contrasts quickly..?
if isfield(EXP,'convec') && isfield(EXP,'conname')
  error('mypms_mycon is depricated. just use regular contrast estimation module.')
%   EXP=myspm_mycon(EXP);
end

end

function [stdout,stderr] = myunix(cmd,~)
% [stdout,stderr] = myunix(cmd)
%
% (cc) 2017, sgKIM

% if strcmp(cmd(1:7),'tkmedit') || strcmp(cmd(1:8),'tksurfer')
%  setenv('LD_LIBRARY_PATH','glnx86:/afs/cbs.mpg.de/software/matlab/9.1/sys/os/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/bin/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64')
% else
%  setenv('LD_LIBRARY_PATH','/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64');
% end

setenv('LD_LIBRARY_PATH','/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/matlab/9.1/extern/lib/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/runtime/glnxa64:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/afs/cbs.mpg.de/software/matlab/9.1/sys/java/jre/glnxa64/jre/lib/amd64/server:/usr/lib/fsl/5.0:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/zlib1g/ubuntu-xenial-amd64:/afs/cbs.mpg.de/software/lib/matlabfix/ubuntu-xenial-amd64');
disp(['$ ',cmd]);
[stdout,stderr]=unix(cmd);
if ~isempty(stderr), disp(stderr); end
end
