function JOB = myy_compcor(JOB)
% JOB = myy_compcor(JOB)
% does:
%  extracting covariates from csf/wm(>.99) signals
%  and saving eigenvalues, plots, and mean gm(>.99) signal
%
% JOB requires:
%  .dir_data
%  .name_epi
%  .TR_sec
% (.param_mask) [gm_prob_thres, wm/csf_prob_thres]
% (.num_pcs) 16 (default)
%  .t1w_suffix 't1w' (default) for c?${t1w_suffix}.nii
% (.nofigure)
% asumming processes are done by myspm_fmriprep12_func.m
%
% (cc) 2015, 2016, 2018, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

path0=pwd;
if ~isfield(JOB,'prefix'), JOB.prefix='o';  end
if ~isfield(JOB,'t1w_suffix'), JOB.t1w_suffix='t1w'; end
if ~isfield(JOB,'overwrite'), JOB.overwrite=0; end
if ~isfield(JOB,'bpf'),    JOB.bpf = [0 Inf]; end
if ~isfield(JOB,'num_pcs'), JOB.num_pcs = 16; end
if ~isfield(JOB,'detrend'), JOB.detrend = 1; end
if ~isfield(JOB,'varnorm'), JOB.varnorm = 1; end

global fig_dpi
if ~isfield(JOB,'fig_dpi'), fig_dpi=300; else fig_dpi=JOB.fig_dpi; end
if ~isfield(JOB,'path1')
  path1=fullfile(JOB.dir_data);
else
  path1=JOB.path1;
end
cd(path1);
if ~isfield(JOB,'fname_epi') && isfield(JOB,'name_epi')
  JOB.fname_epi = [path1,'/',JOB.name_epi];
elseif isfield(JOB,'fname_epi') && ~isfield(JOB,'name_epi')
  [~,name1,ext1]=myfileparts(JOB.fname_epi);
  JOB.name_epi=[name1,ext1];
end
[~,name1,~]=myfileparts(JOB.fname_epi);
output_suffix=sprintf('_n%db%0.2f-%0.2f',JOB.num_pcs, JOB.bpf);
JOB.output_suffix = output_suffix;
fname_out=[name1,output_suffix,'_eigenval.txt'];

if ~exist(fname_out,'file')
  %% 1. bring segmentations to functional space
  if ~isfield(JOB,'param_mask'),JOB.param_mask=[.95 .999 .99];end
  tprob=cell(1,3);
  fnames_out=cell(1,3);
  if ~isfield(JOB,'skipcoreg')
    need2run=0;
    for c=1:3
      tprob{c} = sprintf('%0.2f',JOB.param_mask(c));
      fnames_out{c} = [name1,'_c',num2str(c),JOB.t1w_suffix,'_',tprob{c},'.nii'];
      need2run=need2run||~exist(fnames_out{c},'file');
      name_tpm{c} = ['c',num2str(c),JOB.t1w_suffix,'.nii'];
    end
    job1=[];
    job1.prefix=[name1,'_'];
    job1.fname_fixed  = [path1,'/mmean',JOB.name_epi]; % B1-bias-corrected EPI
    job1.interp=0;
    if isfield(JOB,'name_meanepi')
      job1.fname_fixed = [path1,'/',JOB.name_meanepi];
    end
    job1.fname_moving = [path1,'/bm',JOB.t1w_suffix,'.nii'];
    %  if strcmpi(JOB.t1w_suffix,'uni') % if this is MP2RAGE uniform image:
    %   job1.fname_moving = [path1,'/bm',JOB.t1w_suffix,'.nii'];
    %  end
    for c=1:3
      job1.fname_others{c} = [path1,'/',name_tpm{c}];
    end
    if need2run
      if isfield(JOB,'fn_reg') % if ANTs registration is given:
        for c=1:3
          job=[];
          job.fname_fixed = job1.fname_fixed;
          job.fname_moving = job1.fname_others{c};
          job.fname_out = [name1,'_',name_tpm{c}];
          job.interpolation = 'NearestNeighbor';
          job.transforms = {JOB.fn_reg,1}; % MNI to T1w
          myants_antsApplyTransforms(job)
        end
      else
        job1.overwrite=1;
        myspm_coreg(job1); % without changing header!
      end
    end
    %   end
    
  else
    % for EPIs already in MNI152 space
    name_tpm={'c1tpm.nii','c2tpm.nii','c3tpm.nii'};
    fnames_out={};
    for c=1:3
      unix(['ln -sf c',num2str(c),'tpm.nii ',name1,'_c',num2str(c),'tpm.nii']);
      tprob{c} = sprintf('%0.2f',JOB.param_mask(c));
      fnames_out{c} = [name1,'_c',num2str(c),JOB.t1w_suffix,'_',tprob{c},'.nii'];
    end
  end
  
  % After unwarpping, there are always NaN at image boundary.
  % For partial FOVs, it's possible that those NaN are placed in the middle
  % of the brain.
  epi=load_untouch_nii(JOB.fname_epi,1);
  epi_nan=isnan(epi.img);
  
  C=[];
  for c=1:3
    nii_c = load_uns_nii([name1,'_',name_tpm{c}]);
    nii_c.img = (nii_c.img > JOB.param_mask(c)) & ~epi_nan;
    if c==1
      C = single(~~nii_c.img);
    else
      if sum(C(~~nii_c.img(:)))
        error('GM/WM/CSF mask overlaps!');
      end
      C(~~nii_c.img)=c;
    end
    save_untouch_nii(nii_c, fnames_out{c});
  end
  if ~isfield(JOB,'skipcoreg') && ~(isfield(JOB,'nofigure') && JOB.nofigure)
    [p5,f5,e5]=fileparts(job1.fname_moving);
    
    hf=figure('visible','off');
    cfg=struct('voxsize', nii_c.hdr.dime.pixdim(2:4),'colormap',[0 0 0; eye(3).*0.8],...
      'contourcolor','w','num_contour',2);
    fixed = load_untouch_nii(job1.fname_fixed,1);
    imageorth(C,cfg, fixed.img);
    export_fig([name1,'_tpm_masks.png'],['-r',num2str(fig_dpi)])
    close(hf);
  end
  JOB.fname_masks  = fnames_out(2:3);
  if isfield(JOB,'onlycsf')
    JOB.fname_masks  = fnames_out(3);
  end
  JOB.fname_gmmask = fnames_out{1};
  JOB.tprob = tprob;
  
  %% 2. let's run
  JOB = y_CompCor_PC(JOB);
else
  JOB.fname_cc=[name1,output_suffix,'_eigenvec.txt'];
end
disp('Done')
ls(fname_out)
cd(path0);

end


function JOB = y_CompCor_PC(JOB)
% FORMAT [PCs] = y_CompCor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR, IsVarianceNormalization)
% Input:
%   ADataDir    -  The data direcotry
%   Nuisance_MaskFilename   -  The Mask file for nuisance area, e.g., the combined mask of WM and CSF
%                           -  Or can be cells, e.g., {'CSFMask','WMMask'}
%	OutputName  	-	Output filename
%   PCNum - The number of PCs to be output    ,
%   IsNeedDetrend   -   0: Dot not detrend; 1: Use Matlab's detrend
%                   -   DEFAULT: 1 -- Detrend (demean) and variance normalization will be performed before PCA, as done in Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage 37, 90-101.
%   Band            -   Temporal filter band: matlab's ideal filter e.g. [0.01 0.08]. Default: not doing filtering
%   TR              -   The TR of scanning. (Used for filtering.)
%   IsVarianceNormalization - This will perform variance normalization (subtract mean and divide by standard deviation)
%                   -   DEFAULT: 1 -- Detrend (demean) and variance normalization will be performed before PCA, as done in Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage 37, 90-101.
% Output:
%   PCs - The PCs of the nuisance area (e.g., the combined mask of WM and CSF) for CompCor correction
%__________________________________________________________________________
% Written by YAN Chao-Gan (ycg.yan@gmail.com) on 130808.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA

%[JOB.PCs, JOB.eigval, JOB.GMs] = y_CompCor_PC(JOB.fname_epi, JOB.fname_masks, '', JOB.num_pcs, ...
% JOB.detrend, FilterBand, TR_sec, JOB.varnorm, output_suffix, JOB.fname_gmmask, JOB.fname_epi, tprob, JOB);

global fig_dpi

ADataDir                = JOB.name_epi;
%fname_epi               = JOB.fname_epi;
Nuisance_MaskFilename   = JOB.fname_masks;
gm_mask                 = JOB.fname_gmmask;
%tprob                   = JOB.tprob;
PCNum                   = JOB.num_pcs;
IsNeedDetrend           = JOB.detrend;
Band                    = JOB.bpf;
TR                      = JOB.TR_sec;
IsVarianceNormalization = JOB.varnorm;
output_suffix           = JOB.output_suffix;

if ~exist('CUTNUMBER','var')
  CUTNUMBER = 20;
end

fprintf('\nExtracting principle components for CompCor Correction:\t"%s"', ADataDir);
[AllVolume,VoxelSize,theImgFileList, Header] = y_ReadAll(ADataDir);
% AllVolume=single(AllVolume);
[nDim1,nDim2,nDim3,nDimTimePoints]=size(AllVolume);
BrainSize=[nDim1 nDim2 nDim3];
AllVolume=(reshape(AllVolume,[],nDimTimePoints).');

% global (gm) signal
[GMMaskData,~,~]=y_ReadRPI(gm_mask);
GM = AllVolume(:,find(GMMaskData(:)));
imageintensity = AllVolume(:,find(GMMaskData(:)));
gm = nanmean(GM,2); % before computing signal change(%)
gm0 = repmat(nanmean(GM,1),[nDimTimePoints,1]);
GM = (GM-gm0)./(eps+gm0)*100;

% wm and csf
if ischar(Nuisance_MaskFilename)
  [MaskData,MaskVox,MaskHead]=y_ReadRPI(Nuisance_MaskFilename);
elseif iscell(Nuisance_MaskFilename)
  MaskData = 0;
  for iMask=1:length(Nuisance_MaskFilename)
    [MaskDataTemp,MaskVox,MaskHead]=y_ReadRPI(Nuisance_MaskFilename{iMask});
    MaskData = MaskData + MaskDataTemp;
  end
  MaskData = MaskData~=0;
end
MaskDataOneDim=reshape(MaskData,1,[]);
AllVolume=AllVolume(:,find(MaskDataOneDim));

% Detrend
if ~(exist('IsNeedDetrend','var') && IsNeedDetrend==0)
  %DEFAULT: 1 -- Detrend (demean) and variance normalization will be performed before PCA, as done in Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage 37, 90-101.
  fprintf('\n# Detrending...');
  AllVolume=detrend(AllVolume);
end

% Filtering
if exist('Band','var') && ~isempty(Band)
  fprintf('\n# Filtering...');
  AllVolume = y_IdealFilter(AllVolume, TR, Band);
end

%Variance normalization
if ~(exist('IsVarianceNormalization','var') && IsVarianceNormalization==0)
  %DEFAULT: 1 -- Detrend (demean) and variance normalization will be performed before PCA, as done in Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage 37, 90-101.
  AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1)) ...
    ./repmat(std(AllVolume),size(AllVolume,1),1);
  AllVolume(isnan(AllVolume))=0;
end
[path1,~,~] = fileparts(ADataDir);
[~,name1,~]=fileparts(JOB.name_epi);
if PCNum
  % SVD
  [U,S,V] = svd(AllVolume,'econ');
  eigval = diag(S);
  xvar=cumsum(eigval)/sum(eigval)*100;
  %% SCREE PLOT
  if ~isfield(JOB,'nofigure')
    
    hf=figure('position',[2237         234         560         634],'visible','off');
    subplot(311)
    plot(eigval,'k'); xlabel('Order of eigenvalues'); ylabel('Eigenvalue')
    ylim0=ylim; hold on; line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    %title([num2str(eigval(PCNum)),'@',num2str(PCNum),'-th PC of WM/CSF voxels'])
    title('BOLD timeseries in WM/CSF voxels')
    xlim([1 50]);
    text(PCNum,eigval(PCNum), ...
      ['(',num2str(PCNum),', ',num2str(eigval(PCNum),3),')'], ...
      'color','r');
    
    subplot(312)
    ddy=gaussblur([0; 0; diff(diff(eigval))],3);
    upto = min([50 numel(ddy)]);
    plot(ddy(1:upto),'k'); xlabel('Order of eigenvalues');
    ylabel('Change of slope');
    ylim0=ylim; hold on;
    line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    title('Smoothed with FWHM = 3 orders');
    xlim([1 50]);
    [~,jj]=max(ddy(1:upto));
    line([jj jj]', [ylim0(1) ylim0(2)]','color','b');
    text(PCNum,ylim0(1), ...
      ['(',num2str(PCNum),', ',num2str(ddy(PCNum),3),')'], ...
      'color','r','fontsize',15);
    text(jj,ylim0(2), ...
      ['(',num2str(jj),', ',num2str(ddy(jj),3),')'], ...
      'color','b','fontsize',15);
    
    subplot(313);
    plot(xvar,'k'); xlabel('Order of eigenvalues');
    ylabel({'Cumulative','expalined','variance(%)'})
    ylim0=ylim; hold on; line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    line([jj jj]', [ylim0(1) ylim0(2)]','color','b');
    
    %title([num2str(xvar(PCNum)),'% with ',num2str(PCNum),'PCs',])
    text(PCNum,ylim0(1), ...
      ['(',num2str(PCNum),', ',num2str(xvar(PCNum),3),')'], ...
      'color','r','fontsize',15);
    text(jj,ylim0(2), ...
      ['(',num2str(jj),', ',num2str(xvar(jj),3),')'], ...
      'color','b','fontsize',15);
    xlim([1 50]);
    [~,name1,~]=fileparts(JOB.name_epi);
    %     screen2png(fullfile(path1,[name1,'_screeplot',output_suffix,'.png']),fig_dpi);
    export_fig(fullfile(path1,[name1,'_screeplot',output_suffix,'.png']),['-r',num2str(fig_dpi)]);
    close(hf);
  end
  PCs = U(:,1:PCNum);
  PCs = double(PCs);
  save(fullfile(path1,[name1,output_suffix,'_eigenval.txt']), ...
    'eigval', '-ASCII', '-DOUBLE','-TABS')
else
  % mean
  PCs = mean(AllVolume,2);
  eigval = [];
end
[~,name1,~]=fileparts(JOB.name_epi);
save(fullfile(path1,[name1,output_suffix,'_eigenvec.txt']), ...
  'PCs', '-ASCII', '-DOUBLE','-TABS');
JOB.fname_cc=[name1,output_suffix,'_eigenvec.txt'];
save(fullfile(path1,[name1,'_gm.txt']), 'gm', '-ASCII', '-DOUBLE','-TABS');
JOB.fname_gs=fullfile(path1,[name1,'_gm.txt']);
fprintf('\nFinished Extracting principle components for CompCor Correction.\n');
end

function [Data, VoxelSize, FileList, Header] = y_ReadAll(InputName)
%function [Data, VoxelSize, FileList, Header] = y_ReadAll(InputName)
% Read NIfTI files in all kinds of input formats.
% Will call y_ReadRPI.m, which reads a single file.
% ------------------------------------------------------------------------
% Input:
% InputName - Could be the following format:
%                  1. A single file (.img/hdr, .nii, or .nii.gz), give the path and filename.
%                  2. A directory, under which could be a single 4D file, or a set of 3D images
%                  3. A Cell (nFile * 1 cells) of filenames of 3D image file, or a single file of 4D NIfTI file.
% Output:
% Data - 4D matrix of image data. (If there is no rotation in affine matrix, then will be transformed into RPI orientation).
% VoxelSize - the voxel size
% FileList - the list of files
% Header - a structure containing image volume information (as defined by SPM, see spm_vol.m)
% The elements in the structure are:
%       Header.fname - the filename of the image.
%       Header.dim   - the x, y and z dimensions of the volume
%       Header.dt    - A 1x2 array.  First element is datatype (see spm_type).
%                 The second is 1 or 0 depending on the endian-ness.
%       Header.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       Header.pinfo - plane info for each plane of the volume.
%              Header.pinfo(1,:) - scale for each plane
%              Header.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*Header.pinfo(1,j) + Header.pinfo(2,j)
%              Header.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%__________________________________________________________________________
% Written by YAN Chao-Gan 130624.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com
global fig_dpi

if iscell(InputName)
  if size(InputName,1)==1
    InputName=InputName';
  end
  FileList = InputName;
elseif (7==exist(InputName,'dir'))
  DirImg=dir(fullfile(InputName,'*.img'));
  if isempty(DirImg)
    DirImg=dir(fullfile(InputName,'*.nii.gz'));
  end
  if isempty(DirImg)
    DirImg=dir(fullfile(InputName,'*.nii'));
  end
  
  FileList={};
  for j=1:length(DirImg)
    FileList{j,1}=fullfile(InputName,DirImg(j).name);
  end
elseif (2==exist(InputName,'file'))
  FileList={InputName};
else
  error(['The input name is not supported by y_ReadAll: ',InputName]);
end

fprintf('\nReading images from "%s" etc.\n', FileList{1});

if length(FileList) == 0
  error(['No image file is found for: ',InputName]);
elseif length(FileList) == 1
  [Data, VoxelSize, Header] = y_ReadRPI(FileList{1});
elseif length(FileList) > 1 % A set of 3D images
  [Data, VoxelSize, Header] = y_ReadRPI(FileList{1});
  Data = zeros([size(Data),length(FileList)],'single');
  for j=1:length(FileList)
    [DataTemp] = y_ReadRPI(FileList{j});
    Data(:,:,:,j) = single(DataTemp);
  end
  
end
end

function [Data_Filtered] = y_IdealFilter(Data, SamplePeriod, Band)
% FORMAT    [Data_Filtered] = y_IdealFilter(Data, SamplePeriod, Band)
% Input:
% 	Data		    -	2D data matrix (nDimTimePoints * nTimeSeries)
% 	SamplePeriod	-   Sample period, i.e., 1/sample frequency. E.g., TR
%   Band            -   The frequency for filtering, 1*2 Array. Could be:
%                   [LowCutoff_HighPass HighCutoff_LowPass]: band pass filtering
%                   [0 HighCutoff_LowPass]: low pass filtering
%                   [LowCutoff_HighPass 0]: high pass filtering
% Output:
%	Data_Filtered       -   The data after filtering
%-----------------------------------------------------------
% Written by YAN Chao-Gan 120504 based on REST.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com

global fig_dpi

sampleFreq 	 = 1/SamplePeriod;
sampleLength = size(Data,1);
paddedLength = 2^nextpow2(sampleLength);
LowCutoff_HighPass = Band(1);
HighCutoff_LowPass = Band(2);

% Get the frequency index
if (LowCutoff_HighPass >= sampleFreq/2) % All high stop
  idxLowCutoff_HighPass = paddedLength/2 + 1;
else % high pass, such as freq > 0.01 Hz
  idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1);
end

if (HighCutoff_LowPass>=sampleFreq/2)||(HighCutoff_LowPass==0) % All low pass
  idxHighCutoff_LowPass = paddedLength/2 + 1;
else % Low pass, such as freq < 0.08 Hz
  idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1);
end

FrequencyMask = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1) = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;

%Remove the mean before zero padding
Data = Data - repmat(mean(Data),size(Data,1),1);

Data = [Data;zeros(paddedLength -sampleLength,size(Data,2))]; %padded with zero

Data = fft(Data);

Data(FrequencyMask==0,:) = 0;

Data = ifft(Data);

Data_Filtered = Data(1:sampleLength,:);
end

function [Data, VoxelSize, Header] = y_ReadRPI(FileName, VolumeIndex)
%function [Data, VoxelSize, Header] = y_ReadRPI(FileName, VolumeIndex)
% Read NIfTI image in RPI orientation -- for NIfTI files without rotation in affine matrix!!!
% Will call y_Read.m, which does not adjust orientation.
% ------------------------------------------------------------------------
% Input:
% FileName - the path and filename of the image file (*.img, *.hdr, *.nii, *.nii.gz)
% VolumeIndex - the index of one volume within the 4D data to be read, can be 1,2,..., or 'all'.
%               default: 'all' - means read all volumes
% Output:
% Data - 3D or 4D matrix of image data in RPI orientation (if there is no rotation in affine matrix).
% VoxelSize - the voxel size
% Header - a structure containing image volume information (as defined by SPM, see spm_vol.m)
% The elements in the structure are:
%       Header.fname - the filename of the image.
%       Header.dim   - the x, y and z dimensions of the volume
%       Header.dt    - A 1x2 array.  First element is datatype (see spm_type).
%                 The second is 1 or 0 depending on the endian-ness.
%       Header.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       Header.pinfo - plane info for each plane of the volume.
%              Header.pinfo(1,:) - scale for each plane
%              Header.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*Header.pinfo(1,j) + Header.pinfo(2,j)
%              Header.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%__________________________________________________________________________
% Written by YAN Chao-Gan 130624.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com

global fig_dpi

if ~exist('VolumeIndex', 'var')
  VolumeIndex='all';
end

[Data,Header] = y_Read(FileName,VolumeIndex);

if sum(sum(Header.mat(1:3,1:3)-diag(diag(Header.mat(1:3,1:3)))~=0))==0 % If the image has no rotation (no non-diagnol element in affine matrix), then transform to RPI coordination.
  if Header.mat(1,1)>0 %R
    Data = flipdim(Data,1);
    Header.mat(1,:) = -1*Header.mat(1,:);
  end
  if Header.mat(2,2)<0 %P
    Data = flipdim(Data,2);
    Header.mat(2,:) = -1*Header.mat(2,:);
  end
  if Header.mat(3,3)<0 %I
    Data = flipdim(Data,3);
    Header.mat(3,:) = -1*Header.mat(3,:);
  end
end
temp = inv(Header.mat)*[0,0,0,1]';
Header.Origin = temp(1:3)';

VoxelSize = sqrt(sum(Header.mat(1:3,1:3).^2));
end

function [Data, Header] = y_Read(FileName, VolumeIndex)
%function [Data, Header] = y_Read(FileName, VolumeIndex)
% Read NIfTI file Based on SPM's nifti
% ------------------------------------------------------------------------
% Input:
% FileName - the path and filename of the image file (*.img, *.hdr, *.nii, *.nii.gz)
% VolumeIndex - the index of one volume within the 4D data to be read, can be 1,2,..., or 'all'.
%               default: 'all' - means read all volumes
% Output:
% Data - 3D or 4D matrix of image data
% Header - a structure containing image volume information (as defined by SPM, see spm_vol.m)
% The elements in the structure are:
%       Header.fname - the filename of the image.
%       Header.dim   - the x, y and z dimensions of the volume
%       Header.dt    - A 1x2 array.  First element is datatype (see spm_type).
%                 The second is 1 or 0 depending on the endian-ness.
%       Header.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       Header.pinfo - plane info for each plane of the volume.
%              Header.pinfo(1,:) - scale for each plane
%              Header.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*Header.pinfo(1,j) + Header.pinfo(2,j)
%              Header.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%__________________________________________________________________________
% Written by YAN Chao-Gan 130624 based on SPM's NIfTI.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com

global fig_dpi

if ~exist('VolumeIndex', 'var')
  VolumeIndex='all';
end

[pathstr, name, ext] = fileparts(FileName);

if isempty(ext)
  FileName = fullfile(pathstr,[name '.nii']);
  if ~exist(FileName,'file')
    FileName = fullfile(pathstr,[name '.hdr']);
  end
  if ~exist(FileName,'file')
    FileName = fullfile(pathstr,[name '.nii.gz']);
    [pathstr, name, ext] = fileparts(FileName);
  end
end

if ~exist(FileName,'file')
  error(['File doesn''t exist: ',fullfile(pathstr,[name ext])]);
end

FileNameWithoutGZ = FileName;
if strcmpi(ext,'.gz')
  gunzip(FileName);
  FileName = fullfile(pathstr,[name]);
end

Nii  = nifti(FileName);
V = spm_vol(FileName);

if(~strcmpi(VolumeIndex,'all'))
  Data = squeeze(double(Nii.dat(:,:,:,VolumeIndex)));
  Header = V(VolumeIndex);
else
  Data = double(Nii.dat);
  Header = V(1);
end
Header.fname=FileNameWithoutGZ;

if strcmpi(ext,'.gz')
  delete(FileName);
end
end
