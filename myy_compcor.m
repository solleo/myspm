function EXP = myy_compcor(EXP)
% EXP = myy_compcor(EXP)
% does:
%  extracting covariates from csf/wm(>.99) signals
%  and saving eigenvalues, plots, and mean gm(>.99) signal
%
% EXP requires:
%  .dir_base
%  .subjID
%  .name_epi
% (.dir_figure)
% (.param_mask) [gm_prob_thres, wm/csf_prob_thres]
% (.num_pcs) 16 (default)
% (.t1w_suffix) 't1w' (default)
% (.nofigure)
% asumming processes are done by myspm_fmriprep12_func.m
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

path0=pwd;
if ~isfield(EXP,'prefix'), EXP.prefix='o';  end
if ~isfield(EXP,'t1w_suffix'), EXP.t1w_suffix='t1w'; end
subjID = fsss_subjID(EXP.subjID);

for i=1:numel(subjID)
  subjid = subjID{i};
  EXP.subjid = subjid;
  path1=fullfile(EXP.dir_base,subjid);
  cd(path1);
  EXP.fname_epi = fullfile(EXP.dir_base, subjid, EXP.name_epi);
  
  %% 1. bring segmentations to functional space
  if ~isfield(EXP,'param_mask')
    EXP.param_mask=[.99 .99];
  end
  if numel(EXP.param_mask) < 2
    error('EXP.param_mask: Enter two prob threshold for GM and [WM/CSF]');
  end
  EXP.param_mask(3) = EXP.param_mask(2);
  
  need2run=0;
  tprob=cell(1,3);
  fnames_out=cell(1,3);
  for c=1:3
    tprob{c} = sprintf('%0.2f',EXP.param_mask(c));
    fnames_out{c} = ['oc',num2str(c),EXP.t1w_suffix,'_',tprob{c},'.nii'];
    need2run=need2run||~exist(fnames_out{c},'file');
    name_tpm{c} = ['c',num2str(c),EXP.t1w_suffix,'.nii'];
  end
  if need2run
    if ~exist(['o',name_tpm{1}],'file')||...
       ~exist(['o',name_tpm{2}],'file')||...
       ~exist(['o',name_tpm{3}],'file')
      exp1=[];
      exp1.name_fixed  = ['mean',EXP.name_epi];
      if isfield(EXP,'name_meanepi')
        exp1.name_fixed = EXP.name_meanepi;
      end
      if isfield(EXP,'name_t1w')
        exp1.name_moving = 'Brain.nii';
      else
        exp1.name_moving = EXP.name_t1w;
      end
      exp1.name_others = name_tpm;
      myspm_coreg(exp1);
    end
    for c=1:3
      unix(['fslmaths o',name_tpm{c},' ',...
        ' -thr ',tprob{c},' -bin ',fnames_out{c},' -odt char']);
    end
  end
  EXP.fname_masks  = fnames_out(2:3);
  EXP.fname_gmmask = fnames_out{1};
  EXP.tprob = tprob;
  
  %% 2. let's run
  if ~isfield(EXP,'bpf1'),    EXP.bpf1 = [0 Inf]; end
  if ~isfield(EXP,'num_pcs'), EXP.num_pcs = 16; end
  if ~isfield(EXP,'detrend'), EXP.detrend = 1; end
  if ~isfield(EXP,'varnorm'), EXP.varnorm = 1; end
  output_suffix=sprintf('_n%dd%dv%db%0.2f-%0.2f', ...
    EXP.num_pcs, EXP.detrend, EXP.varnorm, EXP.bpf1);
  EXP.output_suffix = output_suffix;
  EXP = y_CompCor_PC(EXP);
  
  % copy figures
  if isfield(EXP,'dir_figure') && ~isfield(EXP,'nofigure')
    [~,~]=mkdir(EXP.dir_figure);
    copyfile([path1,'/cc_plot',output_suffix,'.png'], ...
      [EXP.dir_figure,'/cc_plot',output_suffix,'_',subjid,'.png']);
    if exist([path1,'/cc_eigval',output_suffix,'.png'],'file')
      copyfile([path1,'/cc_eigval',output_suffix,'.png'], ...
        [EXP.dir_figure,'/cc_eigval',output_suffix,'_',subjid,'.png']);
    end
  end
  
end
cd(path0);

end


function EXP = y_CompCor_PC(EXP)
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

%[EXP.PCs, EXP.eigval, EXP.GMs] = y_CompCor_PC(EXP.fname_epi, EXP.fname_masks, '', EXP.num_pcs, ...
% EXP.detrend, FilterBand, TR_sec, EXP.varnorm, output_suffix, EXP.fname_gmmask, EXP.fname_epi, tprob, EXP);

ADataDir                = EXP.name_epi;
fname_epi               = EXP.fname_epi;
Nuisance_MaskFilename   = EXP.fname_masks;
gm_mask                 = EXP.fname_gmmask;
tprob                   = EXP.tprob;
PCNum                   = EXP.num_pcs;
IsNeedDetrend           = EXP.detrend;
Band                    = EXP.bpf1;
TR                      = EXP.TR_sec;
IsVarianceNormalization = EXP.varnorm;
output_suffix           = EXP.output_suffix;

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
gm = mean(GM,2); % before computing signal change(%)
gm0 = repmat(mean(GM,1),[nDimTimePoints,1]);
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
  SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER); % spatial segment (columns of AllVolume)
  for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
      Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
      Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = detrend(AllVolume(:,Segment));
    fprintf('.');
  end
end

% Filtering
if exist('Band','var') && ~isempty(Band)
  fprintf('\n# Filtering...');
  SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
  for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
      Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
      Segment = (iCut-1)*SegmentLength+1 : size(AllVolume,2);
    end
    AllVolume(:,Segment) = y_IdealFilter(AllVolume(:,Segment), TR, Band);
    fprintf('.');
  end
end

%Variance normalization
if ~(exist('IsVarianceNormalization','var') && IsVarianceNormalization==0)
  %DEFAULT: 1 -- Detrend (demean) and variance normalization will be performed before PCA, as done in Behzadi, Y., Restom, K., Liau, J., Liu, T.T., 2007. A component based noise correction method (CompCor) for BOLD and perfusion based fMRI. Neuroimage 37, 90-101.
  AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1)) ...
    ./repmat(std(AllVolume),size(AllVolume,1),1);
  AllVolume(isnan(AllVolume))=0;
end
[path1,~,~] = fileparts(ADataDir);
if PCNum
  % SVD
  [U,S,V] = svd(AllVolume,'econ');
  eigval = diag(S);
  xvar=cumsum(eigval)/sum(eigval)*100;
  if ~isfield(EXP,'nofigure')
    hf=figure('position',[2237         234         560         634]);
    subplot(311)
    plot(eigval,'b'); xlabel('Order of eigenvalues'); ylabel('Eigenvalue')
    ylim0=ylim; hold on; line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    title([num2str(eigval(PCNum)),'@',num2str(PCNum),'-th PC'])
    xlim([1 50]);
    
    % "Scree plot" method
    subplot(312)
    ddy=gaussblur([0; 0; diff(diff(eigval))],3);
    plot(ddy(1:50),'b'); xlabel('Order of eigenvalues'); 
    ylabel({'Smoothed (fwhm=3)','change of slope'});
    ylim0=ylim; hold on; line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    title([num2str(eigval(PCNum)),'@',num2str(PCNum),'-th PC'])
    xlim([1 50]);
    
    subplot(313);
    plot(xvar,'b'); xlabel('Order of eigenvalues'); 
    ylabel('Cumulative expalined variance(%)')
    ylim0=ylim; hold on; line([PCNum,PCNum]', [ylim0(1) ylim0(2)]','color','r');
    title([num2str(xvar(PCNum)),'% with ',num2str(PCNum),'PCs',])
    xlim([1 50]);
    
    screen2png(['cc_eigval',output_suffix,'.png']);
    close(hf);
  end
  PCs = U(:,1:PCNum);
  PCs = double(PCs);
  save(fullfile(path1,['cc_wmcsf',tprob{3},output_suffix,'_eigenval.txt']), ...
    'eigval', '-ASCII', '-DOUBLE','-TABS')
else
  % mean
  PCs = mean(AllVolume,2);
  eigval = [];
end
save(fullfile(path1,['cc_wmcsf',tprob{2},output_suffix,'_eigenvec.txt']), ...
  'PCs', '-ASCII', '-DOUBLE','-TABS')
save(fullfile(path1,['cc_gm',tprob{1},'.txt']), 'gm', '-ASCII', '-DOUBLE','-TABS')
fprintf('\nFinished Extracting principle components for CompCor Correction.\n');

[~,name1,~] = fileparts(EXP.name_epi);
[~,res] = mydir(fullfile(path1,['art_regression_outliers_and_movement_',name1,'*']));
if isempty(res)
  exp=EXP;
  exp.subjID = {EXP.subjid};
  myspm_art(exp);
end
if ~isfield(EXP,'nofigure')
  
  [~,res] = mydir(fullfile(path1,['art_regression_outliers_and_movement_',name1,'*']));
  load (res{1},'R');
  load(fullfile(path1,['cc_gm',tprob{1},'.txt']), 'gm');
  load(fullfile(path1,['cc_wmcsf',tprob{2},output_suffix,'_eigenvec.txt']),'PCs');
  
  % create figure
  hf=figure('position',[2237         168         706        1009]);
  rp=R(:,end-6:end-1);
  subplot(611); plot(R(:,end)); ylabel('mov_{art}(mm)');
  ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
  title(fname_epi);
  subplot(6,1,[2:4]); imagesc(GM'); ylabel(['GM>',tprob{1},' voxels']);
  set(gca,'ydir','nor'); caxis([-5 5]); hb=colorbar; ylabel(hb,'Change from mean(%)');
  colormap(sgcolormap('CKM'));
  subplot(615); plot(gm); ylabel(['mean GM>',tprob{1}]);
  ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
  if PCNum
    subplot(616); plot(PCs(:,[1:min(6,PCNum)]));
    ylabel({['Top ',num2str(PCNum),' PCs'], ['from wm+csf>',tprob{2}]});
  else
    subplot(616); plot(PCs);
    ylabel({'Mean wm+csf',['>',tprob{2}]});
  end
  ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
  title(sprintf('BPF=[%0.2f,%0.2f] Hz',Band));
  xlabel('TR');
  screen2png(['cc_plot',output_suffix,'.png']);
  close(hf);
  
  if isfield(EXP,'imageintensity')
    hf=figure('position',[2237         168         706        1009]);
    load('rp_arest410.txt');
    subplot(611); plot([0; l2norm(diff(rp_arest410))]); ylabel('||dm/dt||_2(mm)');
    ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
    title(fname_epi)
    subplot(6,1,[2:4]); imagesc(imageintensity'); ylabel(['GM>',tprob{1},' voxels']);
    colormap hot; set(gca,'ydir','nor'); hb=colorbar; ylabel(hb,'Image intensity');
    subplot(615); plot(gm); ylabel(['mean GM>',tprob{1}]);
    ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
    if PCNum
      subplot(616); plot(PCs(:,[1:6]));
      ylabel({['Top ',num2str(PCNum),' PCs'], ['from wm+csf>',tprob{2}]});
    else
      subplot(616); plot(PCs);
      ylabel({'Mean wm+csf',['>',tprob{2}]});
    end
    ha=colorbar; set(ha,'visible','off'); xlim([1 nDimTimePoints]);
    title(sprintf('BPF=[%0.2f,%0.2f] Hz',Band));
    xlabel('TR');
    screen2png(['cc_imgval_plot',output_suffix,'.png']);
    close(hf);
  end
end
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
  Data = zeros([size(Data),length(FileList)]);
  
  %   if prod([size(Data),length(FileList),8]) < 1024*1024*1024 %If data is with two many volumes, then it will be converted to the format 'single'.
  %     for j=1:length(FileList)
  %       [DataTemp] = y_ReadRPI(FileList{j});
  %       Data(:,:,:,j) = DataTemp;
  %     end
  %   else
  Data = single(Data);
  for j=1:length(FileList)
    [DataTemp] = y_ReadRPI(FileList{j});
    Data(:,:,:,j) = single(DataTemp);
  end
  %   end
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