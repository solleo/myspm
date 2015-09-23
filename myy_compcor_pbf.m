function EXP = myy_compcor_pbf(EXP)

if ~isfield(EXP,'fname_dcm')
  EXP.fname_dcm = ['/scr/vatikan1/skim/Tonotopy/main/dicom/SLET_3T/',...
    '0014cmrr_mbep2d_lemon_32_rest.dcm'];
end
hdr = spm_dicom_headers(EXP.fname_dcm);
TR_sec = hdr{1}.RepetitionTime/1000
if ~isfield(EXP,'bpf')
  FilterBand = EXP.bpf;
else
  FilterBand = [0.01 0.1];
end
IsVarianceNormalization=1;

if ~isfield(EXP,'NumPCs')
  EXP.NumPCs = 6;
end

% extracting covariates from csf/wm signals after detending, filtering, and
% variance normalization
[PCs] = y_CompCor_PC(fname_func, fnames_masks, '', EXP.NumPCs, ...
  IsNeedDetrend, FilterBand, TR_sec, IsVarianceNormalization);

% gray matter: GM mean? or PC1?

end




function [PCs] = y_CompCor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR, IsVarianceNormalization)
% FORMAT [PCs] = y_CompCor_PC(ADataDir,Nuisance_MaskFilename, OutputName, PCNum, IsNeedDetrend, Band, TR, IsVarianceNormalization)
% Input:
%   ADataDir    -  The data direcotry
%   Nuisance_MaskFilename   -  The Mask file for nuisance area, e.g., the combined mask of WM and CSF
%                           -  Or can be cells, e.g., {'CSFMask','WMMask'}
%	OutputName  	-	Output filename
%   PCNum - The number of PCs to be output
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

if ~exist('CUTNUMBER','var')
    CUTNUMBER = 10;
end

if ~exist('PCNum','var')
    PCNum = 5;
end


fprintf('\nExtracting principle components for CompCor Correction:\t"%s"', ADataDir);
[AllVolume,VoxelSize,theImgFileList, Header] = y_ReadAll(ADataDir);
[nDim1 nDim2 nDim3 nDimTimePoints]=size(AllVolume);
BrainSize = [nDim1 nDim2 nDim3];

AllVolume=reshape(AllVolume,[],nDimTimePoints)';

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
    %AllVolume=detrend(AllVolume);
    fprintf('\n\t Detrending...');
    SegmentLength = ceil(size(AllVolume,2) / CUTNUMBER);
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
    fprintf('\n\t Filtering...');
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
    AllVolume = (AllVolume-repmat(mean(AllVolume),size(AllVolume,1),1))./repmat(std(AllVolume),size(AllVolume,1),1);
    AllVolume(isnan(AllVolume))=0;
end

% SVD
[U S V] = svd(AllVolume,'econ');

PCs = U(:,1:PCNum);
%Save the results
%[pathstr, name, ext] = fileparts(OutputName);
PCs = double(PCs);
%save([fullfile(pathstr,[name]), '.mat'], 'PCs')
%save([fullfile(pathstr,[name]), '.txt'], 'PCs', '-ASCII', '-DOUBLE','-TABS')
[path1,~,~] = fileparts(ADataDir);
save(fullfile(path1,['gmcsf_pcs',num2str(PCNum),'.txt']), 'PCs', '-ASCII', '-DOUBLE','-TABS')

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
    
    if prod([size(Data),length(FileList),8]) < 1024*1024*1024 %If data is with two many volumes, then it will be converted to the format 'single'.
        for j=1:length(FileList)
            [DataTemp] = y_ReadRPI(FileList{j});
            Data(:,:,:,j) = DataTemp;
        end
    else
        Data = single(Data);
        for j=1:length(FileList)
            [DataTemp] = y_ReadRPI(FileList{j});
            Data(:,:,:,j) = single(DataTemp);
        end
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

FrequencySetZero_Index = find(FrequencyMask==0);

%Remove the mean before zero padding
Data = Data - repmat(mean(Data),size(Data,1),1);

Data = [Data;zeros(paddedLength -sampleLength,size(Data,2))]; %padded with zero

Data = fft(Data);

Data(FrequencySetZero_Index,:) = 0;

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