function EXP = myspm_fmriprep12 (EXP)
% EXP = myspm_fmriprep12 (EXP)
%
% EXP requires:
%  .fwhm_mm
%  .fname_epi
%  .fname_t1w
%  .fname_vdm
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

% this is SPM12-included batch process using SPM12
a=spm('version');
if ~strcmp(a(4:5),'12')
  error('Run this function on SPM12!');
end

% find slice timing in msec and repetition time in sec from a example DICOM
[~,res] = mydir('/scr/vatikan1/skim/Tonotopy/main/dicom/SLET_3T/0014cmrr_mbep2d_lemon_32_rest.dcm');
EXP.fname_dcm= res{1};
hdr = spm_dicom_headers(EXP.fname_dcm);
slice_order = hdr{1}.Private_0019_1029;
TR_sec = hdr{1}.RepetitionTime/1000
EXP.TR_sec = TR_sec;
ref_slice_msec = TR_sec*1000/2        % in msec (when slice_order in given in msec)
hdr = load_nii_hdr(EXP.fname_epi);
NumFrames = hdr.dime.dim(5)

if ~isfield(EXP,'vox_mm')
  vox_mm = round(mean(hdr.dime.pixdim(2:4))*10)/10;
  EXP.vox_mm = vox_mm;
end
if ~isfield(EXP,'fwhm_mm')
  fwhm_mm = round(EXP.vox_mm*2.5*10)/10;
  EXP.fwhm_mm = fwhm_mm;
end

spm('Defaults','fmri')
spm_jobman('initcfg');
matlabbatch={};
ls(EXP.fname_epi);
for t=1:NumFrames
  matlabbatch{1}.spm.temporal.st.scans{1}{t} = [EXP.fname_epi,',',num2str(t)];
end
matlabbatch{1}.spm.temporal.st.nslices = NumFrames;
matlabbatch{1}.spm.temporal.st.tr = TR_sec;
matlabbatch{1}.spm.temporal.st.ta = 0; %TR_sec-(TR_sec/NumFrames); %     OR timing = [0 TR] when previous inputs are specified in milliseconds
matlabbatch{1}.spm.temporal.st.so = slice_order;
matlabbatch{1}.spm.temporal.st.refslice = ref_slice_msec;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

matlabbatch{2}.spm.spatial.realignunwarp.data.scans(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.realignunwarp.data.pmscan = {[EXP.fname_vdm,',1']};
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.rtm = 0;
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.einterp = 2;
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
matlabbatch{2}.spm.spatial.realignunwarp.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.jm = 0;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.sot = [];
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.rem = 1;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.noi = 5;
matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.mask = 1;
matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

ls(EXP.fname_t1w);
matlabbatch{3}.spm.spatial.preproc.channel.vols = {[EXP.fname_t1w,',1']};
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{3}.spm.spatial.preproc.channel.write = [0 1];
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,1'};
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,2'};
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,3'};
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,4'};
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,5'};
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'/scr/vatikan1/skim/matlab/spm12/tpm/TPM.nii,6'};
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 1];

matlabbatch{4}.cfg_basicio.file_dir.cfg_fileparts.files(1) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));

matlabbatch{5}.spm.util.imcalc.input(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.input(4) = cfg_dep('Segment: Bias Corrected (1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
matlabbatch{5}.spm.util.imcalc.output = 'Brain';
matlabbatch{5}.spm.util.imcalc.outdir(1) = cfg_dep('Get Pathnames: Directories (unique)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','up'));
matlabbatch{5}.spm.util.imcalc.expression = '(i1 + i2 + i3) .* i4';
matlabbatch{5}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{5}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{5}.spm.util.imcalc.options.mask = 0;
matlabbatch{5}.spm.util.imcalc.options.interp = 1;
matlabbatch{5}.spm.util.imcalc.options.dtype = 4;
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Image Calculator: Imcalc Computed Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{6}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{6}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{7}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{7}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{7}.spm.spatial.normalise.write.woptions.vox = [1 1 1]*EXP.vox_mm;
matlabbatch{7}.spm.spatial.normalise.write.woptions.interp = 4;

matlabbatch{8}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{8}.spm.spatial.smooth.fwhm = [1 1 1]*EXP.fwhm_mm;
matlabbatch{8}.spm.spatial.smooth.dtype = 0;
matlabbatch{8}.spm.spatial.smooth.im = 0;
matlabbatch{8}.spm.spatial.smooth.prefix = 's';
matlabbatch{9}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{9}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Image Calculator: Imcalc Computed Image', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{9}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
matlabbatch{9}.spm.spatial.normalise.write.woptions.vox = [1 1 1];
matlabbatch{9}.spm.spatial.normalise.write.woptions.interp = 4;

[p1,~,~] = fileparts(EXP.fname_epi);
if isempty(p1)
  p1=pwd;
end
fname_matlabbatch=[p1,'/spm12_matlabbatch_fmriprep.mat'];
save(fname_matlabbatch,'matlabbatch');
spm_jobman('run', matlabbatch);
end
