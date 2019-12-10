function myspm_wds(fname_epi_orig)
% JOB requires:
%
% .fname_epi_orig


%{ 
Let's try Patel's way first as described in
Patel & Bullmore, 2016, "A wavelet-based estimator of the degrees of
freedom in denoised fMRI time series for probabilistic testing of
functional connectivity and brain graphs", NeuroImage.

Core Image Processing included the following steps: 
(i) slice ac- quisition correction; 
(ii) rigid-body head movement correction to the first frame of data; 
(iii) obliquity transform to the structural image; 
(iv) affine co-registration to the skull-stripped structural image using 
a gray matter mask; 
(v) standard space transform to the MNI152 template in MNI space; 
(v) spatial smoothing (6 mm full width at half maximum); 
(vi) a within-run intensity normalization to a whole-brain median of 1000. 
All spatial transforms were applied in one step to avoid incremental blurring 
of the data that can occur from multiple independent transforms.

Denoising steps included: 
(vii) wavelet despiking (performed voxel- wise with the BrainWavelet Toolbox); 
(viii) confound signal regression including the 6 motion parameters estimated 
in (ii), their first order temporal derivatives, and ventricular cerebrospinal 
fluid (CSF) signal; and 
(ix) a wavelet “band-pass” filter. In this last step, the Maximal Over- 
lap Discrete Wavelet Transform (MODWT) was used to produce a set of scales 
(frequency bands), and coefficients from scales representing fre- 
quency bands of interest were recomposed to produce frequency- 
filtered time series.

%}

[p1,f1,~]=myfileparts(fname_epi_orig);
epi_suffix=f1;
cd(p1)

% 1. (vii) WDS
fname_epi=['ua',epi_suffix,'.nii'];
fname_wds=['WD',fname_epi];
if ~exist(fname_wds,'file')
 WaveletDespike(fname_epi, ['ua',epi_suffix]);
 unix(['gunzip ua',epi_suffix,'*.nii.gz']);
 unix(['mv ua',epi_suffix,'_wds.nii ',fname_wds]);
 unix(['mv ua',epi_suffix,'_EDOF.nii WDua',epi_suffix,'_EDOF.nii']);
 unix(['mv ua',epi_suffix,'_noise.nii WDua',epi_suffix,'_noise.nii']);
 unix(['mv ua',epi_suffix,'_SP.txt WDua',epi_suffix,'_SP.txt']);
end
ls(fname_wds)

% 2. need to swap matrices and change vox-to-world matrix...
job1=[];
job1.fname_moving = fname_wds;
job1.fname_fixed  = ['meanua',epi_suffix,'.nii'];
job1.fname_other  = ['WDua',epi_suffix,'_EDOF.nii'];
job1.interp = 0;
myspm_coreg4d(job1)
fname_wds=['o',fname_wds];
setenv('FSLOUTPUTTYPE','NIFTI')
unix(['fslmaths ',fname_wds,' -Tmean mean',fname_wds]);

% 3. Parcellation
[dir_subj,~,~]=myfileparts(p1);
%fsss_native_to_BN(dir_subj)
fsss_native_to_BN_vol(dir_subj)


% 4. (viii) regression 
% 4-1. motion parameters (6+2)
cd(p1)
fname_rp=['rp_a',epi_suffix,'.txt'];
rp=load(fname_rp);
rp=[...
 rp(:,1:3), l2norm([0 0 0; diff(rp(:,1:3))]) , ...
 rp(:,4:6), l2norm([0 0 0; diff(rp(:,4:6))]) ];
% 4-2. CSF regressors (1) 
job2=[];
job2.dir_data=p1;
job2.fname_epi=fname_wds;
hdr=load_untouch_header_only([epi_suffix,'.nii']);
job2.TR_sec=hdr.dime.pixdim(5);
disp(['TR = ',num2str(job2.TR_sec),' sec']);
job2.num_pcs=1;
job2.onlycsf=1;
job2=myy_compcor(job2);
cc=load(job2.fname_cc);
% 4-3. get residual
load([dir_subj,'/BN_vol/ME.mat'])
load([dir_subj,'/BN_vol/WD.mat'])
M=[ones(size(ME(:,1))) rp cc];
rWD = WD-M*(pinv(M'*M)*M'*WD);
maxT=size(M,1);
save([dir_subj,'/BN_vol/r1WD.mat'],'rWD')

%%
%figure
M=[ones(size(ME(:,1))) rp cc];
rWD = WD-M*(pinv(M'*M)*M'*WD);
clf
subplot(421)
imagesc(zscore(M(:,2:end))'); caxis([-3 3])
title('Move8 + CSF');
ylabel('Nodes'); xlabel('TR')
subplot(423)
imagesc(zscore(ME)'); caxis([-3 3])
title('Raw'); 
subplot(425)
imagesc(zscore(WD)'); caxis([-3 3])
title('Raw:WD')
subplot(427)
imagesc(zscore(rWD)'); caxis([-3 3])
title('Raw:WD:reg')

subplot(422)
plot(l2norm(zscore(M(:,2:end)))); xlim([1 maxT])
ylabel('L_2norm'); xlabel('TR');
title('Move8 + CSF')
subplot(424)
plot(l2norm(zscore(ME))); xlim([1 maxT])
title('Raw');
subplot(426)
plot(l2norm(zscore(WD))); xlim([1 maxT])
title('Raw:WD')
subplot(428)
plot(l2norm(zscore(rWD))); xlim([1 maxT])
title('Raw:WD:reg')

ha=axes('position',[0.02 0.75 0.05 0.2]);
hcb=colorbar(ha,'location','east');
title(hcb,' Z-socre'); caxis([-3 3])
axis off
%%
vol=zeros(246,1,1,239);
for t=1:239
 vol(:,1,1,t)=ME(t,:);
end
nii=make_nii(vol);
save_nii(nii,'ME.nii')
WaveletDespike('ME.nii','_')
nii=load_nii('__wds.nii.gz');
WD_ME=squeeze(nii.img);

%%
%{
rME = ME-M*(pinv(M'*M)*M'*ME);
vol=zeros(246,1,1,239);
for t=1:239
 vol(:,1,1,t)=rME(t,:);
end
nii=make_nii(vol);
save_nii(nii,'rME.nii')
WaveletDespike('rME.nii','_')
nii=load_nii('__wds.nii.gz');
WDrME=squeeze(nii.img);
figure
subplot(421)
imagesc(zscore(M(:,2:end))'); caxis([-3 3])
title('Move8 + CSF');
ylabel('Nodes'); xlabel('TR')
subplot(423)
imagesc(zscore(ME)'); caxis([-3 3])
title('Raw'); 
subplot(425)
imagesc(zscore(rME)'); caxis([-3 3])
title('Raw:regD')
subplot(427)
imagesc(zscore(WDrME)'); caxis([-3 3])
title('Raw:reg:WD')

subplot(422)
plot(l2norm(zscore(M(:,2:end)))); xlim([1 maxT])
ylabel('L_2norm'); xlabel('TR');
title('Move8 + CSF')
subplot(424)
plot(l2norm(zscore(ME))); xlim([1 maxT])
title('Raw');
subplot(426)
plot(l2norm(zscore(rME))); xlim([1 maxT])
title('Raw:reg')
subplot(428)
plot(l2norm(zscore(WDrME))); xlim([1 maxT])
title('Raw:reg:WD')

ha=axes('position',[0.02 0.75 0.05 0.2]);
hcb=colorbar(ha,'location','east');
title(hcb,' Z-socre'); caxis([-3 3])
axis off
%}


%% 5-1. Computing WCC
cd(dir_subj)
load ([dir_subj,'/BN_vol/ME.mat']);
[WCC,ncoef] = wavcorr(ME, 'warning',0, 'mincoef',100);
WCC_unc=WCC;
ncoef(~ncoef)=[];
J=numel(ncoef); % Let's see WCC upto scale 5:
% 5-2. Correction of effective DOF:
load ([dir_subj,'/BN_vol/eDOF.mat']);
Z=WCC*0;
for j=1:J
 % pair-wise df is defined by minimum of df's of nodes (NI(2016)142:p.20)
 D=zeros(246);
 for a=1:246
  for b=1:246
   D(a,b)=min(eDOF(j,a), eDOF(j,b));
  end
 end
 yeta0=mode(D(:))
 if yeta0>3
 for a=1:245
  for b=(a+1):246
   z = atanh(WCC(a,b,j))*sqrt(D(a,b)-3);
   Z(a,b,j)=z;
   Z(b,a,j)=z;
   WCC(a,b,j) = tanh(z/sqrt(yeta0-3));
   WCC(b,a,j) = WCC(a,b,j);
  end
 end
 else
  J=J-1;
  continue
 end
end
WCC(:,:,(J+1):end)=[];
ncoef((J+1):end)=[];

T=readtable('~/Dropbox/BN_Atlas/BN_Atlas_246_COT.xlsx');
band_mHz=zeros(size(WCC,3),2);
idx=[1:2:246 246:-2:2];
figure('position',[8         632        1218         401])
ax1=axeslayout1(8,[2 4],[.1 0]);
ax1.y=ax1.y-0.05;
ax2=ax1;
ax2.w=0.005;
ax2.x=ax2.x+0.011;
ax3=axeslayout1(8,[2 4],[.1 0.22]);
ax3.w=0.15;
ax3.x=ax3.x+0.03;
for j=1:J
 ha1=axespos(ax1,j);
 hold on
 imagesc(WCC(idx,idx,j)); axis image; caxis([-1 1]);
 band_mHz(j,:)=[2^-(j+1), 2^-j]/job2.TR_sec*1000;
 title({['Scale ',num2str(j)],[num2str(round(band_mHz(j,1))),'-',...
  num2str(round(band_mHz(j,2))),' mHz']});
 set(gca,'xtick','','ytick','')
  idx_k=find(diff(T.Label_lobe(idx)));
 idx_k(1)=[];
 for k=1:numel(idx_k)
  line([idx_k(k);idx_k(k)]+0.5,[0.5;246.5],'color','k')
  line([0.5;246.5],[idx_k(k);idx_k(k)]+0.5,'color','k')
 end
 colormap(ha1,bipolar);
 
 ha2=axespos(ax2,j);
 imagesc(T.Label_lobe(idx))
 axis on; set(gca,'xtick','','ytick','');
 colormap(ha2, sgcolormap('hsv-soft'));
 
 axespos(ax3,j+J)
 hold on;
 x=triu(WCC(idx,idx,j).*~eye(246));
 [ci,xi]=hist(x(~~x(:)),100);
 bar(xi(xi<0), ci(xi<0), 1,'facecolor',[.4 .4 1],'linestyle','none');
 bar(xi(xi>0), ci(xi>0), 1,'facecolor',[1 .4 .4],'linestyle','none');
 xlim([-1 1]); xlabel(['WCC:',num2str(j)]); ylabel('# of edges');
end
colormap(bipolar)
%%
screen2png([dir_subj,'/BN_vol/ME_WCC.png'],200)
close(gcf)
save ([dir_subj,'/BN_vol/WCC.mat'],'WCC','ncoef','band_mHz');
%%
figure('position',[8         632        1218         401])
ax1=axeslayout1(8,[2 4],[.1 0]);
ax1.y=ax1.y-0.05;
ax2=ax1;
ax2.w=0.005;
ax2.x=ax2.x+0.011;
ax3=axeslayout1(8,[2 4],[.1 0.22]);
ax3.w=0.15;
ax3.x=ax3.x+0.03;
for j=1:J
 ha1=axespos(ax1,j);
 hold on
 imagesc(WCC_unc(idx,idx,j)); axis image; caxis([-1 1]);
 band_mHz(j,:)=[2^-(j+1), 2^-j]/job2.TR_sec*1000;
 title({['Scale ',num2str(j)],[num2str(round(band_mHz(j,1))),'-',...
  num2str(round(band_mHz(j,2))),' mHz']});
 set(gca,'xtick','','ytick','')
  idx_k=find(diff(T.Label_lobe(idx)));
 idx_k(1)=[];
 for k=1:numel(idx_k)
  line([idx_k(k);idx_k(k)]+0.5,[0.5;246.5],'color','k')
  line([0.5;246.5],[idx_k(k);idx_k(k)]+0.5,'color','k')
 end
 colormap(ha1,bipolar);
 
 ha2=axespos(ax2,j);
 imagesc(T.Label_lobe(idx))
 axis on; set(gca,'xtick','','ytick','');
 colormap(ha2, sgcolormap('hsv-soft'));
 
 axespos(ax3,j+J)
 hold on;
 x=triu(WCC_unc(idx,idx,j).*~eye(246));
 [ci,xi]=hist(x(~~x(:)),100);
 bar(xi(xi<0), ci(xi<0), 1,'facecolor',[.4 .4 1],'linestyle','none');
 bar(xi(xi>0), ci(xi>0), 1,'facecolor',[1 .4 .4],'linestyle','none');
 xlim([-1 1]); xlabel(['WCC_unc:',num2str(j)]); ylabel('# of edges');
end
colormap(bipolar)
screen2png([dir_subj,'/BN_vol/ME_WCC_unc.png'],200)
close(gcf)
%%

end
