function myspm_wds(dir_subj, figs, TR_sec)
% myspm_wds(dir_subj, figs)
%
% (dir_subj) default= pwd
% (figs)     default= 1
if ~exist('figs','var'), figs=1; end

%{
Let's try Patel's way first as described in
Patel & Bullmore, 2016, "A wavelet-based estimator of the degrees of
freedom in denoised fMRI time series for probabilistic testing of
functional connectivity and brain graphs", NeuroImage.

Core Image Processing included the following steps:
(i) slice ac-quisition correction;
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

%% FOR ME-ICA data, nuissance regression (viii) was omitted.

%% FOR the parcellation in the native space, 
%  sptial normalization & smoothing are omitted.

%%
if ~exist('dir_subj','var')
 dir_subj=pwd;
end
cd(dir_subj)
epi_suffix='me';
if ~exist('TR_sec','var')
 hdr=load_untouch_header_only(['spm12/',epi_suffix,'.nii']);
 TR_sec=hdr.dime.pixdim(5);
end
% 1. (vii) WDS
if ~exist('BN_vol/rWD.mat','file') || 1
 load ('BN_vol/ME.mat','ME')
 [NT,NV]=size(ME);
 vol=zeros(NV,1,1,NT);
 for t=1:NT
  vol(:,1,1,t)=ME(t,:);
 end
 nii=make_nii(vol);
 save_nii(nii,'BN_vol/ME.nii')
 WaveletDespike('BN_vol/ME.nii','BN_vol/ME')
 nii=load_nii('BN_vol/ME_wds.nii.gz');
 WD=squeeze(nii.img)';
 SP=load('BN_vol/ME_SP.txt');
 
 % 2. Motion descriptors
 fname_rp=['BN_vol/motion.1D'];
 rp=load(fname_rp);
 % 2-1. frame-wise displacement
 arc_length_change= 50^2 * diff([0 0 0; rp(:,4:6)]) ./ (360);
 FDp = sum(abs( [ [0 0 0; diff(rp(:,1:3))], arc_length_change ] ),2);
 %FDv = rms( [0 0 0; diff(rp(:,1:3))],2 );
 FDs = l2norm( [ [0 0 0; diff(rp(:,1:3))], [0 0 0; diff(rp(:,4:6))] ] );
 
 % 2-3. DVARS l2norm of percent change
 PCme=[zeros(1,NV); diff(ME)./repmat(ME(1,:),(NT-1),1)]*100;
 DVARSme = rms(PCme,2);
 
%  CC=load(['spm12/ua',epi_suffix,'_n3b0.00-Inf_eigenvec.txt']);
%  DVARScc = l2norm([zeros(1,size(CC,2)); diff(CC)./repmat(CC(1,:),(NT-1),1)])*100;
 %%
 figure
 subplot(411)
 plot(FDp); grid minor; ylabel('FD-power [mm]');xlim([1 NT])
 subplot(412)
 plot(DVARSme); grid minor; ylabel('DVARS-ME [%]');xlim([1 NT])
 subplot(413);
 imagesc(zscore(ME)'); caxis([-3 3]);xlim([1 NT]); ylabel('Nodes'); title('Z(ME)')
 subplot(414)
 imagesc(zscore(WD)'); caxis([-3 3]);xlim([1 NT]); ylabel('Nodes'); title('Z(WD)')
 colormap(sgcolormap('paula'))
 
 %%
 WCCme = wavcorr(ME, 'warning',0, 'mincoef',100);
 
 cfg=[];
 cfg.colorbartitle=['PCC_WCC-ME:scale3'];
 cfg.fname_png=['PCC_WCC-ME_j3.png'];
 cfg.caxis=[-1 1];
 fsss_view_BN(squeeze(WCCme(151,:,3)+WCCme(152,:,3))/2, cfg);
 
 
 WCCwd = wavcorr(WD, 'warning',0, 'mincoef',100);
 cfg=[];
 cfg.colorbartitle=['PCC_WCC-Wd:scale3'];
 cfg.fname_png=['PCC_WCC-WD_j3.png'];
 cfg.caxis=[-1 1];
 fsss_view_BN(squeeze(WCCwd(151,:,3)+WCCwd(152,:,3))/2, cfg);
 
 %%
%  % 2-3. get residual
%  M=[ones(NT,1) rp ];
%  rWD = WD-M*(pinv(M'*M)*M'*WD);
%  save('BN_vol/rWD.mat','rWD')
%  %%
%  if figs
%   hf=figure('position',[40     1   704   562]);
%   subplot(431)
%   imagesc(zscore(M(:,2:end))'); caxis([-3 3])
%   title('Move6 + ||dMdt|| ','fontsize',10)
%   ylabel('Params'); xlabel('TR')
%   subplot(434)
%   imagesc(zscore(ME)'); caxis([-3 3])
%   title('Raw','fontsize',10)
%   ylabel('Nodes');
%   subplot(437)
%   imagesc(zscore(WD)'); caxis([-3 3])
%   title('Raw:WD','fontsize',10)
%   subplot(4,3,10)
%   imagesc(zscore(rWD)'); caxis([-3 3])
%   title('Raw:WD:reg','fontsize',10)
%   
%   subplot(432)
%   plot(l2norm(zscore(M(:,2:end)))); xlim([1 NT])
%   ylabel('L_2norm'); xlabel('TR');
%   title('Move6 + ||dMdt|| ','fontsize',10)
%   subplot(435); hold on;
%   plot(l2norm(zscore(ME)),'r'); xlim([1 NT])
%   ylim([0 max(l2norm(zscore(ME)))])
%   title('Raw','fontsize',10)
%   subplot(438); hold on
%   plot(l2norm(zscore(ME)),'color',[.5 .5 .5]); xlim([1 NT])
%   plot(l2norm(zscore(WD)),'r'); xlim([1 NT])
%   ylim([0 max(l2norm(zscore(ME)))])
%   title('Raw:WD','fontsize',10)
%   subplot(4,3,11); hold on;
%   plot(l2norm(zscore(WD)),'color',[.5 .5 .5]);
%   plot(l2norm(zscore(rWD)),'r'); xlim([1 NT])
%   title('Raw:WD:reg','fontsize',10)
%   
%   ha=axes('position',[0.02 0.73 0.05 0.2]);
%   hcb=colorbar(ha,'location','east');
%   title(hcb,' Z-socre'); caxis([-3 3])
%   axis off
%   
%   ha=axes('position',[0.05 0.13 0.02 0.137]);
%   imagesc(SP); title('SP%','fontsize',9)
%   caxis([0 100]);
%   
%   subplot(436); hold on;
%   SurfStatPlot(l2norm(zscore(M(:,2:end))), l2norm(zscore(ME)))
%   title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end))), l2norm(zscore(ME))))],'fontsize',10)
%   xlabel('||Params||'); ylabel('||Raw||')
%   
%   subplot(439); hold on;
%   SurfStatPlot(l2norm(zscore(M(:,2:end))), l2norm(zscore(WD)))
%   title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end))), l2norm(zscore(WD))))],'fontsize',10)
%   xlabel('||Params||'); ylabel('||Raw:reg||')
%   
%   subplot(4,3,12); hold on;
%   SurfStatPlot(l2norm(zscore(M(:,2:end))), l2norm(zscore(rWD)))
%   title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end))), l2norm(zscore(rWD))))],'fontsize',10)
%   xlabel('||Params||'); ylabel('||Raw:reg:WD||')
%  screen2png('BN_vol/ME_WD_denoise.png',200);
%  close(hf);
%  end
% else
%  load('BN_vol/rWD.mat','rWD')
% end

%% 2. Computing WCC
[WCC,ncoef] = wavcorr(rWD, 'warning',0, 'mincoef',100);
ncoef(~ncoef)=[];
J=numel(ncoef); % Let's see WCC upto scale 5:

%% FIGURE: R matrix
T=readtable('~/Dropbox/BN_Atlas/BN_Atlas_246_COT.xlsx');
band_mHz=zeros(size(WCC,3),2);
idx=[1:2:246 246:-2:2];

if figs
figure('position',[8         560        1234/5*J         473])
ax1=axeslayout1(J*2,[2 J],[.1 0]);
ax1.y=ax1.y-0.05;
ax3=axeslayout1(J*2,[2 J],[.22 0.22]);
for j=1:J
 ha1=axespos(ax1,j);
 hold on
 imagesc(WCC(idx,idx,j)); axis image; caxis([-1 1]);
 band_mHz(j,:)=[2^-(j+1), 2^-j]/TR_sec*1000;
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
 set(gca,'ydir','rev');
 axis([0 246 0 246]); box on
 
 axespos(ax3,j+J)
 hold on;
 x=triu(WCC(idx,idx,j).*~eye(246));
 [ci,xi]=hist(x(~~x(:)),100);
 bar(xi(xi<0), ci(xi<0), 1,'facecolor',[.4 .4 1],'linestyle','none');
 bar(xi(xi>0), ci(xi>0), 1,'facecolor',[1 .4 .4],'linestyle','none');
 xlim([-1 1]); xlabel(['WCC:',num2str(j)]); ylabel('# of edges');
end
screen2png([dir_subj,'/BN_vol/ME_WCC.png'],200)
close(gcf)
end

%%  5-2. Correction of effective DOF:
nii=load_nii('BN_vol/ME_EDOF.nii.gz');
eDOF=squeeze(nii.img)';
Z=WCC*0;
P=WCC*0;
WCC_cor=WCC*0;
Yeta0=ncoef*0;
for j=1:J
 % pair-wise df is defined by minimum of df's of nodes (NI(2016)142:p.20)
 D=zeros(246);
 for a=1:246
  for b=1:246
   D(a,b)=min(eDOF(j,a), eDOF(j,b));
  end
 end
 yeta0=min(D(:));
 Yeta0(j)=yeta0;
 if yeta0>3
  for a=1:245
   for b=(a+1):246
    z = atanh(WCC(a,b,j))*sqrt(D(a,b)-3);
    Z(a,b,j) = z;
    Z(b,a,j) = z;
    if z>0
     p = (1-normcdf(z,0,1))*2;
    else
     p = (1-normcdf(-z,0,1))*2;
    end
    P(a,b,j) = p;
    P(b,a,j) = p;
    WCC_cor(a,b,j) = tanh(z/sqrt(yeta0-3));
    WCC_cor(b,a,j) = WCC_cor(a,b,j);
   end
  end
 else
  J=J-1;
  continue
 end
end

%% FIGURE: Z
if figs
figure('position',[8         560        1234/5*J         473])
ax1=axeslayout1(J*2,[2 J],[.1 0]);
ax1.y=ax1.y-0.05;
ax3=axeslayout1(J*2,[2 J],[.22 0.22]);
for j=1:J
 ha1=axespos(ax1,j);
 hold on
 imagesc(Z(idx,idx,j)); axis image; 
 CAXIS=prctile(reshape(abs(Z(idx,idx,j)),[],1),99);
 CAXIS=[-CAXIS CAXIS];
 caxis(CAXIS);
 band_mHz(j,:)=[2^-(j+1), 2^-j]/TR_sec*1000;
 title({['Scale ',num2str(j)],[num2str(round(band_mHz(j,1))),'-',...
  num2str(round(band_mHz(j,2))),' mHz']});
 set(gca,'xtick','','ytick','')
 idx_k=find(diff(T.Label_lobe(idx)));
 idx_k(1)=[];
 for k=1:numel(idx_k)
  line([idx_k(k);idx_k(k)]+0.5,[0.5;246.5],'color','k')
  line([0.5;246.5],[idx_k(k);idx_k(k)]+0.5,'color','k')
 end
 colormap(ha1,sgcolormap('BCWYR256'));
 set(gca,'ydir','rev');
 axis([0 246 0 246]); box on
 
 axespos(ax3,j+J)
 hold on;
 x=triu(Z(idx,idx,j).*~eye(246));
 [ci,xi]=hist(x(~~x(:)),100);
 bar(xi(xi<0), ci(xi<0), 1,'facecolor',[.4 .8 1],'linestyle','none');
 bar(xi(xi>0), ci(xi>0), 1,'facecolor',[1 .8 .4],'linestyle','none');
 xlim(CAXIS); xlabel(['WCC-Z:',num2str(j)]); ylabel('# of edges');
end
screen2png([dir_subj,'/BN_vol/ME_WCC-Z.png'],200)
close(gcf)
end
%% Save all
save ([dir_subj,'/BN_vol/WCC.mat'],'WCC','ncoef','band_mHz','Yeta0','Z','P','WCC_cor');

end
