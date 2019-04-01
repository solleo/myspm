function myspm_wds_sg(dir_subj, figs, TR_sec)
% myspm_wds(dir_subj, figs)
%
% (dir_subj) default= pwd
% (figs)     default= 1
if ~exist('figs','var'), figs=1; end

%{
Now try my way:
(1) 8 motion parameters, 3 aCompCor regressors
(2) then WDS
%}

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
if ~exist('BN_vol/WDr.mat','file')
 load ('BN_vol/ME.mat','ME')
 % 1. regression
 % 1-1: motion parameters (6+2)
 fname_rp=['spm12/rp_a',epi_suffix,'8.txt'];
 rp=load(fname_rp);
 % 1-2. aCompCor (3)
 cc=load(['spm12/ua',epi_suffix,'_n3b0.00-Inf_eigenvec.txt']);
 % 1-3. get residual
 M=[ones(size(rp,1),1) rp cc];
 rME = ME-M*(pinv(M'*M)*M'*ME);
 save('BN_vol/rME.mat','rME')
 
 % 2. WDS
 [NT,NV]=size(rME);
 vol=zeros(NV,1,1,NT);
 for t=1:NT
  vol(:,1,1,t)=rME(t,:);
 end
 nii=make_nii(vol);
 save_nii(nii,'BN_vol/rME.nii')
 WaveletDespike('BN_vol/rME.nii','BN_vol/rME')
 nii=load_nii('BN_vol/rME_wds.nii.gz');
 WDr=squeeze(nii.img)';
 SP=load('BN_vol/rME_SP.txt');
 
 % 3. just BPF?
 
%  
%  
%  fWDr = myy_filter(WDr, TR_sec, [0.009, 0.080]);
%  imagesc(corr(fWDr)); axis image; caxis([-1 +1])
 %%
 if figs
  hf=figure('position',[40     1   704   562]);
  subplot(431)
  imagesc(zscore(M(:,2:end))'); caxis([-3 3])
  title('Move6 + ||dMdt|| + CC3','fontsize',10)
  ylabel('Params'); xlabel('TR')
  subplot(434)
  imagesc(zscore(ME)'); caxis([-3 3])
  title('Raw','fontsize',10)
  ylabel('Nodes');
  subplot(437)
  imagesc(zscore(rME)'); caxis([-3 3])
  title('Raw:reg','fontsize',10)
  subplot(4,3,10)
  imagesc(zscore(WDr)'); caxis([-3 3])
  title('Raw:reg:WD','fontsize',10)
  
  subplot(432)
  plot(l2norm(zscore(M(:,2:end)))); xlim([1 NT])
  ylabel('L_2norm'); xlabel('TR');
  title('Move6 + ||dMdt|| + CC3','fontsize',10)
  subplot(435); hold on;
  plot(l2norm(zscore(ME)),'r'); xlim([1 NT])
  ylim([0 max(l2norm(zscore(ME)))])
  title('Raw','fontsize',10)
  subplot(438); hold on
  plot(l2norm(zscore(ME)),'color',[.5 .5 .5]); xlim([1 NT])
  plot(l2norm(zscore(rME)),'r'); xlim([1 NT])
  ylim([0 max(l2norm(zscore(ME)))])
  title('Raw:reg','fontsize',10)
  subplot(4,3,11); hold on;
  plot(l2norm(zscore(rME)),'color',[.5 .5 .5]);
  plot(l2norm(zscore(WDr)),'r'); xlim([1 NT])
  title('Raw:reg:WD','fontsize',10)
  
  ha=axes('position',[0.02 0.73 0.05 0.2]);
  hcb=colorbar(ha,'location','east');
  title(hcb,' Z-socre'); caxis([-3 3])
  axis off
  
  ha=axes('position',[0.05 0.13 0.02 0.137]);
  imagesc(SP); title('SP%','fontsize',9)
  caxis([0 100]);
  
  subplot(436); hold on;
  SurfStatPlot(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(ME)))
  title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(ME))))],'fontsize',10)
  xlabel('||Params||'); ylabel('||Raw||')
  
  subplot(439); hold on;
  SurfStatPlot(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(rME)))
  title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(rME))))],'fontsize',10)
  xlabel('||Params||'); ylabel('||Raw:reg||')
  
  subplot(4,3,12); hold on;
  SurfStatPlot(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(WDr)))
  title(['r = ',num2str(corr(l2norm(zscore(M(:,2:end-3))), l2norm(zscore(WDr))))],'fontsize',10)
  xlabel('||Params||'); ylabel('||Raw:reg:WD||')
  
  screen2png('BN_vol/ME_WD_denoise_sg.png',200);
  close(hf);
 end
 save('BN_vol/WDr.mat','WDr')
else
 load('BN_vol/WDr.mat','WDr')
end

%% 2. Computing WCC
[WCC,ncoef] = wavcorr(WDr, 'warning',0, 'mincoef',100);
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
 screen2png([dir_subj,'/BN_vol/ME_WCC_sg.png'],200)
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
 screen2png([dir_subj,'/BN_vol/ME_WCC_sg-Z.png'],200)
 close(gcf)
end
%% Save all
save ([dir_subj,'/BN_vol/WCC_sg.mat'],'WCC','ncoef','band_mHz','Yeta0','Z','P','WCC_cor');

%%
if figs
 for j=1:3
  cfg=[];
  cfg.colorbartitle=['PCC_WCC',num2str(j)];
  cfg.fname_png=['BN_vol/WCC',num2str(j),'_PCC_cor.png'];
  cfg.caxis=[-.5 .5];
  cfg.fstemplate='fsaverage5';
  fsss_view_BN(squeeze(WCC_cor(151,:,j)+WCC_cor(152,:,j))/2, cfg);
  
  cfg=[];
  cfg.colorbartitle=['PCC_WCC',num2str(j)];
  cfg.fname_png=['BN_vol/WCC',num2str(j),'_PCC_unc.png'];
  cfg.caxis=[-.5 .5];
  cfg.fstemplate='fsaverage5';
  fsss_view_BN(squeeze(WCC(151,:,j)+WCC(152,:,j))/2, cfg);
 end
end

end
