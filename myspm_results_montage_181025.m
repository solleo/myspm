function JOB = myspm_results_montage(JOB)
% JOB = myspm_results_montage(JOB)
%
% It creates summary montage for a given GLM directory
%
% JOB requires:
% (.dir_glm)      'Nx1' directory to save SPM results (default=pwd)
% (.slicedim)   [1x1] slice dimension (default=3, axial in RAS)
% (.fname_t1w)
% (.pitch)
%
% (cc) 2018, sgKIM. solleo@gmail.com, https://ggooo.wordpress.com/

if nargin < 1,  JOB=[];  end
if ~isfield(JOB,'dir_glm'), JOB.dir_glm=pwd; end

load([JOB.dir_glm,'/SPM.mat']);
% See if the even contrasts are flipped ones:
if sum(SPM.xCon(1).c + SPM.xCon(2).c) == 0
  bothsigns=1;
  K=1:2:numel(SPM.xCon);
else
  bothsigns=0;
  K=1:numel(SPM.xCon);
end
[~,modelname]=myfileparts(JOB.dir_glm);
JOB.fnames_png={};
for i=1:numel(K)
  k=K(i);
  cfg=JOB;
  cfg.SPM=[SPM.xCon(k).STAT,'_',zeropad(k,4)];
  if bothsigns
    cfg.cntrst={SPM.xCon(k).name,SPM.xCon(k+1).name};
  else
    cfg.cntrst={SPM.xCon(k).name,''};
  end
  cfg.title=[modelname,':',SPM.xCon(k).name];
  if isfield(JOB,'title_suffix')
    cfg.title=[cfg.title,' ',JOB.title_suffix];
  end
  if isfield(cfg,'fname_struct') && ~isfield(cfg,'fname_t1w')
    cfg.fname_t1w=cfg.fname_struct;
  end
  JOB.fnames_png{i}=imagemontage_20181025(cfg);
end
end

function [fname,cfg]=imagemontage_20181025(cfg)
% cfg requires:
%  .dir_glm
% (.cntrst)  {1x2} dir_glm/sigclus_"+var".nii
%  .SPM       dir_glm/spm"T_0001".nii
%  .title
% (.mni_xyz)
%  .slicedim
% (.radius)
% (.seedcoord)
% (.fname_seed)
% (.fname_t1w)
% (.baseimagerange)
%
% (cc) 2016, sgKIM

fname='';
if ~nargin, help(mfilename),return;end
[~,pathname]=myfileparts(cfg.dir_glm);
if ~isfield(cfg,'radius'), radius=10; else radius=cfg.radius; end
pwd0=pwd;
cd(cfg.dir_glm);
if ~isfield(cfg,'dir_fig')
  dir_fig=pwd;
else
  dir_fig=cfg.dir_fig;
end
if ~isfield(cfg,'fname_t1w')
  fname_t1w='~/MNI152_T1_1mm.nii';
else
  fname_t1w=cfg.fname_t1w;
end
%cmap=sgcolormap('GRAY-BLUE-RED-LIGHT');
cmap=sgcolormap('GRAY-BLUE-RED-BRIGHT');
if isfield(cfg,'colormap')
  cmap=colormap;
end
t1w=load_uns_nii(fname_t1w);
if strcmp(fname_t1w,'~/MNI152_T1_1mm.nii')
  baseimagerange=[2500 11000];
else
  imrange=[mode(t1w.img(:)) max(t1w.img(:))];
  baseimagerange=[imrange(1)+diff(imrange)*0.1 imrange(2)];
end
if isfield(cfg,'baseimagerange')
  baseimagerange=cfg.baseimagerange;
end


%% Read T(F)-map and sigclus-map(s)
if isfield(cfg,'custommap')
  fname_nii=cfg.custommap;
else
  fname_nii=['spm',cfg.SPM,'.nii'];
end
if ~isfield(cfg,'cntrst')
  cfg.cntrst={'',''};
end
over=load_untouch_nii_like(fname_nii,fname_t1w);
fname1=['sigclus_',cfg.cntrst{1},'.nii'];
fname2=['sigclus_',cfg.cntrst{2},'.nii'];
if exist(fname1,'file')
  hdr=load_untouch_header_only(fname1);
  clus1=load_untouch_nii_like(fname1, fname_t1w); % interpolation for visualization
else
  clus1=over;
  clus1.img=clus1.img*0;
end
if exist(fname2,'file')
  hdr=load_untouch_header_only(fname2);
  clus2=load_untouch_nii_like(fname2, fname_t1w); % interpolation for visualization
else
  clus2=over;
  clus2.img=clus2.img*0;
end
if exist(fname1,'file') || exist(fname2,'file')
  thres_u = str2double(hdr.hist.descrip(22:26));
  thres_k = str2double(hdr.hist.descrip(33:end));
  if isfield(cfg,'caxis')
    caxis1=cfg.caxis;
  else
    caxis1=[thres_u max(abs(over.img(:)))];
  end
  clus=clus1.img>0.5;
  clus(clus2.img>0.5) = 1;
  clus=double(clus);
  T_thres = over.img .* ~~clus;
else
  warning('No significant clusters/voxels found. Nothing to display.')
  return
%   T_thres = over.img;
%   caxis1=cfg.caxis;
end

%% Read a table
tab = myspm_readtable;
if ~isempty(tab)
  idx=ismember(tab.contrast_name, cfg.cntrst{1});
  mni_xyz=[tab.x_mm(idx) tab.y_mm(idx) tab.z_mm(idx)];
  idx=ismember(tab.contrast_name, cfg.cntrst{2});
  mni_xyz=[mni_xyz; tab.x_mm(idx) tab.y_mm(idx) tab.z_mm(idx)];
end
if isfield(cfg,'mni_xyz'), mni_xyz=cfg.mni_xyz; end
mni_ijk=round(xyz2ijk(mni_xyz,over));
if ~isfield(cfg,'slicedim')
  slicedim=ones(1,size(mni_ijk,1))*3;
else
  slicedim=cfg.slicedim;
end
tab_mni_ijk=mni_ijk(1,:);
if ~isfield(cfg,'seedcoord')
  seedcoord=mni_xyz(1,:);
else
  seedcoord=cfg.seedcoord;
end
seedcoord_ijk = xyz2ijk(seedcoord, over);
num_clus=size(mni_xyz,1);
%% seed ROI
if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
  roi = load_untouch_nii_like(cfg.fname_seed, fname_t1w);
  roi.img = round(roi.img);
elseif isfield(cfg,'roi')
  roi = cfg.roi;
end
if numel(slicedim) == 1
  slicedim=repmat(slicedim,[1 num_clus]);
end
%%
n=size(mni_xyz,1);
figpos=[-1721 127 126.8929+129.6429*n 201];
if isfield(cfg,'width_b')
  figpos(3)=cfg.width_b;
end
figure('color','k', 'position',figpos);
base=double(t1w.img);
%% rotation images?
if isfield(cfg,'pitch')&&cfg.pitch
  base = imrotate3(base,cfg.pitch,[0 1 0],'cubic','crop','FillValues',0);
  T_thres = imrotate3(T_thres,cfg.pitch,[0 1 0],'cubic','crop','FillValues',0);
end

%% main slices
ax=axeslayout1(num_clus+1,[1,num_clus+1]);
% if ~isfield(cfg,'bbox'), cfg.bbox=[-78 -112 -70; 78 76 85]; end
% bbox_ijk=round(xyz2ijk(cfg.bbox,over));
for i=1:num_clus
  if isfield(cfg,'subcortzoom') && cfg.subcortzoom(i)
    ax3=ax;
    ax3.x=ax3.x+0.001;
    ax3.w=ax3.w-0.01;
  else
    ax3=ax;
  end
  switch slicedim(i)
    case 3
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['z = ',num2str(mni_xyz(i,3)),' mm'];
    case 2
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['y = ',num2str(mni_xyz(i,2)),' mm'];
    case 1
      cfg1 = struct('ijk',mni_ijk(i,:));
      label1=['x = ',num2str(mni_xyz(i,1)),' mm'];
  end
  j=i+1;
  ha=axespos(ax3,j);
  if isfield(cfg,'seedcoord')
    cfg1.crosshair=seedcoord_ijk;
    cfg1.sphere.center_vox = seedcoord_ijk;
    cfg1.sphere.radius_vox = radius/mean(t1w.hdr.dime.pixdim(2:4));
  end
  cfg1.xdir='rev';
  cfg1.baseimagerange=baseimagerange;
  cfg1.caxis=caxis1;
  cfg1.slicedim = slicedim(i);
  if i == 1
    cfg1.isR=1;
  end
  if isfield(cfg,'isR')
    cfg1.isR=cfg.isR(i);
  end
  if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
    cfg1.roi=roi.img; % seed mask contour
  elseif isfield(cfg,'roi')
    cfg1.roi=double(roi);
  end
  if isfield(cfg,'baseedge')
    cfg1.baseedge = cfg.baseedge;
  end
  if isfield(cfg,'subcort')
    cfg1.subcort=cfg.subcort(i);
  end
  if isfield(cfg,'subcortbox')
    cfg1.subcortbox=cfg.subcortbox(i);
  end
  if isfield(cfg,'subcortzoom')
    cfg1.subcortzoom=cfg.subcortzoom(i);
  end
  if isfield(cfg,'subcortannot')
    cfg1.subcortannot = cfg.subcortannot(i);
  end
  if isfield(cfg,'roicontourcolor')
    cfg1.roicontourcolor = cfg.roicontourcolor;
  end
  cfg1.colormap=cmap;
  
  % execute imageover1
  imageover1(base, T_thres, cfg1);
  if isfield(cfg,'axis')
    axis(cfg.axis) %% Zoom
  end
  if isfield(cfg,'kill')&&cfg.kill(i)
    cla
  end
  
  % adjust axes size: coronal section looks TOO low compared to sagittal section
  switch slicedim(i)
    case 1
      set(ha,'position',[ax3.x(j) 0.05 ax3.w ax3.h])
    case 2
      set(ha,'position',[ax3.x(j) 0.28 ax3.w ax3.h*.6])
    case 3
      set(ha,'position',[ax3.x(j)+ax3.w*0.05 ax3.y(j) ax3.w*.9 ax3.h])
  end
  
  % coordinate label
  coordinate_label_y=0.85;
  coordinate_label_fontsize=7;
  hal=axes('position',[ax3.x(j) coordinate_label_y ax3.w 0.1]); axis off;
  text(hal,0.5,0,label1,'fontsize',coordinate_label_fontsize, ...
    'backgroundcolor','none', 'color','w', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
end

%% Saggittal navigator
ax=axeslayout1(num_clus+1,[1,num_clus+1]);
ax3=ax;
ax3.y(1)=0.05;
axespos(ax3,1)
if ~isfield(cfg,'sagnav_xyz')
  cfg1 = struct('ijk',tab_mni_ijk, 'slicedim',1, 'caxis',caxis1);
else
  cfg1 = struct('ijk',round(xyz2ijk(cfg.sagnav_xyz,t1w)), ...
    'slicedim',1, 'caxis',caxis1);
end
N=size(mni_ijk,1);
if isfield(cfg,'lines_y')
  cfg1.lines_y = cfg.lines_y;
else
  cfg1.lines_y = slicedim==2;
end

for j=1:N
  if slicedim(j)==2 && cfg1.lines_y(j)
    cfg1.slices(j) = mni_ijk(j,2);
  elseif slicedim(j)==3
    cfg1.slices(j) = mni_ijk(j,3);
  else
    cfg1.slices(j) = nan;
  end
end
if isfield(cfg,'fname_seed') && ~isfield(cfg,'roi')
  cfg1.roi=roi.img; % seed mask contour
elseif isfield(cfg,'roi')
  cfg1.roi=double(roi);
end
cfg1.baseimagerange = baseimagerange;
if isfield(cfg,'baseedge')
  cfg1.baseedge = cfg.baseedge;
end
if isfield(cfg,'nolines'), cfg1.nolines=cfg.nolines; end
cfg1.colormap=cmap;
imageover1(base, T_thres, cfg1);

% coordinate label
if cfg1.slicedim == 1
  xyz_sag=round(ijk2xyz(cfg1.ijk, t1w));
  label1=['x = ',num2str(xyz_sag(1,1)),' mm'];
  hal=axes('position',[ax3.x(1) coordinate_label_y ax3.w 0.1]); axis off;
  text(hal,0.5,0,label1,'fontsize',coordinate_label_fontsize, ...
    'backgroundcolor','none', 'color','w', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle');
end
if ~isfield(cfg,'notitle')
  if ~isfield(cfg,'title')
    cfg.title=pathname;
  end
  Title=cfg.title;
  if ~isfield(cfg,'titlefontsize'), cfg.titlefontsize=13; end
  % TITLE label
  hal=axes('position',[0 0.94 ax3.w 0.1]); axis off;
  if isfield(cfg,'titlepositiony')
    set(hal,'position',[0 cfg.titlepositiony 1 0.1]);
  end
  text(hal,0.01,0, Title, 'fontsize',cfg.titlefontsize, ...
    'color','w', 'interp','none', 'BackgroundColor', 'none', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment','Middle');
end

%% colorbar
ax2=ax;
ax2.y(1)=0.8;
ax2.h=0.2;
pos1=[ax2.w/1.7857  0.22  ax2.w/3.125  0.0311];
pos2=[ax2.w/8.3333  0.22  ax2.w/3.125  0.0311];

H=axespos(ax2,1); hold on; axis off
set(H,'color','k', 'visible','off')
cb1=colorbar('peer',H,'north','color','w');
set(cb1,'position',pos1, 'color','w')
caxis([caxis1(1) caxis1(2)]);
xticks={num2str(caxis1(1),2),num2str(caxis1(2),2)};
colormap(H,cmap(129:end,:)) % red-yellow
xtick=[caxis1(1)+diff(caxis1)*0.1 caxis1(2)-diff(caxis1)*0.1];
set(cb1, 'xtick', xtick, 'xticklabel', xticks, 'fontsize',10);
H=axespos(ax2,1); hold on; axis off;
cb2=colorbar('peer',H,'north','color','w');
set(cb2,'position',pos2, 'color','w')
caxis2=[-caxis1(2) -caxis1(1)];
caxis(caxis2);
colormap(H,cmap(65:128,:)) % blue-cyan
xtick=[caxis2(1)+diff(caxis2)*0.15 caxis2(2)-diff(caxis2)*0.15];
xticks={num2str(caxis2(2),2),num2str(caxis2(1),2)};
set(cb2, 'xtick',xtick, 'xticklabel',xticks , 'fontsize',10, 'xdir','rev');
if min(over.img(:))>=0
  set(cb2,'visible','off')
  set(cb1,'position',[ax2.w/6  0.22  ax2.w/3.125*2  0.0311]);
end

%% save
if isfield(cfg,'suffix')
  suffix=cfg.suffix;
else
  suffix='';
end
if ~isfield(cfg,'fname_png')
  fname=[dir_fig,'/fig_',pathname,'_',cfg.cntrst{1},suffix,'.png'];
else
  fname=cfg.fname_png;
end
[d1,~,~]=fileparts(fname);
[~,~]=mkdir(d1);
if ~isfield(cfg,'dpi')
  dpi=150;
else
  dpi=cfg.dpi;
end
screen2png(fname,dpi);

%% crop the bottom
if ~isfield(cfg,'donotcropit')
  I=imread(fname);
  ymargin=round(200/900*dpi);
  imwrite(I(1:end-ymargin,:,:), fname);
end
cd(pwd0);
end
