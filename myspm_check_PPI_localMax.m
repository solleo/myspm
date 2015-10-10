function myspm_check_PPI_localMax(EXP)
% myspm_check_PPI_localMax(EXP)
%
% EXP
%  .dir_glm
%  .dir_ppi
% checks original coordinates and adjusted local maxima
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

dir0 = EXP.dir_base;
EXP.subjID = fsss_subjID(EXP.subjID);
NumSubj = numel(EXP.subjID);

figure('color','k', 'position',[1921, 1, 1024, 1176]);
[~,ax] = getLayout(NumSubj+1);
for i = 1 : NumSubj
  subjid = EXP.subjID{i};  
  fname_strc = fullfile(dir0,subjid,EXP.name_strc);
  fname_func = fullfile(dir0,subjid,EXP.dir_glm,EXP.name_func);
  hdr_strc = load_nii_hdr(fname_strc);
  hdr_func = load_nii_hdr(fname_func);
  if sum(hdr_strc.dime.dim(2:4) == hdr_func.dime.dim(2:4)) ~= 3
    [~,b,c] = fileparts(EXP.name_func);
    fname_strc_in_func = [dir0,'/',subjid,'/',b,'_in_func',c];
    if ~exist(fname_strc_in_func,'file')
      unix(['mri_convert --like ',fname_func,' ',fname_strc,' ',fname_strc_in_func]);
    end
    fname_strc = fname_strc_in_func;
  end
  base = load_nii(fname_strc);
  func = load_nii(fname_func);
  if ~isfield(EXP,'dir_ppi')
    EXP.dir_ppi = EXP.dir_glm;
  end
  load(fullfile(dir0,subjid,EXP.dir_ppi,'PPI.mat'),'PPI');
  
  d = size(base.img);
  coord_xyz =  PPI.xY.xyz';
  coord_ijk =  round(xyz2ijk(coord_xyz,fname_func));
  coord_ijk0 = round(xyz2ijk(EXP.voi.coord,fname_func));
  
  axespos(ax,i); hold on;
  img1 = base.img(:,:,coord_ijk(3))';
  img2 = func.img(:,:,coord_ijk(3))';
  img2(isnan(img2))=0;
  img2(img2<EXP.thres) = 0;
  imagesc(zeroone(img1)+zeroone(img2));
  cmap=sgcolormap('GRAY-RED');
  if isfield(EXP,'isneg')&&EXP.isneg, cmap=sgcolormap('GRAY-BLUE'); end
  colormap(cmap);
  caxis([0 2])
  line([d(1)-coord_ijk0(1);d(1)-coord_ijk0(1)],[0;d(2)+1],'color','w');
  line([0;d(1)+1], [coord_ijk0(2);coord_ijk0(2)],'color','w');
  if sum(coord_ijk0==coord_ijk)<3,
    vx=d(1)-coord_ijk(1); vy=coord_ijk(2); vr=EXP.voi.radius;
    rectangle('position',[vx-vr, vy-vr, vr*2, vr*2],'curvature',[1 1],'edgecolor','k');
  end
  axis([1 d(1) 1 d(2)]); axis off; axis image; set(gca,'ydir','nor');
  if isfield(EXP,'masksubjid')&&EXP.masksubjid
    subjname=['subj',pad(i,3)];
  else
    subjname=subjid;
  end
  text(23,10,[subjname,':T>',num2str(round(EXP.thres*10)/10)], 'backgroundColor','w')
  
  drawnow;
end
axespos(ax,NumSubj+1);
set(gca,'color','k'); axis off;
text(0,0.3,EXP.voi.name,  'interp','none','color','w');
text(0,0.2,EXP.dir_glm,   'interp','none','color','w');
text(0,0.1,EXP.name_func, 'interp','none','color','w');
screen2png([EXP.dir_fig,'/',EXP.voi.name,'_',EXP.dir_glm,'_',EXP.name_func,'.png']);

 close(gcf)
end