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
  coord_ijk0 = round(xyz2ijk(EXP.voi_coord,fname_func));
  
  axespos(ax,i); hold on;
  img1 = base.img(:,:,coord_ijk(3))';
  img2 = func.img(:,:,coord_ijk(3))';
  img2(isnan(img2))=0;
  img2(img2<EXP.thres) = 0;
  imagesc(zeroone(img1)+zeroone(img2));
  cmap=hot(64*2);
  colormap([bone(64);cmap(70-20:129-16,:)]);
  caxis([0 2])
  scatter(d(1)-coord_ijk0(1),coord_ijk0(2),200,'b+','linewidth',2); % initial coord
  scatter(d(1)-coord_ijk(1), coord_ijk(2),100,'co','linewidth',2);  % nearest local maximum
  axis off; axis image; set(gca,'ydir','nor');
  text(23,10,[subjid,':T>',num2str(round(EXP.thres*10)/10)], 'backgroundColor','w')
  drawnow;
end
axespos(ax,NumSubj+1);
set(gca,'color','k'); axis off;
text(0,0.3,EXP.voi_name,  'interp','none','color','w');
text(0,0.2,EXP.dir_ppi,'interp','none','color','w');
text(0,0.1,EXP.name_func, 'interp','none','color','w');
screen2png([EXP.dir_fig,'/',EXP.voi_name,'_',EXP.dir_ppi,'_',EXP.name_func,'.png']);
close(gcf)
end