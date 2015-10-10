function myspm_check_localMax(EXP)
% myspm_check_localMax(EXP)
%
% checks original coordinates and adjusted local maxima
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https://ggooo.wordpress.com

dir0 = EXP.dir_base;
EXP.subjID = fsss_subjID(EXP.subjID);
NumSubj = numel(EXP.subjID);

figure;
[~,ax] = getLayout(NumSubj);
for i = 1 : NumSubj
  subjid = EXP.subjID{i};

  fname_strc = [dir0,'/',subjid,'/',EXP.name_strc];
  fname_func = [dir0,'/',subjid,'/',EXP.name_func];
  
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
  
  load([dir0,subjid,'/glm_FC+FD+BC+BD_RP/PPI_',voi(k).name,'xFCBD.mat']);
  d = size(base.img);
  coord_xyz = PPI.xY.xyz';
  coord_ijk = round(xyz2ijk(coord_xyz,fname4));
  
  axespos(ax,i); hold on;
  img1 = base.img(:,:,coord_ijk(3))';
  img2 = func.img(:,:,coord_ijk(3))';
  img2(isnan(img2))=0;
  imagesc(img2);
  scatter(d(1)-coord_ijk0(1),coord_ijk0(2),100,'kx','linewidth',2); % initial coord
  scatter(d(1)-coord_ijk(1), coord_ijk(2),100,'wo','linewidth',2);  % nearest local maximum
  axis image; set(gca,'ydir','nor');%,'xdir','rev');
  
  
end
