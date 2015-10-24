function EXP = myspm_rsfc_seed_vox2vox (EXP)
% EXP = myspm_rsfc_seed_vox2vox (EXP)
%
% computes correlation matrix between voxels in region A vs. voxels in
% region B
%
% EXP requires:
%  .SEEDPAIRS;
%  .dir_base
%
% (cc) 2015, sgKIM   mailto://solleo@gmail.com
global overwrite
if isempty(overwrite), overwrite=0; end
if ~isfield(EXP,'aparcname'), EXP.aparcname='aparc.a2009s+aseg'; end
dir0=pwd;
if ~isfield(EXP,'fwhm'), EXP.fwhm=0; end
fwhm = sprintf('%0.1f',EXP.fwhm);
%%
EXP=[];
EXP.dir_base='/scr/vatikan3/APConn/rest12.410/';
EXP.aparcname='aparc+aseg';
side={'lh-','rh-'};
EXP.SEEDPAIRS={  {'transversetemporal','superiortemporal'} };
EXP.subjID = 2001;
EXP.fwhm=4.6;
fwhm = sprintf('%0.1f',EXP.fwhm);

subjID = fsss_subjID(EXP.subjID);
for j=1:numel(EXP.SEEDPAIRS)
   %seedpairs = EXP.SEEDPAIRS{j};
%   dir_mni=[EXP.dir_base,'/mni152/s',fwhm,'_native_mean_',seedname];
%   [~,~]=mkdir(dir_mni);
  %%
  for n=1:numel(subjID)
    subjid=subjID{n};
%     if ~exist([dir_mni,'/wcbeta2_',subjid,'.nii'],'file')||overwrite
      dir1=[EXP.dir_base,'/',subjid];
      cd (dir1);
      
      % smoothing
      exp1=[];
      exp1.fnames=[dir1,'/fruarest410.nii'];
      exp1.fwhm = EXP.fwhm;
      sfname = [dir1,'/s',fwhm,'fruarest410.nii'];
      if ~exist(sfname,'file');
        myspm_smooth(exp1);
      end
      P = spm_vol(sfname);
      y = load_untouch_nii(sfname,1); % first frame
      
      hf=figure('position',[1064         605         730         535]);
      for s=1:2
        seedpairs = { [side{s},EXP.SEEDPAIRS{j}{1}], [side{s},EXP.SEEDPAIRS{j}{2}]};
      % Find seed coordinates (IJK)
      aparc = load_untouch_nii([dir1,'/o',EXP.aparcname,'.nii']);
      [~,~,roiname,roi_label] ...
        = load_roiname([dir1,'/o',EXP.aparcname,'1.txt']);
      
      for a=1:2
        label{a} = roi_label(ismember(roiname, seedpairs{a}));
        idx{a} = ismember(aparc.img, label{a});
        IJK{a} = find3(idx{a});
        Y{a} = spm_get_data(P, IJK{a}');
      end
      x = pca([IJK{1}; IJK{2}],1); % it should be from ventral to dorsal
      
      paramaxis{1} = x(1:size(IJK{1},1));
      paramaxis{2} = x(size(IJK{1},1)+1:end);
      C = corr(Y{1},Y{2});
      
      subplot(2,1,s)
      hold on;
      vecx=[]; vecy=[];
      for c=1:size(Y{1},2)
        %scatter(paramaxis{2}-paramaxis{1}(c),C(c,:),'k.')
        vecx=[vecx; paramaxis{2}-paramaxis{1}(c)];
        vecy=[vecy;C(c,:)'];
      end
      jhist(vecx,vecy); colormap(flipud(gray)); ylim([-1 1])
      xlabel('dist between voxels along the PC#1 of STG coords (mm)'); ylabel('corr');
      title([subjid,':',seedpairs{1},' vs. ',seedpairs{2}])
      end
      screen2png(['/scr/vatikan3/APConn/figures.local/PAC_vs_STG_corr/corr_over_dist_',subjid,'.png']);
      close(hf);
      
      % and.. and classify????
  end
  
  
  %%
end
            
      

  end


      