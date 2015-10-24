function EXP = myspm_rsfc_seed (EXP)
% EXP = myspm_rsfc_seed (EXP)
%
% EXP requires:
%  .parc.names    {Jx1}
%  .parc.radius   [1x1]
%  .aparcname     'Nx1'   'aparc.a2009s+aseg' (default) or else
% or
%  .sphere.xyz    [Jx3]
%  .sphere.names  {Jx1}
%  .sphere.radius [1x1]
% or
%  .mask.names    {Jx1}
%  .mask.query    'Nx1'   string to find filename, ${subj} will be replaced by subject's id
%
%  .dir_base
% (.ismni)        '1x1'   space to compute rsfc: 0=native (default) or 1=mni152
% (.ispc1)        '1x1'   seed timeseries: 0=mean (default) or 1=first principle component
% (.isglm)        '1x1'   way to compute rsfc: 0=corr (default) or 1=glm
%
% (cc) 2015, sgKIM.

global overwrite
if isempty(overwrite), overwrite=0; end
if ~isfield(EXP,'aparcname'), EXP.aparcname='aparc.a2009s+aseg'; end
dir0=pwd;
if ~isfield(EXP,'fwhm'), EXP.fwhm=0; end
fwhm = sprintf('%0.1f',EXP.fwhm);
if ~isfield(EXP,'ismni152'), EXP.ismni152=0; end;
spaces = {'native','mni152'};
space = spaces{EXP.ismni152+1};
if ~isfield(EXP,'ispc1'), EXP.ispc1=0; end;
seedYs = {'mean','pc1'};
seedY  = seedYs{EXP.ispc1+1};
if ~isfield(EXP,'isglm'), EXP.isglm=0; end;
computes = {'corr','glm'};
compute = computes{EXP.isglm+1};
subjID = fsss_subjID(EXP.subjID);
if isfield(EXP,'parc')
  SEEDNAMES = EXP.parc.names;
elseif isfield(EXP,'sphere')
  SEEDNAMES = EXP.sphere.names;
elseif isfield(EXP,'mask')
  SEEDNAMES = EXP.mask.names;
end
if ~iscell(SEEDNAMES), SEEDNAMES={SEEDNAMES}; end

for j=1:numel(SEEDNAMES)
  seedname = SEEDNAMES{j};
  if isfield(EXP,'mask')&&isfield(EXP.mask,'radius')
    radius = EXP.mask.radius;
  elseif isfield(EXP,'parc')&&isfield(EXP.parc,'radius')
    radius = EXP.parc.radius;
  end  
  if exist('EXP.mask.radius','var')
    seedname_out=[seedname,'_r',num2str(radius),'mm'];
  else
    seedname_out = seedname;
  end
  dir_mni=[EXP.dir_base,'/mni152/s',fwhm,'_',space,'_',seedY,'_',compute,'_',seedname];
  [~,~]=mkdir(dir_mni);  
  
  for n=1:numel(subjID)
    subjid=subjID{n};
    dir2=[EXP.dir_base,'/',subjid,'/s',fwhm,'_',space,'_',seedY,'_',compute,'_',seedname];
    [~,~]=mkdir(dir2);  
    if ~exist([dir_mni,'/wcbeta2_',subjid,'.nii'],'file') || ...
        ~exist([dir_mni,'/wccorr_',subjid,'.nii'],'file') || ...
        overwrite
      dir1=[EXP.dir_base,'/',subjid];
      cd (dir1);
      %% normalization for mni152 space
      if strcmp(space,'mni152') && ~exist([dir1,'wcfruarest410.nii'],'file')
        % coreg: func -> strc
        exp3=EXP;
        exp3.subjID={subjid};
        exp3.prefix='c';
        exp3.name_fixed  =[dir1,'/t1w.nii'];
        exp3.name_moving =[dir1,'/meanuarest410.nii'];
        exp3.name_others ={[dir1,'/fruarest410.nii']};
        myspm_coreg(exp3)
        
        % deform: strc -> mni
        exp4=[];
        exp4.fname_moving=[dir1,'/cfruarest410.nii'];
        exp4.fname_def   =[dir1,'/y_t1w.nii'];
        exp4.vox_mm = 2;
        myspm_deform(exp4)
        
        % deform: strc -> mni
        exp4=[];
        exp4.fname_moving=[dir1,'/aparc.a2009s+aseg.nii'];
        exp4.fname_def   =[dir1,'/y_t1w.nii'];
        exp4.vox_mm = 2;
        exp4.interp = 0;
        myspm_deform(exp4)
      end
      
      %% smoothing
      if strcmp(space,'mni152')
        fname='wcfruarest410.nii';
      else
        %fname='mfruarest410.nii';
        fname='fruarest410.nii';
        if isfield(EXP,'name_epi')
          fname=EXP.name_epi;
        end
      end
      
      exp1=[];
      exp1.fnames=[dir1,'/',fname];
      exp1.fwhm = EXP.fwhm;
      sfname = [dir1,'/s',fwhm,fname];
      if ~exist(sfname,'file');
        myspm_smooth(exp1);
      end
      y = load_untouch_nii(sfname,1);
      
      %% Find seed coordinates (IJK)
      if isfield(EXP,'mask')
        % I need to register mask in 7T-t1w to 3T-fmri space
        query = EXP.mask.fnames_query;
        idx1 = findstr(query,'${');
        idx2 = findstr(query,'}');
        fname=[query(1:idx1-1),subjid,query(idx2+1:end)];
        mask = load_untouch_nii(fname);
        idx1 = ~~mask.img;
        IJK  = find3(idx1);
        if isfield(EXP.mask,'radius')
          radius = EXP.mask.radius;
          ijk0 = mean(IJK); % center of mass of mask
          idx1 = reshape(l2norm(find3(y.img*0+1)...
            -repmat(ijk0,[numel(y.img),1])) ...
            < radius,size(y.img));
          IJK = find3(idx1);
        end
      elseif isfield(EXP,'sphere') % if coordinates are given
        ijk0  = round(xyz2ijk(EXP.sphere.xyz{j},y));
        idx1 = reshape(l2norm(find3(y.img*0+1)-repmat(ijk0,[numel(y.img),1])) ...
          < EXP.sphere.radius,size(y.img));
        IJK = find3(idx1);
      elseif isfield(EXP,'parc')
        if strcmp(space,'mni152'),
          aparc_prefix='w';
          exp0 = EXP;
          exp0.subjID = {subjid};
          myspm_aparc(exp0);
        else
          aparc_prefix='o';
        end
        aparc = load_untouch_nii([dir1,'/',aparc_prefix,EXP.aparcname,'.nii']);
        [~,~,roiname,roi_label] ...
          = load_roiname([dir1,'/',aparc_prefix,EXP.aparcname,'_reidx.txt']);
        label1 = roi_label(ismember(roiname, seedname));
        idx1 = ismember(aparc.img, label1);
        IJK = find3(idx1);
        if isfield(EXP.parc,'radius')
          radius = EXP.parc.radius;
          ijk0 = mean(IJK); % center of mass of parcellation
          idx1 = reshape(l2norm(find3(y.img*0+1)...
            -repmat(ijk0,[numel(y.img),1])) ...
            < radius,size(y.img));
          IJK = find3(idx1);
        end
      end
      
      %% Compute PC#1 and mean
      P = spm_vol(sfname);
      Y = spm_get_data(P, IJK');
      [U,S,V] = svd(Y,'econ');
      meanY = mean(Y,2);
      
      % adjust signs to be positive with mean
      maxU = min([5,size(U,2)]);
      for u1=1:maxU
        if corr(U(:,u1),meanY)<0
          U(:,u1) = -U(:,u1);
        end
      end
      
      %% plot: where in the seed?
      brain = load_untouch_nii([dir1,'/oBrain.nii']);
      K = linspace(min(IJK(:,3)), max(IJK(:,3)), 5);
      J = fliplr(linspace(min(IJK(:,2)), max(IJK(:,2)), 5));
      I = linspace(min(IJK(:,1)), max(IJK(:,1)), 5);
      
      hf=figure('position',[1921 1 1024 1176], 'color','k');
      [~,ax]=axeslayout(5*2*6,[5*2,6],[0 0]);
      [~,ax2]=axeslayout(30,[5,6],[0.1 0.1]);
      
      vals=[U(:,1:4), meanY];
      % PC#1 to #4 and meanY
      for u1=1:maxU
        k=1; ii=1;
        axespos(ax2,1+6*(u1-1));
        plot(vals(:,u1), 'g'); xlim([1 size(Y,1)]); xlabel('TR');
        ylim0=ylim;
        if u1<5
          text(1,ylim0(end)*.9, ['PC#',num2str(u1)], 'color','w');
        else
          text(1,ylim0(end)*.9, 'meanY','color','w');
        end
        set(gca,'xcolor','w','color','k','ycolor','w');
        for k=2:6 % slices
          for ii=1:2
            axespos(ax,k+12*(u1-1)+6*(ii-1));
            corrvol = y.img*0;
            corrvol(idx1) = corr(Y,vals(:,u1));
            if ii==1
              cfg=struct('slicedim',3, 'ijk', [0 0 round(K(k-1))], 'thres', 0);
              h=imageover1(brain.img, corrvol, cfg); set(h,'xdir','rev');
            else
              cfg=struct('slicedim',1, 'ijk', [round(I(k-1)) 0 0], 'thres', 0);
              h=imageover1(brain.img, corrvol, cfg); %set(h,'xdir','rev')
            end
            if k==2 && ii==1
              xlim0=xlim; ylim0=ylim;
              text(xlim0(end), ylim0(end)*0.8, ['max corr=',num2str(max(corrvol(:)))], 'color','w');
            end
          end
        end
      end
      screen2png([dir_mni,'/',subjid,'_pcsNmeans.png']);
      close(hf);
      
      switch compute
        case 'glm'
          %% temporally filtered data, and no psychological condition;
          % > thus 2nd-level glm without HRF modelling
          exp2 = EXP;
          for i=1:size(P,1)
            exp2.fnames{i,1} = [sfname,',',num2str(i)];
          end
          exp2.dir_glm = dir2;
          if strcmp(seedY,'mean')
            exp2.model    = 1 + term(zscore(meanY),seedname);
          else
            exp2.model    = 1 + term(zscore(U(:,1)),seedname);
          end
          exp2.design   = 'mreg';
          exp2.cidx     = 2;
          exp2.noresult = 1;
          exp2 = rmfield(exp2,'fwhm'); % smoothing should be done before running glm in order to extract PC from smoothed data
          myspm_glm(exp2)
          %dir2 = exp2.dir_glm;
          fname2 = 'beta_0002';
          
        case 'corr'
          %% now just compute correlation from the smoothed images
          nii = load_untouch_nii(sfname);
          d = size(nii.img);
          nii_out = y;

          if strcmp(seedY,'mean')
            %dir2=[dir1,'/',compute,'/s',fwhm,'_',space,'_mean_',seedname_out];
            R = reshape(corr(reshape(nii.img,[],d(4))', meanY),d(1:3));
          else
            %dir2=[dir1,'/',compute,'/s',fwhm,'_',space,'_pc1_',seedname_out];
            R = reshape(corr(reshape(nii.img,[],d(4))', U(:,1)),d(1:3));
          end
          %[~,~]=mkdir(dir2);
          nii_out.img = R;
          fname2 = 'corr';
          save_untouch_nii(nii_out, [dir2,'/',fname2,'.nii']);
      end
      
      if strcmp(space,'native')
        %% now coreg it into indi-T1w
        exp3=EXP;
        exp3.subjID={subjid};
        exp3.prefix='c';
        exp3.name_fixed  =[dir1,'/t1w.nii'];
        exp3.name_moving =[dir1,'/meanuarest410.nii'];
        exp3.name_others ={[dir2,'/',fname2,'.nii']};
        myspm_coreg(exp3)
        
        %% now deform it to mni space
        exp4=[];
        exp4.fname_moving=[dir2,'/c',fname2,'.nii'];
        exp4.fname_def   =[dir1,'/y_t1w.nii'];
        exp4.vox_mm = 2;
        myspm_deform(exp4)
        fname2=['wc',fname2];
      end
      
      %% and link it to the mni directory
      unix(['ln -sf ',dir2,'/',fname2,'.nii ',dir_mni,'/',fname2,'_',subjid,'.nii']);
    end
  end
end
cd(dir0);
end