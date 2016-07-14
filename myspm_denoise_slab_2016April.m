function EXP = myspm_denoise_slab(EXP)
% EXP = myspm_denoise_slab(EXP)
%
% EXP requires:
%  .name_epi;
%  .cov_idx    1x1 for rp+cc+gs+scrub
%  .num
% (.nofigure)
%
% (cc) 2015, sgKIM.

% 0. set parameters
path0=pwd;
subjID = fsss_subjID(EXP.subjID);
TR_sec = EXP.TR_sec;
if ~isfield(EXP,'bpf2')
  EXP.bpf2 = [0.009, 0.08];
end
FilterBand = EXP.bpf2;
sumstat=[];

if ~isfield(EXP,'tprob')
  if isfield(EXP,'param_mask')
    if numel(EXP.param_mask)<3
      EXP.param_mask(3) = EXP.param_mask(2);
    end
    for c=1:3
      EXP.tprob{c} = sprintf('%0.2f',EXP.param_mask(c));
    end
  else
    EXP.tprob={'0.99','0.99','0.99'};
  end
end
tprob = EXP.tprob;
if ~isfield(EXP,'num_pcs'), EXP.num_pcs=6;   end;
if ~isfield(EXP,'detrend'), EXP.detrend=1;    end;
if ~isfield(EXP,'varnorm'), EXP.varnorm=1;    end;
if ~isfield(EXP,'bpf1'),    EXP.bpf1=[0,Inf]; end;
if ~isfield(EXP,'param_cc')
  EXP.param_cc = sprintf('n%df%0.2f-%0.2f', EXP.num_pcs, EXP.bpf1);
end

if ~isfield(EXP,'param_art')
  if isfield(EXP,'global_threshold') && isfield(EXP,'motion_threshold')
    EXP.param_art = sprintf('%0.1fstd_%0.1fmm', ...
      EXP.global_threshold, EXP.motion_threshold);
  else
    EXP.param_art='3.0std_0.5mm';
  end
end
if ~isfield(EXP,'t1w_suffix')
  EXP.t1w_suffix='t1w';
end
if ~isfield(EXP,'fname_gmmask')
  EXP.fname_gmmask=['fsc1',EXP.t1w_suffix,'_',tprob{1},'.nii'];
end

subjID = fsss_subjID(EXP.subjID);
num_subj = numel(subjID);
for n=1:num_subj
  EXP.subjid = subjID{n};
  idx=strfind(EXP.name_epi,'?');
  if ~exist(['cc_',EXP.name_epi(1:idx-1),num2str(EXP.num_runs),'_eigvec',EXP.param_cc,'.txt'],'file')
    exp1 = EXP;
    exp1.subjID =subjID{n};
    myy_compcor_slab(exp1);
  end
  for r=1:EXP.num_runs
    idx=strfind(EXP.name_epi,'?');
    name_epi = [EXP.name_epi(1:idx-1),num2str(r),EXP.name_epi(idx+1:end)];
    [~,name1,~]=fileparts(name_epi);
    subjid = subjID{n};
    path1=[fullfile(EXP.dir_base,subjid),'/'];
    EXP.fname_epi = fullfile(path1, name_epi);
    path2=[fullfile(EXP.dir_base,subjid),'/fig/'];
    cd(path1);
    [~,~]=mkdir(path2);
    nii = load_untouch_nii(EXP.fname_epi);
    d = nii.hdr.dime.dim(2:5);
    y = double(reshape(nii.img,[],d(4))');
    if ~isfield(EXP,'name_cc')
      if ~isfield(EXP,'bpf1'),    EXP.bpf1 = [0 Inf]; end
      if ~isfield(EXP,'num_pcs'), EXP.num_pcs = 6; end
      EXP.name_cc= ...
        ['cc_',EXP.name_epi(1:idx-1),num2str(r),'_eigvec_',EXP.param_cc,'.txt'];
    end
    
    if EXP.num_pcs >0
      cc  = load([path1,EXP.name_cc]); % text file
    else
      cc = [];
    end
    if isfield(EXP,'name_art')
      fname_art=[path1,EXP.name_art];
    else
      [~,res] = mydir(fullfile(path1,['art_out_mov_',name1,'*']));
      fname_art=res{1};
    end
    load(fname_art,'R');
    rmparam = R(:,(end-6):(end-1));
    rp= [rmparam, [0; l2norm(diff(rmparam,1))]];
    clear R
    fname_gs = [path1,'cc_gm',tprob{1},'.txt'];
    if exist('fname_gs','file')
      gs = load(fname_gs);
    else
      gs = ones(size(rmparam,1),1);
    end
    output_suffix=['_',EXP.param_cc,'_',EXP.param_art, ...
      '_b',num2str(FilterBand(1),2),'-',num2str(FilterBand(2),2)];
    
    %% plot#1: timeseries
    hf=figure('position',[ 1959          48         958        1092]);
    ax1=subplot(6,7,[1:6]); plot(rp(:,end));
    xlim([1,d(4)]); h=colorbar; set(h,'visible','off');
    ylabel('||dm/dt||_2 (mm)');
    
    sumstat.meanz=zeros(5,size(y,1));
    sumstat.meanabsz=zeros(5,size(y,1));
    fy = myy_filter(y, TR_sec, FilterBand);
    zyres = zscore(fy)';
    
    subplot(6,7,7+[1:6]); imagesc(zyres);
    sumstat.meanz(1,:) = mean(zyres,1);
    sumstat.meanabsz(1,:) = mean(abs(zyres),1);
    caxis([-2 2]); hc=colorbar; title(hc,'z-score');
    
    R={};
    V={};
    rndidx = randi(size(fy,2),[256,1]);
    R{1}=corr(fy(:,rndidx));
    V{1}=triuval(R{1});
    V{1}(isnan(V{1}))=eps;
    
    % setting user-defined regressors to compare
    M=[];
    if ~isfield(EXP,'covset');
      EXP.covset=[1 2 3 4];
      if isfield(EXP,'nogs')
        EXP.covset=[1 2 4 3];
      end
      % or [1 2 4 3] for rp+cc+art+gs
    end
    cov_vals={[ones(d(4),1), linspace(-1,1,d(4))',rp], cc, gs, rmparam};
    cov_names={'trend+rigidmotion','+compcor','+globalsignal','+scrubbing'};
    cov_names_short={'+td+rp',['+cc_th',tprob{2},'_n',num2str(EXP.num_pcs)],...
      ['+gs_',tprob{1}],['+scrb_',num2str(size(rmparam,2))]};
    cov_names_shorter={'+td+rp','+cc','+gs','+scrb'};
    tag1{1}='orig';
    tag2{1}='orig';
    for j=1:4
      COV{j}  = cov_vals{EXP.covset(j)};
      Mdesc{j} = cov_names{EXP.covset(j)};
      tag1{j+1} = cov_names_short{EXP.covset(j)};
      tag2{j+1} = cov_names_shorter{EXP.covset(j)};
    end
    output_suffix_short = output_suffix;
    output_suffix=[output_suffix [tag2{:}]];
    
    title(ax1,[subjid,'/',EXP.name_epi,':',output_suffix(2:end)],'interp','none');
    sumstat.tag=tag1;
    
    for i=1:4
      M=double([M COV{i}]);
      yres = myy_filter(y-M*(pinv(M'*M)*M'*y), TR_sec, FilterBand);
      subplot(6,7,7+7*i+[1:6]);
      zyres=zscore(yres)';
      imagesc(zyres);
      sumstat.meanz(1+i,:) = mean(zyres,1);
      sumstat.meanabsz(1+i,:) = mean(abs(zyres),1);
      caxis([-2 2]); hc=colorbar; title(hc,'z-score');
      if i==4, xlabel('TR'); end
      if i==2, ylabel(['Subsampled voxels with GM>',tprob{1}]); end
      
      subplot(6,7,7+7*i+7);
      imagesc(zscore(M));
      set(gca,'fontsize',10)
      title(Mdesc{i})
      R{i+1}=corr(yres(:,rndidx)+eps);
      V{1+i}=triuval(R{i+1});
      V{1+i}(isnan(V{1+i}))=eps;
    end
    
    colormap(sgcolormap('CKM'));
    name_figure='_timeseries.png';
    screen2png([path2,'resy',output_suffix,name_figure]);
    close(hf);
    if isfield(EXP,'dir_figure')
      [~,~]=mkdir(EXP.dir_figure);
      copyfile([path2,'resy',output_suffix,name_figure], ...
        [EXP.dir_figure,'/resy',output_suffix,name_figure(1:end-4),'_',subjid,'.png']);
    end
    if isfield(EXP,'plotuntil')&&EXP.plotuntil==1
      continue
    end
    
    %% plot#2: distribution
    hf=figure('position',[1972 33 933 1061]);
    subplot(5,3,[1,2]);
    estker=zeros(5,256);
    Momenta=zeros(5,4);
    for i=1:5
      xi=linspace(-1,1,256);
      estker(i,:) =ksdensity(V{i}, xi);
      Momenta(i,1) = mean(V{i});
      Momenta(i,2) = var(V{i});
      Momenta(i,3) = skewness(V{i});
      Momenta(i,4) = kurtosis(V{i});
    end
    sumstat.momenta=Momenta;
    hold on;
    ColorOrder=[1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 1 0.5 1; 0.5 1 1];
    colororder='rgbmc';
    for i=1:5
      h_lines(i)=plot(xi,estker(i,:),'linewidth',7-i,'color',colororder(i));
    end
    tag=tag1;
    h=legend(h_lines,tag, 'location','eastoutside');
    grid on; box on; set(h,'interp','none')
    title([subjid,'/',EXP.name_epi,':',output_suffix_short(2:end)],'interp','none');
    xlabel('Pearson correlation coefficient');ylabel({'estimated','kernel density'})
    
    tag=tag2;
    panelk=[4 7 10 13];
    MomentumDesc={'mean','variance','skewness','kurtosis'};
    for k=2:4
      subplot(5,3,panelk(k)); hold on;
      for j=1:5
        h=barh(j,Momenta(j,k),'facecolor',ColorOrder(j,:));
      end
      set(gca,'ytick',1:5, 'yticklabel',tag, 'ydir','rev');
      xlabel(MomentumDesc{k})
      box on; grid on;
    end
    subplot(5,3,panelk(1));
    barh([max(rp)],'facecolor',[.7 .7 .7]);
    yticklabel={'dx','dy','dz','rx','ry','rz','|dmdt|'};
    title('max')
    set(gca,'ydir','rev','yticklabel',yticklabel)
    xlabel('mm'); grid on; box on;
    
    tag=tag2;
    for i=1:5
      subplot(5,3,i*3)
      h=qqplot(V{i});
      title('');
      set(h(1),'MarkerEdgeColor',ColorOrder(i,:))
      set(h(3),'color','k','linestyle','-');
      [~,p]=kstest(zscore(V{i}));
      sumstat.kspval(i)=p;
      text(-4.7,0.8,['ks:gof=',num2str(p*100,3),'%']);
      grid on; box on;xlabel('norm-q');
      ylabel([tag{i},'-q'],'interp','none');
      axis([-5 5 -1 1]);
    end
    
    % save summary file
    fname_sum=[path1,'resy',output_suffix,'.mat'];
    save(fname_sum, 'sumstat');
    
    % save figure
    name_figure='_corrdist.png';
    screen2png([path2,'resy',output_suffix,name_figure]);
    close(hf);
    if isfield(EXP,'dir_figure')
      copyfile([path2,'resy',output_suffix,name_figure],...
        [EXP.dir_figure,'/resy',output_suffix,name_figure(1:end-4),'_',subjid,'.png']);
    end
    
    if isfield(EXP,'plotuntil')&&EXP.plotuntil==2
      continue
    end
    %% plot#4. sample timeseries (combine with plot#3?)
    if isfield(EXP,'sampletimeseries');
      % 1. find voxel indices for precuneous, med-prefrontal, precentral,
      aparc = load_untouch_nii(fname_aparc);
      %aseglabels=[12130 11130]; % precuneous
      %aseglabels=[12109 11109]; % post cingulate gyrus
      % 11146 left central sulcus
      % 12146 right central sulcus
      % 11107 12107 mid-anterior cingulate gyrus & sulcus
      % 11116 12116 superior frontal gyrus
      roiname={'Precuneous','Central sulcus','Anterior cingulate gyrus'};
      AsegLabels={[11130, 12130], [11146], [11107, 12107]}; %, [11116, 12116]};
      %ind1 = crrmap.ind; % seed voxel
      roi=[];
      roi(1).ijk = crrmap.ijk_med;
      roi(1).ind = crrmap.ind;
      l=2;
      for li=2:numel(AsegLabels)
        mask = ismember(aparc.img, AsegLabels{li});
        % z-coordinate should be on the same plane as z = crrmap.ijk_med(3)
        mask_sub = find3(mask);
        if ~numel(mask_sub(mask_sub(:,3)==crrmap.ijk_med(3)))
          warning('no voxel on this slice!');
          continue;
        else
          mask_sub = mask_sub(mask_sub(:,3)==crrmap.ijk_med(3),:);
          roi(l).ijk = round(median(mask_sub));
          if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
            % if the median is not in the label, find nearest precuneous voxel
            roi(l).dist = l2norm(mask_sub - repmat(roi(l).ijk,[size(mask_sub,1),1]));
            [~,ii] = min(roi(l).dist);
            roi(l).ijk = mask_sub(ii,:);
            % double check
            if ~ismember(aparc.img(roi(l).ijk(1),roi(l).ijk(2),roi(l).ijk(3)), AsegLabels{l})
              error('??!!')
            end
          end
          roi(l).ind = sub2ind(d(1:2), roi(l).ijk(1), roi(l).ijk(2)); % index for 2-D slice
          l=l+1;
        end
      end
      
      sc=2.5;
      corrthres=0.25;
      hf=figure('position',[1972, 33, d(2)*4*sc, d(1)*5*sc], 'color','k');
      t1w = load_untouch_nii([EXP.dir_base,subjid,'/oBrain.nii']);
      t1w = t1w.img(:,:,crrmap.ijk(3))';
      [~,ax] =getLayout(20,[5,4]);
      [~,ax2] =getLayout(20,[5,4],[0.10 0.15]);
      crrmap.ind_gm=find(gm.img(:,:,crrmap.ijk_med(3)));
      
      M=[];
      for i=1:5
        if i == 1
          crrmap.yres = myy_filter(crrmap.y, TR_sec, FilterBand);
        else
          M=double([M COV{i-1}]);
          crrmap.yres = myy_filter(crrmap.y-M*((M'*M)\M'*crrmap.y), TR_sec, FilterBand);
        end
        
        % spatial smoothing!
        crrmap.yres_3d = reshape(crrmap.yres',d(1),d(2),[]);
        numScans=size(crrmap.yres_3d,3);
        for t=1:numScans
          crrmap.yres_3d(:,:,t) = gaussblur(crrmap.yres_3d(:,:,t),1.5);
        end
        crrmap.yres = reshape(crrmap.yres_3d,[],numScans)';
        % correlation with mean of precuneous mask on the slice
        crrmap.rsfc = corr(crrmap.yres, mean(crrmap.yres(:,crrmap.idx_prec),2));
        
        % axial slice: correlation map
        h=axespos(ax,(i-1)*4+1);
        img1=double(t1w);
        img2=reshape(crrmap.rsfc', [d(1) d(2)])';
        imagecorr1(img1, img2, corrthres);
        hold on;
        line([crrmap.ijk_med(1);crrmap.ijk_med(1)],[0 d(2)], 'color','w');
        line([0 d(1)], [crrmap.ijk_med(2);crrmap.ijk_med(2)], 'color','w');
        if i==1, text(d(2)*0.3,d(1)*0.95,'FWHM=1.5 pixels', 'color','w', 'fontsize',12); end
        colororder='gmy';
        for l=1:numel(roi)
          scatter(roi(l).ijk(1), roi(l).ijk(2), 100, ...
            'MarkerFaceColor', colororder(l), 'lineWidth',2, 'MarkerEdgeColor','k');
        end
        text(d(1)*0.07, d(2)*0.07, 'R','color','w','fontsize',12)
        
        % individual timeseries
        xlim1=1;      xlim2=min(150, size(crrmap.yres,1));
        if isfield(EXP,'xlim')
          xlim1=EXP.xlim(1);    xlim2=EXP.xlim(2);
        end
        for l=1:numel(roi)
          h=axespos(ax2,(i-1)*4+l+1);
          plot(crrmap.yres(:, roi(l).ind),'color',colororder(l))
          xlim([xlim1 xlim2])
          set(gca,'color','k','xcolor','w','ycolor','w');
          hold on;
          ylim1=ylim;
          text(xlim1+5,ylim1(2)-diff(ylim1)*0.07, roiname{l}, 'color',colororder(l))
          if l==1
            title(tag{i},'fontsize',15,'color','w','interp','none');
          end
          grid on;
          if i==1 && l==1,  xlabel('TR','fontsize',9); end
        end
      end
      
      h=axes('position',[0 0 0.25 0.02]);
      hb=colorbar('peer',h,'location','South');
      colormap(sgcolormap('BLUE-RED'));
      axis off; caxis([1,7]);
      set(hb,'xtick',1:7,'xticklabel',{'-1','-0.75','-0.5','|0.25|','0.5','0.75','1'});
      cbfreeze(hb);
      
      colormap(sgcolormap('CKM'));
      name_figure='_corrmap_inditimeseries.png';
      screen2png([path2,'resy',output_suffix,name_figure]);
      close(hf);
      if isfield(EXP,'dir_figure')
        copyfile([path2,'resy',output_suffix,name_figure],...
          [EXP.dir_figure,'/resy',output_suffix,name_figure(1:end-4),'_',subjid,'.png']);
      end
    end
    
    %% final decision?
    if ~isfield(EXP,'cov_idx')
      disp([subjid,': enter EXP.cov_idx to select regressors and save a residual image']);
      continue
    else
      nii = load_untouch_nii(EXP.fname_epi);
      % datatype=4, bitpix=16
      y = single(reshape(nii.img,[],d(4))');
      M=[];
      for i=1:EXP.cov_idx
        M=single([M COV{i}]);
      end
      % LSE (y-y_hat) and band-pass filtering
      yres = myy_filter(y-M*(pinv(M'*M)*M'*y), TR_sec, FilterBand);
      % don't forget to enter space x time, instead of time x space
      if  (nii.hdr.dime.bitpix == 16) && (d(4)>1000)
        % datatype: 8=int16(short), 16=single, 64=double
        [yres,offset,slope] = single2int16(yres);
        nii.hdr.hist.descrip = sprintf(['motion artifacts regressed out: %s, ', ...
          'single2int16:offset=%f, slope=%f'], output_suffix, offset, slope);
      else
        nii.hdr.hist.descrip = sprintf(['motion artifacts regressed out: %s'], ...
          output_suffix);
      end
      nii.img = reshape(yres',d);
      pcsnum=num2str(EXP.num_pcs);
      if isfield(EXP,'out_prefix')
        out_prefix=EXP.out_prefix;
      else
        out_prefix=['fr',pcsnum];
      end
      fname_out=[path1,out_prefix,name_epi];
      disp(['> saving residual in ',fname_out,'..']);
      save_untouch_nii(nii, fname_out);
      
      if isfield(EXP,'dir_fs')
        % and (liberal) mask the residual images to find where's where
        setenv('FSLOUTPUTTYPE','NIFTI');
        [~,n1,~] = fileparts_gz (fname_out);
        [~,n2,~] = fileparts_gz (name_epi);
        name2=['mean',n2];
        unix(['fslmaths ',n2,' -Tmean ',name2]);
        unix(['bet ',name2,' ',name2,'_brain -R -f 0.1'])
        unix(['fslmaths ',name2,'_brain -bin epimask.nii'])
        unix(['mv ',name2,'_brain* /tmp'])
        unix(['fslmaths ',n1,' -mas epimask m',n1])
      end
    end
  end
end
cd(path0);

end

function [F,ker] = gaussblur(f,FWHM)
%function F = gaussblur(f,FWHM)
%
% f : input 1D series, 2D image or 3D volume
% FWHM: filter size in sampling point (pixel/voxel).
%
% (c) Moo K. Chung, July 2003. March 2004.
% (cc) sgKIM, 2011, 2015
% + now works for 1D, 2D and 3D.
% + kernel size extended to ensure the smoothness as intended by FWHM.
% + relationship between t and FWHM is corrected now.
%
% See Chung et al., 2001. A Unified Statistical Appraoch to
% Deformation-based morphometry, NeuroImage for implementation detail.


% ksize is the number of neighboring voxels we are averaging.
% it depends on FWHM.
% ksize is allowd to be an odd interger.
ksize = round(FWHM*2+1);

% see Chung et al.(2001) for the relationship between t and FWHM.
% FWHM = 4*sqrt(log(2))*sqrt(t)
%
% FWHM = sqrt(8*log(2)) * sigma (mathworld.Wolfram.com)
% sigma = FWHM / sqrt(8*log(2))
% t = (sigma^2)/2
t = (FWHM^2) / (16*log(2));

d = round(ksize/2);

% we truncate and normalize the Gaussian kernel for fast computation
% 1D, 2D, 3D version of kerenl is slightly different
n = ndims(f);
if n==1
  domain_x=(1:ksize)-d;
  ker=exp(-(domain_x.^2)/(4*t)) /sqrt(t/2);  % not in exponential, so doesn't matter due to normalization
elseif n==2
  domain_x=kron(ones(ksize,1),[1:ksize])-d;
  domain_y=kron(ones(1,ksize),[1:ksize]')-d;
  ker=exp(-(domain_x.^2 +domain_y.^2)/(4*t)) /t;
elseif n==3
  domain_x = repmat ( kron(ones(ksize,1),[1:ksize])-d, [1,1,ksize]);
  domain_y = repmat ( kron(ones(1,ksize),[1:ksize]')-d, [1,1,ksize]);
  domain_z = reshape( repmat ( kron(ones(ksize,1),[1:ksize])-d, [ksize, 1, 1]), [ksize,ksize,ksize]);
  ker=exp(-(domain_x.^2 + domain_y.^2 +domain_z.^2)/(4*t)) /t^(3/2);
end;
ker=ker/sum(ker(:));

% smoothing is done by convolution
F=convn(f,ker,'same');

end
