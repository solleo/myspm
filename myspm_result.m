function EXP = myspm_result(EXP)
% EXP = myspm_result(EXP)
%
% It creates SPM result reports automatically for given a directory (where you have SPM.mat)
% if there is a suprathreshold voxels/clusters,
% It generates orthogonal overlay sections for each local maximum.
%
% EXP requires:
%  .dir_glm         'Nx1' directory to save SPM results
% (.mask)        'Nx1'  filename for an explicit mask
% (.thres.desc)     'Nx1'  'FWE','none', or 'cluster'(default)
% (.thres.alpha)    [1x1]  alpha level (default=0.05)
% (.thres.extent)   [1x1]  extent threshold of clusters in voxels (default=0)
% (.thres.clusterInitAlpha)  [1x1] cluster forming height threshold (default=0.001)
% (.thres.clusterInitExtent) [1x1] cluster forming extent (in voxels) threshold (default=10)
% (.fname_struct)   'Nx1' fullpath filename for background anatomical image for orthogonal slices
%                         (defulat='$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz')
% (.titlestr)       {1xNcont} Title text for SPM result report (default={'positive','negative'})
% (.dir_sum)        'Nx1' a summary directory into where you want to copy significant results
% (.append)         [1x1] whether to append results into an existing table (default=0)
% (.print)          [1x1] whether to generate figures and tables (default=1)
% (.atlas)          'Nx1' atlas to find anatomical names: 'fsl' (default) or 'spm12'
% (.MIP_only)       [1x1] prints out only MIP then exit (no slices/plots)
% (.showlessNC)     [1x1] limit # of clusters per contrast
%
% (cc) 2015-2018, sgKIM. solleo@gmail.com, https://ggooo.wordpress.com/
if nargin < 1,  EXP=[];  end
%% Setting default parameters
%if ~isfield(EXP,'mygraph'), EXP.mygraph.y_name='y'; EXP.mygraph.x_name='x'; end
if ~isfield(EXP,'append'),  EXP.append = 0;  end
if ~isfield(EXP,'print'),   EXP.print = 1;   end
if ~isfield(EXP,'thres')||~isfield(EXP.thres,'desc'),  EXP.thres.desc='cluster';  end
if ~isfield(EXP,'dir_glm'), EXP.dir_glm=pwd; end
if ~isfield(EXP,'suppressPlot'), EXP.suppressPlot=0; end
if isfield(EXP,'masking') && ~isfield(EXP,'mask'), EXP.mask = EXP.masking; end

% set table filename
today = datestr(now,'yyyymmmdd');
fname_tab = fullfile(EXP.dir_glm,['spm_',today,'.csv']);
fmt = {   '%s',   '%s',   '%0.3f', '%s','%0.3f','%0.3f','%-0.3f','%-0.3f','%-0.3f', '%0.0f','%3.0f','%3.0f','%3.0f',  '%s','%0.2f'};
tab_fmt = cell2fmt (fmt);
if ~EXP.append && ~~EXP.print
  fid = fopen(fname_tab,'w');
  hdr_fmt = 'measure\tcontrast_name\teffect_size\tstat\tpeak\tpeakZ\tuncor_pval\tcor_pval_peak\tcor_pval_clus\tK_E\tx_mm\ty_mm\tz_mm\tpeak_strc_name\tpeak_strc_prob\n';
  fprintf(fid, hdr_fmt, EXP.thres.desc);
  fclose(fid);
end
% if ~isfield(EXP,'fname_prefix')
%   prefix='';
% else
%   prefix = [EXP.fname_prefix,'_'];
% end
% if ~isfield(EXP,'model_desc')
%   [~,b,~]=fileparts(EXP.dir_glm);
%   EXP.model_desc=b;
% end

% set figure filename
% opts=spm_print('format','pdf'); % PDF doesn't stack????
if ~isfield(EXP,'fname_spm_fig')
  EXP.fname_spm_fig = fullfile(EXP.dir_glm,[EXP.model_desc,'_spm_',today,'.ps']);
  EXP.fname_spm_fig = strrep(EXP.fname_spm_fig,'>','-gt-');
  EXP.fname_spm_fig = strrep(EXP.fname_spm_fig,'<','-lt-');
end
% 
% % delete previous figure files today (rewriting)
% if EXP.print && exist(EXP.fname_spm_fig,'file')
%   delete(EXP.fname_spm_fig);
% end

% setting a thresholding code
if strcmpi('cluster',EXP.thres.desc)
  method=3;
elseif strcmpi('FWE',EXP.thres.desc)
  method=1;
elseif strcmpi('none',EXP.thres.desc)
  method=2;
end

% setting default thresholding parameters
if ~isfield(EXP.thres,'alpha'),    EXP.thres.alpha=0.05;  end
if ~isfield(EXP.thres,'extent'),   EXP.thres.extent=0;    end
if ~isfield(EXP.thres,'override'), EXP.thres.override=0;  end
% the conventional cluster-forming alpha is double-O (0.001),
if ~isfield(EXP.thres,'clusterInitAlpha')
  EXP.thres.clusterInitAlpha  = 0.001; end
if ~isfield(EXP.thres,'clusterInitExtent')
  EXP.thres.clusterInitExtent = 10; end

% read SPM.mat to find numbers of sessions (NumSess), conditions (NumCond), and
% contrasts (NumCntrst)
load([EXP.dir_glm,'/SPM.mat']);
NumSess = size(SPM.xY.VY,2);
NumCond = size(SPM.xX.X,2)/NumSess-7;
if ~isfield(EXP,'NumCntrst')
  EXP.NumCntrst=numel(SPM.xCon);
end
NumCntrst = EXP.NumCntrst;
% set title texts
% if ~isfield(EXP,'titlestr')
for i=1:NumCntrst
  EXP.titlestr{i}=SPM.xCon(i).name;
end
% end

% find structure image for orthogonal slices
if ~isfield(EXP,'fname_struct')
  fsldir=getenv('FSLDIR');
  EXP.fname_struct = fullfile(fsldir,'data','standard','MNI152_T1_1mm.nii.gz');
end
if strcmp(EXP.fname_struct,'conmus')
  EXP.fname_struct='/home/raid2/skim/MNI152_T1_2.5mm.nii';
elseif strcmp(EXP.fname_struct,'conmus1')
  EXP.fname_struct = '/scr/vatikan4/conmus3/GLM/MNI152_T1_1mm_masked_old.nii';
elseif strcmp(EXP.fname_struct,'conmus1_masked')
  EXP.fname_struct = '/scr/vatikan4/conmus3/GLM/MNI152_T1_1mm_masked.nii';
end
[~,name1,ext1]=fileparts(EXP.fname_struct);
if strcmp(ext1,'.gz')
  if ~exist(['/tmp/',name1],'file')
    gunzip(EXP.fname_struct, '/tmp/');
  end
  EXP.fname_struct = ['/tmp/',name1];
  fprintf('Structure image file: ');
  ls(EXP.fname_struct)
end

%% Now create result reports
spm('defaults','fmri');
spm_jobman('initcfg');
EXP.minp=zeros(1,NumCntrst);
TotalNC=0;
for cntrst=1:NumCntrst %numel(EXP.titlestr) % for each contrast
  matlabbatch={};
  matlabbatch{1}.spm.stats.results.spmmat = {[EXP.dir_glm,'/SPM.mat']};
  matlabbatch{1}.spm.stats.results.conspec(1).contrasts = cntrst;
  if method == 3
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'none';
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = EXP.thres.clusterInitAlpha;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thres.clusterInitExtent;
  else
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = EXP.thres.desc;
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = EXP.thres.alpha;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thres.extent;
  end
  if isfield(EXP,'mask')
    matlabbatch{1}.spm.stats.results.conspec(1).mask.image.name = {EXP.masking};
    matlabbatch{1}.spm.stats.results.conspec(1).mask.image.mtype = 0;
  end
  matlabbatch{1}.spm.stats.results.conspec(1).titlestr = EXP.titlestr{cntrst};
  matlabbatch{1}.spm.stats.results.units = 1;
  
  % this is what you need when you select cluster-threshold
  if method ==3
    matlabbatch{1}.spm.stats.results.print = false;
    spm_jobman('run', matlabbatch)
    
    % now read the goddamn table! :D
    TabDat = evalin('base','TabDat;');
    
    % find extent threshold for cluster-p < alpha
    Npeaks=size(TabDat.dat,1);
    K=[];
    for pi=1:Npeaks
      if (TabDat.dat{pi,3} < EXP.thres.alpha)
        K=[K TabDat.dat{pi,5}];
      end
    end
    if K
      EXP.thres.clusextent(cntrst)=min(K);
      matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thres.clusextent(cntrst) -1;
    end
  end
  
  % result table printing
  matlabbatch{1}.spm.stats.results.print = ~~EXP.print;
  
  %% Result Reports with any given threshold
  save([EXP.dir_glm,'/result.mat'],'matlabbatch');
  spm_jobman('run', matlabbatch);
  
  %% Get all variables results report
  xSPM= evalin('base','xSPM;');
  save('xSPM.mat','xSPM') % how large?
  hReg= evalin('base','hReg;');
  TabDat = evalin('base','TabDat;');
  SPM= evalin('base','SPM;');
  
  %% Now find significant peaks
  % now read the goddamn table! :D
  
  % first just export the table (for each contrast)
  save([EXP.dir_glm,'/TabDat',num2str(cntrst),'.mat'], 'TabDat');  % for matlab
  fid = fopen([EXP.dir_glm,'/TabDat',num2str(cntrst),'.csv'],'w'); % for other applications
  fprintf(fid, cell2fmt(TabDat.hdr(1,:)));
  fprintf(fid, cell2fmt(TabDat.hdr(2,:)));
  fprintf(fid, strrep(strrep(strrep(cell2fmt(TabDat.hdr(3,:)),'\it',''),'\rm_',''),'\equiv',''));
  
  C=[7 11 3];
  Npeaks=size(TabDat.dat,1);
  COORDS={};
  PI=[];
  pvals=[];
  tmax=[];zmax=[];
  for pi=1:Npeaks
    pvals=[pvals TabDat.dat{pi,C(method)}];
    tmax=[tmax TabDat.dat{pi,9}];
    zmax=[zmax TabDat.dat{pi,10}];
    Ncols=size(TabDat.dat,2);
    for col=1:Ncols
      fprintf(fid, TabDat.fmt{col},TabDat.dat{pi,col});
      if col<Ncols
        fprintf(fid, '\t');
      else
        fprintf(fid, '\n');
      end
    end
    if (TabDat.dat{pi,C(method)} < EXP.thres.alpha)
      COORDS=[COORDS TabDat.dat{pi,end}];
      PI=[PI pi];
    end
  end
  fclose(fid);
  [~,b]=find( ~isinf(tmax),1,'first');
  cmax=tmax(b);
  EXP.cmax = cmax;
  
  % now return min P for each contrast
  EXP.minp(cntrst)=min([pvals,1]);
  NC=numel(COORDS);
  
  %% Skip the last for only MIP
  if EXP.print
    spm_print(EXP.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
  end
  if isfield(EXP,'MIP_only')
    continue
  end
  
  %% Now save "significant" cluster maps (from spm code somewhere...)
  if NC
    XYZ  = xSPM.XYZ;
    Z       = round(spm_clusters(XYZ));
    num     = max(Z);
    [n, ni] = sort(histc(Z,1:num), 2, 'descend');
    n       = size(ni);
    n(ni)   = 1:num;
    Z       = n(Z);
    fname_sigclus=['sigclus_',EXP.titlestr{cntrst}];
    fname_sigclus = strrep(fname_sigclus,'>','-gt-');
    fname_sigclus = strrep(fname_sigclus,'<','-lt-');
    fname_sigclus = strrep(fname_sigclus,' ','');
    spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
      sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
      [fname_sigclus,'.nii']);
  end
  
  %% create orthogonal sections and reposition onto peaks
  redcmap=[ gray(64);
    [linspace(0.5667,1,14)', linspace(0,0,14)', linspace(0,0,14)'];   %  1~14 [14] dark red
    [linspace(1,1,30)', linspace(0.03333,1,30)', linspace(0,0,30)'];  % 15~44 [30] red to yellow
    [linspace(1,1,20)', linspace(1,1,20)', linspace(0.05,1,20)']      % 45~64 yellow
    ];
  
  bluecmap=[ gray(64);
    [linspace(0,0,14)', linspace(0,0,14)', linspace(0.5667,1,14)'];   %  1~14 [14] dark blue
    [linspace(0,0,30)', linspace(0.03333,1,30)', linspace(1,1,30)'];  % 15~44 [30] blue to cyan
    [linspace(0.05,1,20)', linspace(1,1,20)', linspace(1,1,20)']      % 45~64 cyan
    ];
  if ~isfield(EXP,'showlessNC')
    showlessNC = 30;
  else
    showlessNC = EXP.showlessNC;
  end
  if NC>showlessNC
    warning(['Too many clusters, only ',num2str(showlessNC),' clusters are visualized']);
    NC=showlessNC;
  end
  if NC % if any cluster exists
    for ci=1:NC
      spm_sections(xSPM,hReg,EXP.fname_struct);
      spm_orthviews('reposition',COORDS{ci});
      global st
      st.vols{1}.blobs{1}.max=EXP.cmax;
      spm_orthviews('Xhirs','on');
      spm_orthviews('redraw');
      
      % odd/even index of contrast
      if mod(cntrst,2)
        colormap(redcmap);
      else
        colormap(bluecmap);
      end
      if strcmpi(xSPM.STAT,'F')
        colormap(redcmap);
      end
      
      if EXP.print
        spm_print(EXP.fname_spm_fig);
      end
      
      % and plots
      if ~strcmpi(xSPM.STAT,'F')
        % SCATTER PLOTS for T-contrast
        xXi = find(SPM.xCon(cntrst).c~=0);
        cfg=[]; cfg.Ic=cntrst; cfg.xXi=xXi; cfg.contrast = SPM.xCon(cntrst).c;
        cfg.origoffset=1;
        cfg.x_name = SPM.xCon(cntrst).name;
        cfg.y_name = 'y';
        if isfield(EXP,'y_name'), cfg.y_name=EXP.y_name; end
        if isfield(EXP,'markCorrThres')
          cfg.markCorrThres=EXP.markCorrThres;
          if isfield(EXP,'CorrFDRthres')
            cfg.markCorrThres.CorrFDRThres=EXP.CorrFDRthres;
          end
        end
        cfg.atlas='fsl';
        if isfield(EXP,'atlas')
          if strfind(lower(EXP.atlas),'spm')
            cfg.atlas='spm12';
          elseif strfind(lower(EXP.atlas),'fsl')
            cfg.atlas='fsl';
          end
        end
        cfg.clusterIdx = ci;
        if isfield(SPM,'Sess')
          cfg.isfmri=1;
        else
          cfg.isfmri=0;
        end
        if isfield(EXP,'suppressPlot'), cfg.suppressPlot=EXP.suppressPlot; end
        [~,~,beta,~,STRC,thres,peakxyz]=myspm_graph(xSPM,SPM,hReg,cfg); % this only reads peak
        if ~isempty(thres)
          EXP.CorrFDRthres = thres;
        end
        if EXP.print
          spm_print(EXP.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        end
        %         %% fMRI: bar graph with 90% CI
        %         if any(contains({SPM.xCon.name},'Effect of interest'))
        %           cfg2 = cfg;
        %           cfg2.plotType = 'Contrast estimates and 90% C.I. effect of interest';
        %           myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
        %           if EXP.print
        %             spm_print(EXP.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        %           end
        %         end
        %%
      else % for F-contrast
        % only contrast estimates plot:
        cfg2 = cfg;
        cfg2.Ic = cntrst;
        cfg2.plotType = 'Contrast estimates and 90% C.I.';
        myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
        if EXP.print
          spm_print(EXP.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
        end
      end
      
      %% fMRI: PSHT?
      %         cfg2 = cfg;
      %         cfg2.plotType = 'Event-related responses';
      %         myspm_graph(xSPM,SPM,hReg,cfg2); % this only reads peak
      %         if EXP.print
      %           spm_print(EXP.fname_spm_fig); % add MIP + peak table in sigcluster .PDF
      %         end
      
      %%
      idx_x = find(SPM.xCon(cntrst).c);
      idx_x = idx_x(1);
      % generate a summary table!
      fid = fopen(fname_tab, 'a');
      % 'y-name  x-name  effect  stattype  peak%s\tpeakZ  uncor_pval \t cor_pval(%s) \tK_E\tMNI-x_mm\tMNI-y_mm\tMNI-z_mm\tpeak-Strc-name\npeak-Strc-prob\clus-nStrc-name\nclus-Strc-prob\n'
      fprintf(fid,tab_fmt,...
        cfg.y_name, cfg.x_name, beta(idx_x), xSPM.STAT, ...
        TabDat.dat{PI(ci),9}, TabDat.dat{PI(ci),10},  ...
        TabDat.dat{PI(ci),11}, TabDat.dat{PI(ci),7}, TabDat.dat{PI(ci),3}, ...
        TabDat.dat{PI(ci),5}, ...
        TabDat.dat{PI(ci),12}(1), TabDat.dat{PI(ci),12}(2), TabDat.dat{PI(ci),12}(3), ...
        STRC.strc.name, STRC.strc.prob); %, ...
      fclose(fid);
    end % of each cluster
    TotalNC=TotalNC+NC;
  end % if any cluster exists
end % of each contrast

%% convert ps files to pdf
myps2pdf(EXP.fname_spm_fig)
% gsbin='/usr/bin/ghostscript';
% if ~exist(gsbin,'file')
%   gsbin='/opt/local/bin/gs';
% end
% ps2pdf('psfile',EXP.fname_spm_fig, ...
%   'pdffile',[EXP.fname_spm_fig(1:end-2),'pdf'],...
%   'gscommand',gsbin)
% delete(EXP.fname_spm_fig)

% %% copy 'significant' results to dir_sum
% if TotalNC && isfield(EXP,'dir_sum')
%   [~,name1,~]=fileparts(EXP.dir_glm);
%   [~,~]=mkdir(EXP.dir_sum);
%   today=datestr(now,'yyyymmmdd');
%   src=fullfile(EXP.dir_glm,['spm_',today,'.ps']);   % whole brain table?
%   if ~isfield(EXP,'prefix_sum'), EXP.prefix_sum=''; end
%   trg=[EXP.dir_sum,'/',EXP.prefix_sum,'_spm_',today,'_',name1,'.ps'];
%   system(['cp ',src,' ',trg]);
%
%   src=EXP.fname_spm_fig;
%   [~,b,c] = fileparts(src);
%   if ~isfield(EXP,'prefix_sum'), EXP.prefix_sum=''; end
%   trg=[EXP.dir_sum,'/',EXP.prefix_sum,'_',b,c];
%   system(['cp ',src,' ',trg]);
%
%   fname_sumtab = [EXP.dir_sum,'/summary.csv'];
%   if ~exist(fname_sumtab,'file')
%     fid = fopen(fname_sumtab,'w');
%     fprintf(fid,hdr_fmt, EXP.thres.desc);
%     fclose(fid);
%   end
%   system(['tail -n +2 ',fname_tab,' > /tmp/tab2.csv']);
%   system(['cat ',fname_sumtab,' /tmp/tab2.csv  >/tmp/sumtab.csv']);
%   system(['mv /tmp/sumtab.csv ',fname_sumtab]);
% end
end
