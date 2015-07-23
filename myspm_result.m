function EXP = myspm_result(EXP)
% EXP = myspm_result(EXP)
%
% It creates SPM result reports automatically for given a directory (where you have SPM.mat)
% if there is a suprathreshold voxels/clusters,
% It generates orthogonal overlay sections for each local maximum.
%
% EXP.dir_name = '' (the directory where you have SPM.mat)
% EXP.thresh.desc = either 'FWE','none', or 'cluster'
% EXP.thresh.alpha = 0.05 (default)
% EXP.thresh.extent = 0 (# fo voxels; default)
% EXP.thresh.clusterInitAlpha = 0.001 (default) for a cluster-threshold
% EXP.thresh.clusterInitExtent = 10 (voxels; default)
% EXP.titlestr = {'positive','negative'} (default)
% EXP.fname_struct = '$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz' (by default)
% EXP.dir_sum = '' (a summary directory where you want to copy 'significant'
% results into)
% EXP.append = 0 (default)
% EXP.print  = 1 (default)
% EXP.mygraph.y_name = 'x'
% EXP.mygraph.x_name = 'y' for the scatterplots and summary tables
% EXP.atlas = 'fsl' (default) or 'spm12'
%
% (cc) 2015, sgKIM. solleo@gmail.com, https://ggooo.wordpress.com/

%% Always print out the annotated tables
if ~isfield(EXP,'mygraph')
  EXP.mygraph.y_name='y';
  EXP.mygraph.x_name='x';
end

%% EXP.append?

if ~isfield(EXP,'append')
  EXP.append = 0;
end
if ~isfield(EXP,'print')
  EXP.print = 1;
end
today=datestr(now,'yyyymmmdd');
fname_tab=fullfile(EXP.dir_name,['spm_',today,'.csv']);
%fname_tab=fullfile(EXP.dir_name,['spm_',today,'.xls']);
%fmt={   '%s',   '%s',   '%0.3f', '%s','%0.3f','%0.3f','%-0.3f','%-0.3f','%-0.3f', '%0.0f','%3.0f','%3.0f','%3.0f',  '%s','%0.2f',  '%s','%0.2f'};
fmt={   '%s',   '%s',   '%0.3f', '%s','%0.3f','%0.3f','%-0.3f','%-0.3f','%-0.3f', '%0.0f','%3.0f','%3.0f','%3.0f',  '%s','%0.2f'};
tab_fmt= cell2fmt (fmt);

if EXP.append==0 && EXP.print==1
  src=fullfile(EXP.dir_name,['spm_',today,'.ps']);
  system(['rm -f ',src]);
  fid=fopen(fname_tab,'w');
  %hdr_fmt='y-name\tx-name\teffect\tstat\tpeak\tpeakZ\tuncor_pval\tcor_pval(peak)\tcor_pval(clus)\tK_E\tMNI-x_mm\tMNI-y_mm\tMNI-z_mm\tpeak-strc-name\tpeak-strc-prob\tclus-strc-name\tclus-strc-prob\n';
  hdr_fmt='y-name\tx-name\teffect\tstat\tpeak\tpeakZ\tuncor_pval\tcor_pval(peak)\tcor_pval(clus)\tK_E\tMNI-x_mm\tMNI-y_mm\tMNI-z_mm\tpeak-strc-name\tpeak-strc-prob\n';
  fprintf(fid, hdr_fmt, EXP.thresh.desc);
  fclose(fid);
  %xlswrite(fname_tab, hdr_fmt, 1, 'A1');
end

%% Result Reports with any given threshold
spm('defaults','fmri');
spm_jobman('initcfg');

if strcmpi('cluster',EXP.thresh.desc)
  method=3;
elseif strcmpi('FWE',EXP.thresh.desc)
  method=1;
elseif strcmpi('none',EXP.thresh.desc)
  method=2;
else
  error('EXP.thresh.desc should be either ''cluster'', ''FWE'', ''none''');
end

if ~isfield(EXP.thresh,'alpha')
  EXP.thresh.alpha=0.05;
end

if ~isfield(EXP.thresh,'extent')
  EXP.thresh.extent=0;
end

if ~isfield(EXP.thresh,'override')
  EXP.thresh.override=0;
end

if ~isfield(EXP.thresh,'clusterInitAlpha')
  EXP.thresh.clusterInitAlpha = 0.001;
end
if ~isfield(EXP.thresh,'clusterInitExtent')
  EXP.thresh.clusterInitExtent = 10;
end

%[TODO]: make it adaptable for other contrasts!
if ~isfield(EXP,'titlestr')
  EXP.titlestr={'positive','negative'};
end

if ~isfield(EXP,'fname_struct')
  fsldir=getenv('FSLDIR');
  EXP.fname_struct = fullfile(fsldir,'data','standard','MNI152_T1_1mm.nii.gz');
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

EXP.minp=zeros(1,numel(EXP.titlestr));
TotalNC=0;
for cntrst=1:numel(EXP.titlestr)
  matlabbatch={};
  matlabbatch{1}.spm.stats.results.spmmat = {[EXP.dir_name,'/SPM.mat']};
  matlabbatch{1}.spm.stats.results.conspec(1).contrasts = cntrst;
  if method == 3
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = 'none';
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = EXP.thresh.clusterInitAlpha;% 0.001;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thresh.clusterInitExtent; %10;
  else
    matlabbatch{1}.spm.stats.results.conspec(1).threshdesc = EXP.thresh.desc;
    matlabbatch{1}.spm.stats.results.conspec(1).thresh = EXP.thresh.alpha;
    matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thresh.extent;
  end
  matlabbatch{1}.spm.stats.results.conspec(1).mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
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
      if (TabDat.dat{pi,3} < EXP.thresh.alpha)
        K=[K TabDat.dat{pi,5}];
      end
    end
    if K
      EXP.thresh.clusextent(cntrst)=min(K);
      matlabbatch{1}.spm.stats.results.conspec(1).extent = EXP.thresh.clusextent(cntrst) -1;
    end
  end
  % result table printing
  matlabbatch{1}.spm.stats.results.print = ~~EXP.print;
  spm_jobman('run', matlabbatch)
  
  %% Get all variables results report
  xSPM= evalin('base','xSPM;');
  hReg= evalin('base','hReg;');
  TabDat = evalin('base','TabDat;');
  SPM= evalin('base','SPM;');
  
  %% Now find significant peaks
  % now read the goddamn table! :D
  
  % first just export the table (for each contrast)
  save([EXP.dir_name,'/TabDat',num2str(cntrst),'.mat'], 'TabDat');  % for matlab
  fid = fopen([EXP.dir_name,'/TabDat',num2str(cntrst),'.csv'],'w'); % for other applications
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
    if (TabDat.dat{pi,C(method)} < EXP.thresh.alpha)
      COORDS=[COORDS TabDat.dat{pi,end}];
      PI=[PI pi];
    end
  end
  fclose(fid);
  
  [~,b]=find( ~isinf(zmax),1,'first');
  cmax=tmax(b);
  if ~isfield(EXP,'cmax')
    EXP.cmax = cmax;
  end
  
  % now return min P for each contrast
  EXP.minp(cntrst)=min([pvals,1]);
  
  NC=numel(COORDS);
  if (NC > 50) && ~EXP.thresh.override
    error('Seriously, do you want 50+ figures for 50+ blobs? But if you really want that, set EXP.thresh.override=1');
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
    if numel(EXP.titlestr) > 2
      error('Now code for 2+ contrasts!');
    end
    Sign='+-';
    try fname_sigclus=['sigclus_',Sign(cntrst),EXP.vi.name];
    catch ME
      fname_sigclus=['sigclus_',Sign(cntrst),'1'];
    end
    spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
      sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k), ...
      fname_sigclus);
  end
  
  
  %% create orthogonal sections and reposition onto peaks
  redcmap=[ gray(64);
    [linspace(0.5667,1,14)', linspace(0,0,14)', linspace(0,0,14)'];        %  1~14 [14]
    [linspace(1,1,30)', linspace(0.03333,1,30)', linspace(0,0,30)'];  % 15~44 [30]
    [linspace(1,1,20)', linspace(1,1,20)', linspace(0.05,1,20)']           % 45~64 [
    ];
  
  bluecmap=[ gray(64);
    [linspace(0,0,14)', linspace(0,0,14)', linspace(0.5667,1,14)'];        %  1~14 [14]
    [linspace(0,0,30)', linspace(0.03333,1,30)', linspace(1,1,30)'];  % 15~44 [30]
    [linspace(0.05,1,20)', linspace(1,1,20)', linspace(1,1,20)']           % 45~64 [
    ];
  
  if NC
    for ci=1:NC
      spm_sections(xSPM,hReg,EXP.fname_struct);
      spm_orthviews('reposition',COORDS{ci});
      global st
      st.vols{1}.blobs{1}.max=EXP.cmax;
      spm_orthviews('redraw');
      if sum(SPM.xCon(cntrst).c)<0
        colormap(bluecmap);
      else
        colormap(redcmap);
      end
      if EXP.print
        if isfield(EXP,'fname_spm_fig')
          spm_print(EXP.fname_spm_fig);
        else
          spm_print
        end
      end
      
      % and scatter plot!
      xXi = find(SPM.xCon(cntrst).c~=0);
      cfg=[]; cfg.Ic=cntrst; cfg.xXi=xXi;
      cfg.origoffset=1;
      if isfield(EXP,'mygraph')
        if isfield(EXP.mygraph,'x_name')
          cfg.x_name = EXP.mygraph.x_name;
        end
        if isfield(EXP.mygraph,'y_name')
          cfg.y_name = EXP.mygraph.y_name;
        end
      end
      if isfield(EXP,'markCorrThres')
        cfg.markCorrThres=EXP.markCorrThres;
        if isfield(EXP,'CorrFDRthres')
          cfg.markCorrThres.CorrFDRThres=EXP.CorrFDRthres;
        end
      end
      cfg.atlas='spm12';
      if isfield(EXP,'atlas')
        if strfind(EXP.atlas,'spm')
          cfg.atlas='spm12';
        elseif strfind(EXP.atlas,'fsl')
          cfg.atlas='fsl';
        end
      end
      [Y,y,beta,Bcov,STRC, thres, peakxyz] = myspm_graph(xSPM,SPM,hReg, cfg); % this only reads peak
      %% ---something wrong with the cluster name finding...
      %      %% now read cluster name & prob.
      %       % read sigcluster
      %       nii = load_nii([EXP.dir_name '/' fname_sigclus '.img']);
      %       sigvol = round(nii.img);
      %       % find mni-xyz (mm) coordinates
      %       peakijk = xyz2ijk(peakxyz, nii);
      %       clusidx = sigvol(peakijk(1), peakijk(2), peakijk(3));
      %       xyzs = ijk2xyz(find3(sigvol==clusidx), nii);
      %       if strcmp(cfg.atlas,'fsl')
      %         cSTRC = myfsl_atlasquery(xyzs);
      %       elseif strcmp(cfg.atlas,'spm12')
      %         cSTRC = myspm_NMatlas(xyzs);
      %       end
      %       %% ---something wrong with the cluster name finding...
      %       cSTRC.name='na';
      %       cSTRC.prob=0;
      %% ---something wrong with the cluster name finding...
      
      if ~isempty(thres)
        EXP.CorrFDRthres = thres;
      end
      if EXP.print
        if isfield(EXP,'fname_spm_fig')
          spm_print(EXP.fname_spm_fig);
        else
          spm_print
        end
      end
      
      %%
      if isfield(EXP,'mygraph')
        idx_x = find(SPM.xCon(cntrst).c);
        idx_x = idx_x(1);
                % generate a summary table!
                fid = fopen(fname_tab, 'a');
                % 'y-name  x-name  effect  stattype  peak%s\tpeakZ  uncor_pval \t cor_pval(%s) \tK_E\tMNI-x_mm\tMNI-y_mm\tMNI-z_mm\tpeak-Strc-name\npeak-Strc-prob\clus-nStrc-name\nclus-Strc-prob\n'
                fprintf(fid,tab_fmt,...
                  EXP.mygraph.y_name, EXP.mygraph.x_name, beta(idx_x), xSPM.STAT, ...
                  TabDat.dat{PI(ci),9}, TabDat.dat{PI(ci),10},  ...
                  TabDat.dat{PI(ci),11}, TabDat.dat{PI(ci),7}, TabDat.dat{PI(ci),3}, ...
                  TabDat.dat{PI(ci),5}, ...
                  TabDat.dat{PI(ci),12}(1), TabDat.dat{PI(ci),12}(2), TabDat.dat{PI(ci),12}(3), ...
                  STRC.strc.name, STRC.strc.prob); %, ...
                  %cSTRC.name, cSTRC.prob);
                fclose(fid);
%         data2write={EXP.mygraph.y_name, EXP.mygraph.x_name, beta(idx_x), xSPM.STAT, ...
%           TabDat.dat{PI(ci),9}, TabDat.dat{PI(ci),10},  ...
%           TabDat.dat{PI(ci),11}, TabDat.dat{PI(ci),7}, TabDat.dat{PI(ci),3}, ...
%           TabDat.dat{PI(ci),5}, ...
%           TabDat.dat{PI(ci),12}(1), TabDat.dat{PI(ci),12}(2), TabDat.dat{PI(ci),12}(3), ...
%           STRC.strc.name, STRC.strc.prob};
%         xlswrite(fname_tab, data2write, 1, [char(TotalNC+ci+65),'1']);
      end
    end
    TotalNC=TotalNC+NC;
  end
end

%% copy 'significant' results to dir_sum
if TotalNC && isfield(EXP,'dir_sum')
  [~,name1,~]=fileparts(EXP.dir_name);
  [~,~]=mkdir(EXP.dir_sum);
  today=datestr(now,'yyyymmmdd');
  src=fullfile(EXP.dir_name,['spm_',today,'.ps']);
  trg=[EXP.dir_sum,'/spm_',today,'_',name1,'.ps'];
  system(['cp ',src,' ',trg]);
  
  fname_sumtab = [EXP.dir_sum,'/summary.csv'];
  if ~exist(fname_sumtab,'file')
    fid = fopen(fname_sumtab,'w');
    fprintf(fid,hdr_fmt, EXP.thresh.desc);
    fclose(fid);
  end
  system(['tail -n +2 ',fname_tab,' > /tmp/tab2.csv']);
  system(['cat ',fname_sumtab,' /tmp/tab2.csv  >/tmp/sumtab.csv']);
  system(['mv /tmp/sumtab.csv ',fname_sumtab]);
end
end

%% SUB-FUNCTIONS

function xyz = ijk2xyz(ijk, nii)
% xyz = ijk2xyz(ijk, nii)
% converts world-to-voxel coordinates.
%
% Inputs:
%   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
% Output:
%   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
%   nii   the nii structure read using load_untouch_nii.m (or the filename)
%
% see xyz2ijk.m
% (cc) sgKIM, 2014. solleo@gmail.com


if ischar(nii)
  nii = load_untouch_nii(nii);
end

T=[nii.hdr.hist.srow_x; nii.hdr.hist.srow_y; nii.hdr.hist.srow_z; 0 0 0 1];

if numel(ijk) == 3
  if size(ijk,1) > size(ijk,2)
    ijk=ijk';
  end
else
  if size(ijk,1) < size(ijk,2)
    ijk=ijk';
  end
end

xyz = (T)*[(ijk-1) ones(size(ijk,1),1) ]';
xyz(4,:)=[];
xyz=xyz';

end


function ijk = xyz2ijk(xyz, nii)
%  ijk = xyz2ijk(xyz, nii)
% converts world-to-voxel coordinates.
%
% Inputs:
%   xyz   [Nx3 vector] is world-coordinates (e.g. MNI-coord)
%   nii   the nii structure read using load_untouch_nii.m (or the filename)
% Output:
%   ijk   [Nx3 vector] is 1-based voxel coordinates for MATLAB
%
% see ijk2xyz.m
% (cc) sgKIM, 2014. solleo@gmail.com


if ischar(nii)
  hdr= load_untouch_header_only(nii);
else
  hdr = nii.hdr;
end

if numel(xyz) == 3
  if size(xyz,1) > size(xyz,2)
    xyz=xyz';
  end
else
  if size(xyz,1) < size(xyz,2)
    xyz=xyz';
  end
end

T=[hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; 0 0 0 1];
ijk = inv(T) * [xyz ones(size(xyz,1),1)]';
ijk(4,:)=[];
ijk = ijk'+1;


end