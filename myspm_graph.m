function [Y,y,beta,Bcov,STRC,thres,xyz] = myspm_graph(xSPM,SPM,hReg, cfg)
% This is a highly mutated version of spm_graph to work with myspm_result.m
%
% (cc?)

%
% Graphical display of adjusted data
%
% FORMAT [Y y beta Bcov,strc] = myspm_graph(xSPM,SPM,hReg, cfg)
%
% xSPM   - structure containing SPM, distributional & filtering details
%          about the excursion set
% SPM    - structure containing generic details about the analysis
% hReg   - handle of MIP register or [x y z] coordinates
% Ic     - contrast index
% xXi    - explanatory variable index in the desginmatrix
%
%
% Y      - fitted   data for the selected voxel
% y      - adjusted data for the selected voxel
% beta   - parameter estimates (ML or MAP)
% Bcov   - Covariance of parameter estimates (ML or conditional)
% xyz    - peak coordinate in MNI-space (mm)
%
% See spm_getSPM for details.
%__________________________________________________________________________
%
% spm_graph is a Callback script that uses the structures above to:  (i)
% send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
%
% a) Contrasts of parameter estimates (e.g. activations) and their
% standard error.
%
% b) Fitted and adjusted responses that can be plotted against time, scan,
% or an indicator variable in the design matrix.
%
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses that would have been seen in the absence of other
% effects. The PSTH (peristimulus-time histogram) option provides a finite
% impulse response (FIR) estimate of the trial-specific evoked response as
% a function of peristimulus time.  This is estimated by refitting a
% convolution model to the selected voxel using an FIR basis set.  This is
% simply a set of small boxes covering successive time bins after trial
% onset.  The width of each bin is usually the TR.  This option provides a
% more time-resolved quantitative characterisation of the evoked
% hemodynamic response.  However, it should not be over-interpreted because
% inference is usually made using a simpler and more efficient basis set
% (e.g., canonical hrf, or canonical plus time derivative).
%
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g. in
% constructing explanatory variables such as in a psychophysiological
% interaction). To remove or correct for specific effects, specify an
% appropriate F contrast and simply plot the fitted (and adjusted)
% responses after selecting that F contrast. The vectors Y (fitted) and y
% (adjusted) in the workspace will now be corrected for the effects in the
% reduced design matrix (X0) specified in the contrast manager with the
% column indices (iX0) of the confounds in this adjustment.
%
% Plotting data:
% All data and graphics use filtered/whitened data and residuals. In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrasts previously defined can be plotted. This
% ensures that the parameters plotted are meaningful even when there is
% collinearity among the design matrix subpartitions.
%
% Selecting contrasts used for PPMs will automatically give plots
% based on conditonal estimates.
%
% The structure     contrast.contrast      = cbeta;
%                   contrast.standarderror = SE;
%                   contrast.interval      = 2*CI;
%
% is assigned in base workspace for plots of contrasts and their error.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% from Karl Friston's spm_graph.m
% $Id: spm_graph.m 4262 2011-03-25 13:39:20Z guillaume $

% modified by sgKIM, 2015.

Ic=cfg.Ic;  % index of contrast
xXi=cfg.xXi; % column # of the variable of interest in the design matrix


%-Get Graphics figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


%-Delete previous axis and their pagination controls (if any)
%--------------------------------------------------------------------------
spm_results_ui('Clear',Fgraph,2);


%-Find nearest voxel [Euclidean distance] in point list & update GUI
%--------------------------------------------------------------------------
if isempty(xSPM.XYZmm)
  spm('alert!','No suprathreshold voxels!',mfilename,0);
  Y = []; y = []; beta = []; Bcov = [];
  return
end

if numel(hReg) == 1
  xyz = spm_XYZreg('GetCoords',hReg);
else
  xyz = hReg;
end
[xyz,i] = spm_XYZreg('NearestXYZ',xyz,xSPM.XYZmm);
if numel(hReg) == 1, spm_XYZreg('SetCoords',xyz,hReg); end
XYZ     = xSPM.XYZ(:,i);

Cplot = 'Fitted responses';

%-------------------------------
% Getting structure names
%-------------------------------
if ~isfield(cfg,'atlas')
  EXP.atlas='fsl';
end
if strcmpi(cfg.atlas,'spm12')
  [strc] = myspm_NMatlas(xyz);
  strcname = [strc.name];
  STRC=[]; STRC.strc=strc; STRC.strc_all=strc;
elseif strcmpi(cfg.atlas,'fsl')
  [strc,strc_all] = myfsl_atlasquery(xyz);
  strcname = [strc.name,' (',num2str(round(strc.prob)),'%)'];
  STRC=[]; STRC.strc=strc; STRC.strc_all=strc_all;
end

% determine which contrast
%------------------------------------------------------------------
%Ic    = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});


%TITLE = {strcname SPM.xCon(Ic).name};
[~,a]=find(xSPM.XYZ(1,:)== XYZ(1));
[~,b]=find(xSPM.XYZ(2,:)== XYZ(2));
[~,c]=find(xSPM.XYZ(3,:)== XYZ(3));
idx=intersect(intersect(a,b),c);
if numel(idx) == 1
  if strcmp(xSPM.STAT,'T')
    dftxt=['(',num2str(xSPM.df(end)),')'];
  elseif strcmp(xSPM.STAT,'F')
    dftxt=['(',num2str(xSPM.df(1)),',',num2str(xSPM.df(end)),')='];
  end
  peakStat=['Peak ',xSPM.STAT,dftxt,'= ',num2str(xSPM.Z(idx))];
end
TITLE = {strcname,peakStat};

% if xSPM.STAT == 'P'
%   TITLE = {Cplot SPM.xCon(Ic).name '(conditional estimates)'};
% end
% if ~exist('str','var')
%   str=[];
% end
% if isfield(str,'title')
%   TITLE=str.title;
% end

spm('Pointer','Watch');

%-Extract filtered and whitened data from files
%==========================================================================
try
  y = spm_get_data(SPM.xY.VY, XYZ);
catch
  try
    % remap files in SPM.xY.P if SPM.xY.VY is no longer valid
    %------------------------------------------------------------------
    SPM.xY.VY = spm_vol(SPM.xY.P);
    y = spm_get_data(SPM.xY.VY, XYZ);
  end
end

if ~isempty(y), y = spm_filter(SPM.xX.K,SPM.xX.W*y); end

XYZstr = sprintf(' at [%g, %g, %g]',xyz);

%-Compute residuals
%--------------------------------------------------------------------------
if isempty(y)
  
  % make R = NaN so it will not be plotted
  %----------------------------------------------------------------------
  R   = NaN(size(SPM.xX.X,1),1);
  
else
  
  % residuals (non-whitened)
  %----------------------------------------------------------------------
  R   = spm_sp('r',SPM.xX.xKXs,y);
  
end

%-Get parameter and hyperparameter estimates
%==========================================================================
if xSPM.STAT ~= 'P'
  
  %-Parameter estimates:   beta = xX.pKX*xX.K*y;
  %-Residual mean square: ResMS = sum(R.^2)/xX.trRV
  %----------------------------------------------------------------------
  beta  = spm_get_data(SPM.Vbeta, XYZ);
  ResMS = spm_get_data(SPM.VResMS,XYZ);
  Bcov  = ResMS*SPM.xX.Bcov;
  
else
  % or conditional estimates with
  % Cov(b|y) through Taylor approximation
  %----------------------------------------------------------------------
  beta  = spm_get_data(SPM.VCbeta, XYZ);
  
  if isfield(SPM.PPM,'VB');
    % Get approximate posterior covariance at ic
    % using Taylor-series approximation
    
    % Get posterior SD beta's
    Nk=size(SPM.xX.X,2);
    for k=1:Nk,
      sd_beta(k,:) = spm_get_data(SPM.VPsd(k),XYZ);
    end
    
    % Get AR coefficients
    nsess=length(SPM.Sess);
    for ss=1:nsess,
      for p=1:SPM.PPM.AR_P
        Sess(ss).a(p,:) = spm_get_data(SPM.PPM.Sess(ss).VAR(p),XYZ);
      end
      % Get noise SD
      Sess(ss).lambda = spm_get_data(SPM.PPM.Sess(ss).VHp,XYZ);
    end
    
    % Which block are we in ?
    % this needs updating s.t xSPM contains labels of selected voxels
    v = find((SPM.xVol.XYZ(1,:)==XYZ(1))&(SPM.xVol.XYZ(2,:)==XYZ(2))&(SPM.xVol.XYZ(3,:)==XYZ(3)));
    block_index = SPM.xVol.labels(v);
    Bcov=zeros(Nk,Nk);
    for ss=1:nsess,
      % Reconstuct approximation to voxel wise correlation matrix
      post_R=SPM.PPM.Sess(ss).block(block_index).mean.R;
      if SPM.PPM.AR_P > 0
        dh=Sess(ss).a(:,1)'-SPM.PPM.Sess(ss).block(block_index).mean.a;
      else
        dh=[];
      end
      dh=[dh Sess(ss).lambda(1)-SPM.PPM.Sess(ss).block(block_index).mean.lambda];
      for i=1:length(dh),
        post_R=post_R+SPM.PPM.Sess(ss).block(block_index).mean.dR(:,:,i)*dh(i);
      end
      % Get indexes of regressors specific to this session
      scol=SPM.Sess(ss).col;
      mean_col_index=SPM.Sess(nsess).col(end)+ss;
      scol=[scol mean_col_index];
      
      % Reconstuct approximation to voxel wise covariance matrix
      Bcov(scol,scol) = Bcov(scol,scol) + (sd_beta(scol,1)*sd_beta(scol,1)').*post_R;
    end
    
  else
    Bcov     = SPM.PPM.Cby;
    for j = 1:length(SPM.PPM.l)
      l    = spm_get_data(SPM.VHp(j),XYZ);
      Bcov = Bcov + SPM.PPM.dC{j}*(l - SPM.PPM.l(j));
    end
  end
end
CI    = 1.6449;                 % = spm_invNcdf(1 - 0.05);

spm('Pointer','Arrow');
%-Plot
%==========================================================================

%-Colour specifications and index;
%--------------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 .5 .5];
switch Cplot
  
  %-Plot parameter estimates
  %======================================================================
  case 'Contrast estimates and 90% C.I.'
    
    % compute contrast of parameter estimates and 90% C.I.
    %------------------------------------------------------------------
    cbeta = SPM.xCon(Ic).c'*beta;
    SE    = sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));
    CI    = CI*SE;
    
    contrast.contrast      = cbeta;
    contrast.standarderror = SE;
    contrast.interval      = 2*CI;
    assignin('base','contrast',contrast)
    
    % bar chart
    %------------------------------------------------------------------
    figure(Fgraph)
    subplot(2,1,2)
    cla
    hold on
    
    % estimates
    %------------------------------------------------------------------
    h     = bar(cbeta);
    set(h,'FaceColor',Col(2,:))
    
    % standard error
    %------------------------------------------------------------------
    for j = 1:length(cbeta)
      line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
        'LineWidth',6,'Color',Col(3,:))
    end
    
    title(TITLE,'FontSize',12)
    xlabel('contrast')
    ylabel(['contrast estimate',XYZstr])
    set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
    hold off
    
    % set Y to empty so outputs are assigned
    %------------------------------------------------------------------
    Y = [];
    
    
    %-All fitted effects or selected effects
    %======================================================================
  case 'Fitted responses'
    
    % predicted or adjusted response
    %------------------------------------------------------------------
    %        str   = 'predicted or adjusted response?';
    %        if spm_input(str,'!+1','b',{'predicted','adjusted'},[1 0]);
    
    % fitted (predicted) data (Y = X1*c*c^-1*beta)
    %--------------------------------------------------------------
    Y = SPM.xX.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
    %         else
    %
    %             % fitted (corrected)  data (Y = X1o*beta)
    %             %--------------------------------------------------------------
    %             Y = spm_FcUtil('Yc',SPM.xCon(Ic),SPM.xX.xKXs,beta);
    %
    %         end
    
    % adjusted data
    %------------------------------------------------------------------
    y     = Y + R;
    
    % get ordinates
    %------------------------------------------------------------------
    %         Xplot = {'an explanatory variable',...
    %                  'scan or time',...
    %                  'a user specified ordinate'};
    %         Cx    = spm_input('plot against','!+1','m',Xplot);
    %Cx=1;
    % an explanatory variable
    %------------------------------------------------------------------
    %        if     Cx == 1
    
    %             str  = 'Which explanatory variable?';
    %             i    = spm_input(str,'!+1','m',SPM.xX.name);
    i    = xXi;
    x    = SPM.xX.xKXs.X(:,xXi);
    XLAB = SPM.xX.name{i};
    
    %offset adjustment (by sgKIM)
    %meany = mean(spm_get_data(SPM.xY.P, XYZ));
    yoffset = beta(1);
    %     if numel(beta) ~= 1 % But I should NOT do this if the model was 1
    y0 = y + yoffset;
    Y1 = Y + yoffset;
    %     else
    %       y0 = y ;
    %       Y1 = Y ; %??
    %     end
    try xoffset = mean(evalin('base','EXP.vi.val'));
    catch ME
      xoffset=0;
    end
    
    x = x + xoffset;
    
    % plot
    %------------------------------------------------------------------
    figure(Fgraph)
    subplot(2,1,2)
    cla
    hold on
    [p q] = sort(x);
    %Col(4,:)=[.5 .5 1];
    h=[];
    if all(diff(x(q))) % no duplication of x (thus likely continuous..?)
      h(2)=plot(x(q),Y1(q),'LineWidth',4,'Color',Col(2,:));
      %plot(x(q),y(q),'-','Color',[.8 .8 .8]);
      %h(2)=plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(3,:));
      h(1)=plot(x(q),y0(q),'.','MarkerSize',8, 'Color',Col(3,:)); %offset adjustment (by sgKIM)
      plot(x(q),y0(q),'-','Color',[.8 .8 .8]);
    else % for discrete values
      try h(2)=plot(x(q),Y1(q),'.','MarkerSize',8,'Color',Col(1,:));
      catch ME
        hh=get(gca,'children');
        h(2)=hh(end);
      end
      %h(2)=plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
      try h(1)=plot(x(q),y0(q),'o','MarkerSize',4, 'Color',Col(3,:)); %offset adjustment (by sgKIM)
      catch ME
        hh=get(gca,'children');
        h(1)=hh(end);
      end
      xlim = get(gca,'XLim');
      xlim = [-1 1]*diff(xlim)/4 + xlim;
      set(gca,'XLim',xlim)
    end
    if isfield(cfg,'markCorrThres')
      % ref: http://en.wikipedia.org/wiki/Bonferroni_correction
      % Z=atanh(R) ~ N(0,1/sqrt(n-3)) where n is # of samples
      % As we compute the correlation from 410 volumes...
      
      % This is bonferoni correction:
      %       alpha=cfg.markCorrThres.alpha;
      %       NumVox=size(SPM.xVol.XYZ,2);
      %       NumFrames=cfg.markCorrThres.NumFrames;
      %       thres=tanh(norminv(1-(alpha/NumVox/2),0,1/sqrt(NumFrames-3)));
      %       CA=axis;
      %       h1=line([CA(1:2); CA(1:2)]', [thres thres; -thres -thres]', 'color',[.8 1 .8]);
      if isfield(cfg.markCorrThres, 'CorrFDRThres')
        thres=cfg.markCorrThres.CorrFDRThres;
      else
        % How about FDR correction?
        alphalevel=cfg.markCorrThres.alpha;
        %NumVox=size(SPM.xVol.XYZ,2);
        NumFrames=cfg.markCorrThres.NumFrames;
        SE=1/sqrt(NumFrames-3);
        % that means, we set a threshold for each individual
        NumSubj=size(SPM.xY.P,1);
        thres=zeros(NumSubj,1);
        for subji=1:NumSubj
          r=sort(spm_get_data(SPM.xY.VY(subji), SPM.xVol.XYZ));
          z=atanh(r);
          pvals = min(normcdf(z,0,SE), normcdf(-z,0,SE))*2;
          thres(subji) = r(find(pvals > fdr(pvals,alphalevel),1,'last'));
        end
        thres=mean(thres);
      end
      CA=axis;
      h_alpha=line([CA(1:2); CA(1:2)]', [thres thres; -thres -thres]', 'color',[.2 .8 .2]);
    else
      thres=[];
    end
    
    title(TITLE,'FontSize',12)
    
    if isfield(cfg,'x_name')
      xlabel(cfg.x_name);
    else
      xlabel(XLAB)
    end
    if isfield(cfg,'y_name')
      ylabel([cfg.y_name,XYZstr])
    else
      ylabel(['response',XYZstr])
    end
    [~,a]=min(x);
    [~,b]=max(x);
    if y(a) < y(b)
      LOC='southeast';
    else
      LOC='northeast';
    end
    if exist('h_alpha','var')
      legend([h h_alpha'], {'observed','fitted',['q=',sprintf('%0.2f',cfg.markCorrThres.alpha)]},'location',LOC);
    else
      legend(h, {'observed','fitted'},'location',LOC);
    end
    hold off
    grid on; box on;
    %%
    
end


% Turn hold button off - this will alert the user to press it again
%--------------------------------------------------------------------------
try
  set(get(gcbo,'Userdata'),'Value',0);
catch
end


%-call Plot UI
%--------------------------------------------------------------------------
spm_results_ui('PlotUi',gca)


