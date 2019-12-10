function job = myspm_viewsurf_glms(job)
% job = myspm_viewsurf_glms(job)
% job requires:
%  .dir_glm
%  .fn_bigfig
%  .prefix      e.g. 'spmT' or 'con'
% (.titleprepfix)
% (cc) 2019, sgKIM.

% NOFIGURE = 0; % For debugging
NOFIGURE = 1; % For normal usage

cd(job.dir_glm)
load('SPM.mat','SPM')
if ~isfield(job,'layout')
  job.layout = '1x4'
end
if ~isfield(job,'fn_bigfig')
  job.fn_bigfig = [pwd,'/bigfig.png']
end

% %% RUNNING all projections in parallel first?
% files = {};
% if isfield(job,'mean')
%   % FIND constant term:
%   if ~isfield(SPM,'Sess')
%     meanimage = 'beta_0001';
%   else
%     error('How can I find it?')
%     meanimage = '???';
%   end
%   files = [files [meanimage,'.nii']];
% end
% if ~isfield(job,'idxCntrst')
%   job.idxCntrst = 1:2:numel(SPM.xCon);
% end
% for i = job.idxCntrst
%   fn_src = sprintf('%s_%04i.nii',job.prefix,i);
%   if ~exist(fn_src,'file')
%     continue
%   end
%   files = [files fn_src];
%   if isfield(job,'sigclus') && job.sigclus
%     sign = {'+','-'};
%     for jSign = 1:2
%       fn_sigclus = ['sigclus_',sign{jSign},SPM.xCon(i).name(2:end),'.nii'];
%       if exist(fn_sigclus,'file')
%         files = [files fn_sigclus];
%       end
%     end
%   end
% end
% files2run = {};
% for i = numel(files):-1:1
%   if ~exist([files{i}(1:end-4),'.fsaverage6.mat'],'file')
%     files2run = [files2run files{i}];
%   end
% end
% if ~isempty(files2run)
% %   delete(gcp('nocreate'))
% %   parpool(min([16 numel(files)]))
% %   par
%   for i = 1:numel(files)
%     [~, ~, ~] ...
%       = fsss_mni_to_fsavg(files{i}, ...
%       'fsaverage6', struct('nofigure',NOFIGURE,'interpopt','nearest',...
%       'mapmethod','nnf'));
%   end
% %   delete(gcp('nocreate'))
% end

%%
fname_png = {};
%% Mean Image
if isfield(job,'mean')
  % FIND constant term:
  if ~isfield(SPM,'Sess')
    meanimage = 'beta_0001';
  else
    error('How can I find it?')
    meanimage = '???';
  end
  [~, surfs, Y] ...
    = fsss_mni_to_fsavg([meanimage,'.nii'], ...
      'fsaverage6', struct('nofigure',NOFIGURE,'interpopt','nearest',...
      'mapmethod','nnf'));
  cfg = [];
  cfg.fname_png = [tempname,'.png'];
  fname_png = [fname_png cfg.fname_png];
  cfg.colorbartitle = ['Mean ',job.measurename];
  Y{1}(Y{1}==0) = nan;  Y{2}(Y{2}==0) = nan;
  cfg.layout = job.layout;
  if isfield(job,'masks')
    cfg.masks = job.masks';
  end
  fsss_view(surfs, Y, cfg)
end
%%
if ~isfield(job,'idxCntrst')
  job.idxCntrst = 1:2:numel(SPM.xCon);
end
for i = job.idxCntrst
  fn_src = sprintf('%s_%04i.nii',job.prefix,i);
  if ~exist(fn_src,'file')
    warning('%s not found. skipping it.',fn_src)
    continue
  end
  [~, surfs, Y] = fsss_mni_to_fsavg(fn_src, ...
    'fsaverage6', struct('nofigure',NOFIGURE,'interpopt','nearest',...
      'mapmethod','nnf'));
  
  cfg = [];
  cfg.fname_png = [tempname,'.png'];
  fname_png = [fname_png cfg.fname_png];
  cfg.colorbartitle = ['Effect of ',SPM.xCon(i).name];
  if contains(cfg.colorbartitle,'All Sessions')
    cfg.colorbartitle(end-14:end) = [];
  end
  if contains(lower(SPM.xCon(i).name),'effect')
    cfg.colorbartitle = SPM.xCon(i).name;
  end
  if isfield(job,'measurename')
    cfg.colorbartitle = [cfg.colorbartitle,' on ',job.measurename];
  end
  if isfield(job,'titleprepfix')
    cfg.colorbartitle = [job.titleprepfix,' ',cfg.colorbartitle];
  end
  
  if isfield(job,'caxis')
    cfg.caxis = job.caxis;
  end
  if ~isfield(job,'thresp')
    job.thresp = 0.05;
  end
  if strcmp(SPM.xCon(i).STAT,'T')
    df = SPM.xX.erdf;
    %     cfg.colorbarxlabel = sprintf('T(%d) | unc-p<%.3f', round(df), job.thresp);
    cfg.colorbarxlabel = ['T(',num2str(round(df)),') | unc-p<',...
      num2str(job.thresp)];
    cfg.thres = icdf('t',1-job.thresp,df);
    if job.thresp == 1
      cfg.thres = 0;
    end
    if ~isfield(job,'caxis')
      cfg.caxis = [-10 10];
    else
      cfg.caxis = job.caxis;
    end
  else
    df = [SPM.xCon(i).eidf SPM.xX.erdf];
    cfg.thres = icdf('f', 1-job.thresp, df(1), df(2));
    %     cfg.colorbarxlabel = sprintf('F(%d,%d) | unc-p<%.3f',round(df), job.thresp);
    if job.thresp == 1
      cfg.thres = 0;
    end
    cfg.colorbarxlabel = ['F(',num2str(round(df(1))),...
      ',',num2str(round(df(2))),' | unc-p<',num2str(job.thresp)];
    if ~isfield(job,'caxis')
      cfg.caxis = [0 20];
    else
      cfg.caxis = job.caxis;
    end
  end
  
  if isfield(job,'colorbarxlabelsuffix')
    cfg.colorbarxlabel = [cfg.colorbarxlabel job.colorbarxlabelsuffix];
  end
  
  Y{1}(Y{1}==0) = nan; Y{2}(Y{2}==0) = nan;
  cfg.layout = job.layout;
  if isfield(job,'masks')
    cfg.masks = job.masks';
  end
  
  if isfield(job,'sigclus') && job.sigclus
    sign = {'+','-'};
    sigclus = cell(1,2);
    for jSign = 1:2
      fn_sigclus = ['sigclus_',sign{jSign},SPM.xCon(i).name(2:end),'.nii'];
      if exist(fn_sigclus,'file')
        [~, surfs, sigclus{jSign}] = fsss_mni_to_fsavg(fn_sigclus, ...
          'fsaverage6',struct('nofigure',NOFIGURE,'interpopt','nearest',...
          'mapmethod','nnf'));
        sigclus{jSign} = sigclus{jSign};
      else
        warning('%s not found',fn_sigclus)
        sigclus{jSign} = {false(size(Y{1})), false(size(Y{2}))};
      end
    end
    sigclus = ~~(sigclus{1} + sigclus{2});
    if ~isfield(cfg,'masks')
      cfg.masks = sigclus;
    else
      cfg.masks = ~~(reshape(cfg.masks,1,[]) .* reshape(sigclus,1,[]));
    end
  end
  
  fsss_view(surfs, Y, cfg)
end
imageconcat(fname_png, job.fn_bigfig, 1, struct('scale',0.25))

end