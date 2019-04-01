function EXP = myspm_cntrst (EXP)
% EXP = myspm_cntrst (EXP)
%
% creates CON and SPM images
% and calls myspm_result
%
% EXP requires:
%  .dir_glm
%  .cntrstMtx  (for T-contrasts)
%  .titlestr   (for T-contrasts)
%  .effectOfInterest [1x1] 0 | 1 (default=1)
%  .FcntrstMtx (for F-contrasts)
%  .Ftitlestr  (for F-contrasts)
%
% (cc) 2018, sgKIM. solleo@gmail.com


spm('defaults','fmri');
if ~isfield(EXP,'effectOfInterest')
  EXP.effectOfInterest = 1;
end

%% contrast manager
matlabbatch={};
con=[];
con.spmmat = {[EXP.dir_glm,'/SPM.mat']};
load(con.spmmat{1});

if ~isfield(EXP,'NumSess')
  try
    EXP.NumSess=numel(SPM.Sess);
  catch
    EXP.NumSess=1;
  end
end
NumSess = EXP.NumSess;
SIGN ={'+','-'};
if isfield(EXP,'cntrstMtx')
  NumCnt = size(EXP.cntrstMtx,1);
  for k=1:NumCnt
    if isfield(EXP,'titlestr')
      con.consess{k}.tcon.name = EXP.titlestr{k};
    else
      con.consess{k}.tcon.name = ['Contrast#',num2str(k)];
    end
    con.consess{k}.tcon.convec = EXP.cntrstMtx(k,:);
    if NumSess>1
      con.consess{k}.tcon.sessrep = 'repl';
    else
      con.consess{k}.tcon.sessrep = 'none';
    end
    if isfield(EXP,'sessrep')
      con.consess{k}.tcon.sessrep=EXP.sessrep;
    end
  end
  % else
  %   NumCond = EXP.NumCond;
  %   k=0;
  %   % First single regressor (compared to implicit baseline)
  %   for i=1:NumCond
  %     for s=1:2
  %       k=k+1;
  %       if isfield(EXP,'titlestr')
  %         con.consess{k}.tcon.name = EXP.titlestr{k};
  %       else
  %         con.consess{k}.tcon.name = [SIGN{s}, EXP.COND(i).name];
  %       end
  %       con.consess{k}.tcon.name   = EXP.titlestr{k};
  %       con.consess{k}.tcon.convec = zeros(1,NumCond+6+1);
  %       con.consess{k}.tcon.convec(i) = s*-2+3;
  %       if NumSess>1
  %         con.consess{k}.tcon.sessrep = 'repl';
  %       else
  %         con.consess{k}.tcon.sessrep = 'none';
  %       end
  %       if isfield(EXP,'sessrep')
  %         con.consess{k}.tcon.sessrep=EXP.sessrep;
  %       end
  %     end
  %   end
end
%% Effect of interest for fMRI data
if isfield(SPM,'Sess') && EXP.effectOfInterest
  if ~isfield(EXP,'FcntrstMtx')
    EXP.FcntrstMtx={};
    EXP.Ftitlestr={};
  end
  EXP.Ftitlestr = [EXP.Ftitlestr 'Effect of interest'];
  load(con.spmmat{1});
  nRegInt = numel(SPM.Sess(1).U); % number of regressors of interest
  nxBForder = size(SPM.xBF.bf,2); % order of the response function
  nRegNsn = size(SPM.Sess(1).C.C,2); % number of regressors of nuisance
  % SANITY CHECK:
  if (nRegInt*nxBForder+nRegNsn+1)*numel(SPM.Sess) ~= size(SPM.xX.X,2)
    error('# of regressors does not fit! CANNOT set up Effect of Interest contrast!')
  end
  EXP.FcntrstMtx = [EXP.FcntrstMtx ...
    [eye(nRegInt*nxBForder) zeros(nRegInt*nxBForder,nRegNsn)]];
end

%% F-contrast
if isfield(EXP,'FcntrstMtx')
  for j=1:numel(EXP.FcntrstMtx)
    k=k+1;
    con.consess{k}.fcon.name = EXP.Ftitlestr{j};
    con.consess{k}.fcon.convec = EXP.FcntrstMtx{j};
    if NumSess>1
      con.consess{k}.fcon.sessrep = 'repl';
    else
      con.consess{k}.fcon.sessrep = 'none';
    end
    if isfield(EXP,'sessrep')
      con.consess{k}.fcon.sessrep=EXP.sessrep;
    end
  end
end
%%
if isfield(EXP,'newContrast')
  con.delete = EXP.newContrast;
else
  con.delete = 1;
end
if con.delete == 1
  unix(['rm -f ',EXP.dir_glm,'/spm*']);
  unix(['rm -f ',EXP.dir_glm,'/con*']);
  unix(['rm -f ',EXP.dir_glm,'/sigclus*']);
end
matlabbatch{1}.spm.stats.con = con;

spm_jobman('initcfg')
spm_jobman('run', matlabbatch)

%% Now create result reports
EXP.mygraph.y_name = 'y';
EXP.mygraph.x_name = 'x';
if ~isfield(EXP,'thres')
  EXP.thres.desc  = 'cluster';
  EXP.thres.alpha = 0.05;
end
if ~isfield(EXP,'NOREPORT') && ~isfield(EXP,'noreport')
  EXP = myspm_result(EXP);
end


end
