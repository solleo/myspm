function EXP = myspm_cntrst (EXP)
% each contrast list will create spmT_0001 ... spmT_xxxx

spm('defaults','fmri');

%% contrast manager
matlabbatch={};
con=[];
con.spmmat = {[EXP.dir_glm,'/SPM.mat']};

NumSess = EXP.NumSess;
SIGN ={'+','-'};

if isfield(EXP,'cntrstMtx')
  NumCond = size(EXP.cntrstMtx,1);
  for k=1:NumCond
    if isfield(EXP,'titlestr')
      con.consess{k}.tcon.name = EXP.titlestr{k};
    else
      con.consess{k}.tcon.name = ['condition#',num2str(k)];
    end
    con.consess{k}.tcon.convec = EXP.cntrstMtx(k,:);
    if NumSess>1
      con.consess{k}.tcon.sessrep = 'repl';
    else
      con.consess{k}.tcon.sessrep = 'none';
    end
  end
else
  NumCond = EXP.NumCond;
  k=0;
  % First single regressor (compared to implicit baseline)
  for i=1:NumCond
    for s=1:2
      k=k+1;
      if isfield(EXP,'titlestr')
        con.consess{k}.tcon.name = EXP.titlestr{k};
      else
        con.consess{k}.tcon.name = [SIGN{s}, EXP.COND(i).name];
      end
      con.consess{k}.tcon.name = EXP.titlestr{k};
      con.consess{k}.tcon.convec = zeros(1,NumCond+6+1);
      con.consess{k}.tcon.convec(i) = s*-2+3;
      if NumSess>1
        con.consess{k}.tcon.sessrep = 'repl';
      else
        con.consess{k}.tcon.sessrep = 'none';
      end
    end
  end
end

if isfield(EXP,'newContrast')
  con.delete = EXP.newContrast;
else
  con.delete = 1;
end
if con.delete == 1
  unix(['rm -f ',EXP.dir_glm,'/spm*']);
  unix(['rm -f ',EXP.dir_glm,'/con*']);
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
  myspm_result(EXP);
end

end
