function JOB = myspm_strcNterm(JOB)
% converts 'model' structure from SurfStat into the SPM's structure
%
% (cc) 2016. sgKIM.  mailto://solleo@gmail.com  https://ggooo.wordpress.com/

if isfield(JOB,'model')&&~isfield(JOB,'vi')
  % terms to spm.struct
  varnames = char(JOB.model);
  X = double(JOB.model);
  [numSubj,numReg] = size(X);
  % where is the regressor to test?
  if isfield(JOB,'cidx')
    c = JOB.cidx;
  else
    c = numReg;
  end
  if isempty(varnames{c})
    JOB.vi.name=['var',num2str(c)];
  else
    JOB.vi.name = varnames{c};
  end
  JOB.vi.val  = X(:,c);
  % where is the constant?
  d = find(sum(X)==numSubj);
  % what's left for covariates?
  idx = 1:numReg;
  idx([c,d]) = [];
  if numel(idx)
    for i=1:numel(idx)
      j=idx(i);
      JOB.vn(i).name = varnames{j};
      JOB.vn(i).val  = X(:,j);
    end
  end
elseif ~isfield(JOB,'model')&&isfield(JOB,'vi')
  % spm.struct to terms
  JOB.model = 1 + term(JOB.vi.val, JOB.vi.name);
  JOB.cidx = 2;
  if isfield(JOB,'vn')
    for j=1:numel(JOB.vn)
      JOB.model = JOB.model + term(JOB.vn(j).val, JOB.vn(j).name);
    end
  end
end
if ~isfield(JOB,'model_desc')
  if isfield(JOB,'cidx')
    JOB.model_desc = fsss_model_desc(JOB.model, JOB.cidx);
  else
    JOB.model_desc = fsss_model_desc(JOB.model);
  end
end
end
