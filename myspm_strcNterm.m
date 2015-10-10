function EXP = myspm_strcNterm(EXP)

if isfield(EXP,'model')&&~isfield(EXP,'vi')
  % terms to spm.struct
  varnames = char(EXP.model);
  X = double(EXP.model);
  [numSubj,numReg] = size(X);
  % where is the regressor to test?
  if isfield(EXP,'cidx'),   c = EXP.cidx;
  else                      c = numReg;     end
  EXP.vi.name = varnames{c};
  EXP.vi.val  = X(:,c);
  % where is the constant?
  d = find(sum(X)==numSubj);
  % what's left for covariates?
  idx = 1:numReg;
  idx([c,d]) = [];
  if numel(idx)
    for i=1:numel(idx)
      j=idx(i);
      EXP.vn(i).name = varnames{j};
      EXP.vn(i).val  = X(:,j);
    end
  end
elseif ~isfield(EXP,'model')&&isfield(EXP,'vi')
  % spm.struct to terms
  EXP.model = 1 + term(EXP.vi.val, EXP.vi.name);
  EXP.cidx = 2;
  if isfield(EXP,'vn')
    for j=1:numel(EXP.vn)
      EXP.model = EXP.model + term(EXP.vn(j).val, EXP.vn(j).name);
    end
  end
end
  
EXP.model_desc = fsss_model_desc(EXP.model, EXP.cidx);
end
