function EXP=myspm_create_myreg(EXP)
% creates custom regressors with convolution by HRF
% 
% EXP requires:
% .tr_sec    [1x1] TR in seconds
% .n_scans   [1x1] # of scans (volumes)
% .fname_mat '1xS' e.g. matlabfilename.mat
% .ONSETS_sec {[1 x v1], [1 x v2], ...}
% .duration_sec [1x1]
% .REGNAMES   {'1 x s1', '1 x s2', ...}
%
% (cc) 2017, sgKIM. solleo@gmail.com

% 1. create matrix
onset_tr  = round(ONSETS_sec./EXP.tr_sec);
dur_tr    = round(EXP.duration_sec/EXP.tr_sec);
K         = numel(ONSET_sec);
myregvals = zeros(n_scans,K);
for k=1:numel(K)
 idx =[];
 for j=1:dur_tr; % over the duration
  idx = [idx onset_tr(:,k)+(j-1)];
 end
 X(idx,k) = 1;
 exp1=[];
 exp1.TR_sec = EXP.tr_sec;
 exp1.u = X(:,k);
 
 % and convolute it with HRF
 exp1 = myspm_conv(exp1);
 X(:,k) = exp1.uc;
end

myregvals =X;
myregnames=REGNAMES;
save(fname_mat,'myregnames','myregvals');
end