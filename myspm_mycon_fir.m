function EXP = myspm_mycon_fir(EXP)
% EXP = myspm_mycon_fir(EXP)
%
% this does... what?
%
% EXP requires:
%  .dir_glm
%  .cntname_fir
%  .cntvec_fir
%  .num_delay
%  .TR_sec
%  .length_sec
%  .winsize
%  .num_cond
%
% (cc) 2016, sgKIM

num_delay  = EXP.num_delay;
length_sec = EXP.length_sec;
TR_sec  = EXP.TR_sec;
winsize = EXP.winsize;
postonset_sec=[];
cntMtx=[];
cntName={};
for k=1:size(EXP.cntvec_fir,1)
 if ~mod(winsize,2)
  winsize=winsize-1;
  warning(['Window size should be even: changing as:',num2str(winsize)]);
 end
 for j=1:numel(num_delay)
  for t=1:(num_delay-winsize-1)
   idx=t:(t+winsize-1);
   ref_time=median(idx*TR_sec);
   postonset_sec=[postonset_sec ref_time];
   cntName=[cntName, [EXP.cntname_fir{k},'[',pad(ref_time,3),'s]']];
   
   cntVec=[];
   for c=1:EXP.num_cond
    x=zeros(1,num_delay);
    x(idx) = ones(1,numel(idx))*EXP.cntvec_fir(k,c);
    cntVec=[cntVec x];
   end
   cntMtx=[cntMtx; cntVec];
  end
 end
end
hf=figure; imagesc(cntMtx);
drawnow;
screen2png([EXP.dir_glm,'/fir_',num2str(length_sec),'s_', ...
 num2str(num_delay),'-th.png']); close(hf);

EXP2=EXP;
EXP2.conname = cntName;
EXP2.convec  = cntMtx;

EXP = myspm_mycon (EXP2);
EXP.postonset_sec = postonset_sec;
end
