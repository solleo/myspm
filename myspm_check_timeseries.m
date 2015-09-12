function EXP = myspm_check_timeseries (EXP)
% EXP = myspm_check_timeseries (EXP)
%
% EXP
%  .fname_proc
%  .name_proc
%
% (cc) 2015, sgKIM.   solleo@gmail.com   https:ggooo.wordpress.com

subjid = EXP.subjid;

NumProc = numel(EXP.fname_proc);
P=cell(1,NumProc);
y=cell(1,NumProc);
p=1;
P{p} = spm_vol(EXP.fname_proc{p});
xyz=[round(P{1}(1).dim(1)/2), round(P{1}(1).dim(2)/2), round(P{1}(1).dim(3)/(2))]';
y{p} = spm_get_data(P{p},xyz);
for p=2:NumProc
  P{p} = spm_vol(EXP.fname_proc{p});
  y{p} = spm_get_data(P{p},xyz);
end
dim=[P{1}(1).dim numel(P{1})];
numtp = min(100,dim(4));
if dim(4)>=200
  tp=[round(dim(4)/4):round(dim(4)/4)+(numtp-1)];
else
  tp=1:numtp;
end

%ProNames={'Original','Slice timing correction','Unwarping/realignment'};
ProNames = EXP.name_proc;
figure;
set(gcf,'position',[ 841         206        1075         744]);
subplot(5,NumProc,[1:NumProc]); plot([y{1}(tp) y{2}(tp) y{3}(tp)]);
legend(ProNames,'location','EastOutside')
title(subjid,'fontsize',16);
xlabel('TR'); ylabel('Image intensity');

for p=1:NumProc
  subplot(5,NumProc,NumProc+p)
  x=spm_read_vols(P{p}(1));
  imagesc(x(:,:,xyz(3))'); set(gca,'ydir','nor'); axis image; colormap(bone);
  hold on; scatter(xyz(1),xyz(2),'r+');
  axis off; title(ProNames{p})
  set(gca,'
end

% take 3 principle components
for p=1:NumProc
  Y=zeros(dim(4),prod(dim(1:3)));
  for j=1:dim(4)
    Y(j,:)=reshape(spm_read_vols(P{p}(j)),[1,prod(dim(1:3))]);
  end
  coef = pca(Y');
  PC{p} = coef(:,[1:3]);
end

for k=1:3
  subplot(5,NumProc,[1:NumProc]+(2+k-1)*NumProc);
  plot([PC{1}(tp,k), PC{2}(tp,k), PC{3}(tp,k)]);
  legend(ProNames)
  xlabel('TR'); ylabel(['PC',num2str(k)]);
end
end
