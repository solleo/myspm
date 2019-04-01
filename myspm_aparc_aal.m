function APARC = myspm_aparc_aal(EXP)
% EXP requires:
%  .fname_epi_mni

%%- resample AAL for different resolution
trg=EXP.fname_epi_mni;
dir1=[spm('dir'),'/toolbox/aal/'];
src='/tmp/aal2.nii';
if ~exist(src,'file')
 copyfile([dir1,'/aal2.nii.gz'],'/tmp/');
 unix('gunzip /tmp/aal2.nii.gz');
end
nii1 = load_untouch_nii([dir1,'/aal2.nii.gz']);

exp1=[];
exp1.interp=0;
exp1.fname_moving=src;
exp1.fname_fixed=trg;
myspm_coreg(exp1);

%%- read AAL parcellation
nii=load_untouch_nii('/tmp/oaal2.nii');
APARC.img = nii.img;
M = max(nii.img(:));  xyz=zeros(M,3);  numvox=zeros(M,1);
for i=1:M
 ind = nii.img == i;
 numvox(i) = sum(ind(:));
 % find centroid from 2mm resolution (highest given)
 ijk = mean(find3(nii1.img == i),1);
 xyz(i,:) = ijk2xyz(ijk, nii1);
end

%%- read structure names
T=readtable([spm('dir'),'/toolbox/aal/aal2.nii.txt']);
label=T.Var1;
strnames=T.Var2;
origlabel=T.Var3;
tab=table(label,strnames, numvox, origlabel, xyz);
APARC.tab = tab;
end


