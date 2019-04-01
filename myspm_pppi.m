function EXP = myspm_pppi (EXP)
% EXP = myspm_pppi (EXP)
%
% EXP requires:
% .subjid
% .dir_glm
% .fname_voi
% .voi_name
% .convec
% .conname
%
% (cc) sgKIM, 2017.
cd(EXP.dir_glm)
%% Check whether VOI fits into the mask.nii
fname1=mydir('mask*');
mask=load_uns_nii(fname1);
voi =load_uns_nii(EXP.fname_voi);
if sum(double(voi.img(:))) ~= sum(double(voi.img(:)).*double(mask.img(:)))
 error('The VOI mask includes voxels outside the GLM mask! Use VOI mask that fits to interaction of all GLM masks!');
end

%%
P = [];
P.subject   = EXP.subjid;
P.directory = EXP.dir_glm;
P.VOI       = EXP.fname_voi;
P.Region    = EXP.voi_name;%STR{j};
P.estimate  = 1;
P.contrast  = 0;
P.analysis  = 'psy';
P.extract   = 'eig';
P.method    = 'cond';
P.weighted = 0;
P.contrast = 0;
P.Estimate = 1;
PPPI(P)
% delete ResI files (because I modified my SPM12 to leave ResI by me.
dir_ppi = [P.directory,'/PPI_',P.Region];
myunix(['rm -f ',dir_ppi,'/ResI*.nii']);
% move all the log/mat/txt
[~,~]=mkdir([dir_ppi,'/pppi_log/']);
myunix(['mv ',P.directory,'/',P.subject,'_PP* ',dir_ppi,'/pppi_log/']);
myunix(['mv ',P.directory,'/',P.subject,'*',EXP.voi_name,'* ',dir_ppi,'/pppi_log/']);

% compute contrast files (requires .conname and .convec)
if isfield(EXP,'convec') && isfield(EXP,'conname')
EXP2=EXP;
EXP2.dir_glm=dir_ppi;
h1=figure;
myspm_mycon(EXP2)
close(h1)
end

end
