function myspm_surf (EXP)
% myspm_surf (EXP)
%
% EXP requires:
%  .fname
% (.fstemplate)
% (.dir_fig)
%
% projects spm T-maps onto "fsaverage" surfaec, prints out lateral+medial
% views with+without parcellation boundaries
%
% (cc) 2015, sgKIM

if ~isfield(EXP,'fstemplate')
  fstemplate='fsaverage';
else
  fstemplate=EXP.fstemplate;
end

fname = EXP.fname;
if ~strcmp(fname(1),'/')
  fname=[pwd,'/',fname];
end
if isfield(EXP,'dir_fig')
  dir_fig = EXP.dir_fig;
  [~,~]=mkdir(dir_fig);
else
  dir_fig = pwd;
end
if ~isfield(EXP,'fig_prefix')
  fig_prefix='';
else
  fig_prefix=EXP.fig_prefix;
end

if ~isfield(EXP,'fminmax')
  fminmax=[3 5];
else
  fminmax=EXP.fminmax;
end
fmin=num2str(fminmax(1));
fmax=num2str(fminmax(2));

[path1,name1,~]=fileparts_gz(fname);
cd(path1);
fname2=[path1,'/',name1,'_thres.nii'];
myunix(['FSLOUTPUTTYPE=NIFTI; fslmaths sigclus_+1.nii -add sigclus_-1.nii -bin -mul ',fname,' ',fname2],1)
fname = fname2;

hemi={'lh','rh'};
for s=1:2
  fname_out=[path1,'/',hemi{s},'.',name1,'.mgz'];
  myunix(['mri_vol2surf --mov ',fname,' --mni152reg --trgsubject ',fstemplate,' ', ...
    ' --hemi ',hemi{s},' --o ',fname_out],1);
  suffix={'','.aparc'};
  for a=1:2
    fname_tcl=[tempname,'.tcl'];
    fid=fopen(fname_tcl,'w');
    % find names for lateral/medial
    fname_tiff1=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.lat.tiff'];
    fname_tiff2=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.med.tiff'];
    % and load parcellations
    aparcarg={'','labl_import_annotation aparc.annot \n set labelstyle 1\n '};
    fprintf(fid,[aparcarg{a},' set colscalebarflag 1\n redraw\n',...
      ' save_tiff %s\n rotate_brain_y 180\n redraw\n save_tiff %s\n exit\n'],...
      fname_tiff1,fname_tiff2);
    fclose(fid);
    myunix(['tksurfer ',fstemplate,' ',hemi{s},' inflated -gray ', ...
      ' -fminmax ',fmin,' ',fmax,' -overlay ',fname_out,' -tcl ',fname_tcl],1);
  end
end

end