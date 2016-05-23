function myspm_surf (EXP)
% myspm_surf (EXP)
%
% EXP requires:
%  .fname_tmap
%  .fname_clustermap
% (.fstemplate)
% (.dir_fig)
% (.aparc) additional map with automatic parcellation boundaries: TRUE or FALSE
% (.proj) projection method: "max", "avg", "mid"
%
% projects T-maps (.fname_tmap) onto "fsaverage" surfaces, prints out lateral+medial views
% masked by cluster ID map (.fname_clustermap)
%
% (cc) 2015, sgKIM

if ~exist('EXP','var')
  EXP=[];
end
if ~isfield(EXP,'fstemplate')
  fstemplate='fsaverage';
else
  fstemplate=EXP.fstemplate;
end

if ~isfield(EXP,'fname_tmap')
  EXP.fname_tmap='spmT_0001.img';
end
fname = EXP.fname_tmap;
if ~strcmp(fname(1),'/'), 
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

if ~isfield(EXP,'projarg')
  EXP.projarg='max';
end
switch EXP.projarg
  case 'max'
    projarg=' --projfrac-max -1 2 .1 ';
  case 'avg'
    projarg=' --projfrac-avg .2 .8 .2 ';
  case 'mid'
    projarg=' --projfrac 0.5 ';
  otherwise
    projarg='';
end

if isfield(EXP,'aparc')&&EXP.aparc
  A=2;
else
  A=1;
end

[path1,name1,~]=fileparts_gz(fname);
cd(path1);
fname2=[path1,'/',name1,'_thres.nii'];
[path1,~,~]=fileparts_gz(fname);
if isfield(EXP,'fname_clustermap')
  fnames = EXP.fname_clustermap;
  if iscell(fnames) && (numel(fnames) == 2)
    query1=fnames{1}; query2=fnames{2};
  else
    query1=fnames; query2=tempname;
  end
else
  query1='sigclus_+1.*'; query2='sigclus_-1.*';
end
[~,res1]= mydir([path1,'/',query1]);
[~,res2]= mydir([path1,'/',query2]);
if isempty(res1) || isempty(res2)
if isempty(res1) && ~isempty(res2)
  res1=res2;
elseif ~isempty(res1) && isempty(res2)
  res2=res1;
else
  error('No significant cluster found!');
end
end

unix(['FSLOUTPUTTYPE=NIFTI; fslmaths ',res1{end},' -add ',res2{end}, ...
  ' -bin -mul ',fname,' ',fname2])
fname = fname2;

hemi={'lh','rh'};
for s=1:2
  fname_out=[path1,'/',hemi{s},'.',name1,'.mgz'];
  unix(['mri_vol2surf --mov ',fname,' --mni152reg --trgsubject ', ...
    fstemplate,' --hemi ',hemi{s},' --o ',fname_out,' ',projarg]);
  suffix={'','.aparc'};
  for a=1:A
    fname_tcl='/tmp/cmd.tcl';
    fid=fopen(fname_tcl,'w');
    % find names for lateral/medial
    fname_tiff1=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.view1.tiff'];
    fname_tiff2=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.view2.tiff'];
    fname_tiff3=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.view3.tiff'];
    fname_tiff4=[dir_fig,'/',fig_prefix,name1,suffix{a},'.',hemi{s},'.view4.tiff'];
    % and load parcellations
    aparcarg={'','labl_import_annotation aparc.annot \n set labelstyle 1\n '};
    %' resize_window 1200\n scale_brain 2\n', ... DIDN'T WORK!!
    fprintf(fid,[aparcarg{a},' set colscalebarflag 1\n ', ...
      ' redraw\n save_tiff %s\n', ...
      ' rotate_brain_y 90 \n redraw\n save_tiff %s\n', ...
      ' rotate_brain_y 180\n redraw\n save_tiff %s\n', ...
      ' rotate_brain_y -90\n redraw\n save_tiff %s\n exit\n'],...
      fname_tiff1,fname_tiff2, fname_tiff3, fname_tiff4);
    fclose(fid);
    unix(['tksurfer ',fstemplate,' ',hemi{s},' inflated -gray ', ...
      ' -fminmax ',fmin,' ',fmax,' -overlay ',fname_out,' -tcl ',fname_tcl]);
  end
end
end

function [path,name,ext]=fileparts_gz(filename)
% [path,name,ext]=fileparts_gz(filename)
% (cc) 2014, sgKIM.

[path,name,ext] = fileparts(filename);
if strcmp(ext,'.gz')
  ext = [name(end-3:end),ext];
  name = name(1:end-4);
end
end

function [name,fullname] = mydir(str0)
% [name,fullname] = mydir(str0, SessSorting)
%
% allows you to use mulitple wildcard queries
% if SessSorting=1, sort the filenames so that "S9" comes before "S10"
% (cc) 2015, sgKIM

try txt = dir(str0);
catch exception
  name=[]; fullname=[];
  return
end
name={};
fullname={};
j=0;
for i=1:numel(txt)
  if ~txt(i).isdir
    j=j+1;
    name{i} = txt(i).name;
    [~,name1,ext1] = fileparts(name{i});
    [path1,~,~] = fileparts_gz(str0);
    if isempty(path1) || strcmp(path1,'./')
      fullname{i}=fullfile(pwd,[name1,ext1]);
    else
      fullname{i}=fullfile(path1,[name1,ext1]);
    end
  end
end

if exist('SessSorting','var')
  for i=1:numel(name)
    idx = findstr(name{i},'_');
    sessNum(i) = str2double(name{i}(2:(idx-1)));
  end
  [~,sortidx] = sort(sessNum);
  name = name(sortidx);
  fullname = fullname(sortidx);
end
end 