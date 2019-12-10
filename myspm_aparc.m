% function JOB = myspm_aparc(JOB)
% % JOB = myspm_aparc(JOB)
% %
% % create a look-up-table
% %
% % when prefix='o', coregister automatic parcellation (e.g. a2009s)
% % from freesurfer into functional volumetric space
% %
% % (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com
% 
% % 0. set parameters
% path0=pwd;
% subjID = fsss_subjID(JOB.subjID);
% Segname={'aparc.a2009s+aseg','aparc+aseg'};
% Nonbrain={'unknown','CC','Optic-Chiasm','VentralDC','WM-hypointensities','Cerebellum-Cortex', ...
%  'Ventricle','plexus','vessel','White-Matter','Inf-Lat-Vent','CSF','Brain-Stem'};
% if ~isfield(JOB,'aparcs'), JOB.aparcs = [1,2]; end
% if ~isfield(JOB,'prefix'), prefix='o'; else prefix=JOB.prefix; end
% for k=JOB.aparcs
%  for n=1:numel(subjID)
%   subjid = subjID{n};
%   path1=[fullfile(JOB.dir_base,subjid),'/'];
%   segname = Segname{k};
%   fname_aparc = [path1,segname,'.nii'];
%   if ~exist(fname_aparc,'file')
%    unix(['mri_convert --resample_type nearest ', ...
%     fullfile(fullfile(JOB.dir_fs,subjid),'mri',[segname,'.mgz']), ...
%     ' ',fname_aparc]);
%   end
%   
%   if strcmp(prefix,'o')
%    cd(path1)
%    fname_aparc=[prefix,segname,'.nii'];
%    if ~exist(fname_aparc,'file')
%     job1=[];
%     job1.name_fixed  = JOB.name_epi;
%     job1.name_moving = JOB.name_t1w; % moving name must be identical to the last one
%     job1.name_others = {[path1,segname,'.nii']};
%     job1.interp       = 0;
%     myspm_coreg(job1);
%    end
%   end
%   
%   % simplifying lables for mygraph_adjmtx_density
%   fname_aparc=[prefix,segname,'.nii'];
%   nii = load_uns_nii(fname_aparc);
%   nii1=nii;
%   nii1.img = nii1.img*0;
%   aparc = single(nii.img(:));
%   label_fs = unique(aparc(:));
%   nvox_fs = zeros(size(label_fs));
%   lut = fsss_read_lut;
%   
%   % saving structure names:
%   % new label   old label   structure name   #voxels
%   fid = fopen([prefix,segname,'_reidx.txt'],'w');
%   j=1;
%   for i=1:numel(label_fs)
%    idx = aparc(:)==label_fs(i);
%    nvox_fs(i) = sum(idx);
%    strname_fs{i} = lut.strcname{label_fs(i)==lut.label};
%    nonbrain = 0;
%    for m=1:numel(Nonbrain)
%     if strfind(lower(strname_fs{i}), lower(Nonbrain{m}))
%      nonbrain = 1;
%     end
%    end
%    if ~nonbrain
%     strname = strname_fs{i};
%     if strcmpi(strname(1:3),'ctx')
%      strname = strname(5:end);
%     end
%     fprintf(fid,'%d %d %s %d\n', j, label_fs(i), strname, nvox_fs(i));
%     nii1.img(idx) = j;
%     j=j+1;
%    end
%   end
%   fclose(fid);
%   fname_aparc1=['o',segname,'_reidx.nii'];
%   save_untouch_nii(nii1, fname_aparc1);
%  end
% end
% cd(path0)
% 
% end
