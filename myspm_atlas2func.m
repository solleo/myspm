%function EXP = myspm_atlas2func (EXP)

subjid='2001';
EXP.name_t1w='Brain.nii';
EXP.name_epi='uarest410.nii';
EXP.dir_base='/scr/vatikan3/APConn/rest12.410/';
EXP.prefix_ants='/scr/vatikan3/APConn/ants_t1w_n17/t1w';
EXP.fname_atlas='/scr/vatikan1/skim/matlab/conn/rois/atlas.nii';

path1=[fullfile(EXP.dir_base,subjid),'/'];
EXP.fname_t1w=[path1,EXP.name_t1w];
[~,name1,~]=fileparts(EXP.fname_t1w);
EXP.fname_epi=[path1,'mean',EXP.name_epi];
[~,name3,~]=fileparts(EXP.fname_atlas);

%%

% 1. 
% atlas@mni => @t1w
exp1=[];
exp1.fname_ref = EXP.fname_t1w;
exp1.fname_in  = EXP.fname_atlas;
exp1.fname_out = [path1,name3,'_in_',name1,'.nii'];
exp1.xfm(1).prefix = [EXP.prefix_ants,subjid];
exp1.xfm(1).isInverse = 1; % MNI to t1w
exp1.intp='NearestNeighbor';
myants_applyxfm(exp1);

% atlas@t1w => @func
exp1=[];
exp1.name_fixed  = EXP.fname_epi;
exp1.name_moving = EXP.fname_t1w;
exp1.name_others = [path1,name3,'_in_',name1,'.nii'];
exp1.interp=0;
myspm_coreg(exp1);

% mask it with GM>0.5
EXP.fname_atlas = [path1,'o',name3,'_in_',name1,'.nii'];
gm = load_uns_nii([path1,'oc1t1w.nii']); % in integer format!
gm.img = zeroone(gm.img);
atlas = load_uns_nii(EXP.fname_atlas);
atlas.img = uint16((gm.img>0.5) .* atlas.img);
atlas.hdr.dime.datatype=512;
save_untouch_nii(atlas, EXP.fname_atlas);

%%
% %% But not a good idea...? There are some errors, but also in MNI space too.
% % it should be better to warp it, and mask with GM>0.50
% exp1.name_fixed  = EXP.fname_t1w;
% exp1.name_moving = EXP.fname_epi;
% exp1.prefix='g';
% myspm_coreg(exp1)
% exp1=[];
% exp1.name_def='y_t1w.nii';
% exp1.name_moving='gmeanuarest410.nii';
% exp1.vox_mm=[2 2 2];
% tic
% myspm_deform1(exp1)
% 
