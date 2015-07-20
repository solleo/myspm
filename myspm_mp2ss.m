function EXP = myspm_mp2ss (EXP)
% EXP = myspm_mp2ss (EXP)
%
% [Input]
% EXP
% .fname_t1w
% .fname_inv2
% .fname_out
%
% (cc) 2015. sgKIM. solleo@gmail.com

[p1,~,~] = fileparts(EXP.fname_t1w);
cd (p1);
nii1 = load_untouch_nii (EXP.fname_t1w);
nii2 = load_untouch_nii (EXP.fname_inv2);
nii = nii1;
nii.img = double(nii1.img .* nii2.img);
save_untouch_nii(nii, 'prod.nii');

vbm8_brainExt({'prod.nii'});
nii1 = load_untouch_nii ('p1prod.nii');
nii2 = load_untouch_nii ('p2prod.nii');
nii3 = load_untouch_nii ('p3prod.nii');
nii4 = load_untouch_nii ('mprod.nii');
nii = nii1;
nii.img =  double(~~(nii1.img + nii2.img + nii3.img)) .* double(nii4.img);
save_untouch_nii(nii, 'mp2_brain.nii');

if isfield(EXP,'fname_out')
  movefile('mp2_brain.nii', EXP.fname_out);
end

end
