
%Pre
strMRIPrePath = 'Z:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Adrian Steiner\2 Imaging Data Nifti\28 WC\MR_pre\';
MRI_file = [strMRIPrePath,'t1_mprage_sag_nativ.nii'];
mri = ft_read_mri(MRI_file);

cfg = [];
mri_anon = ft_defacevolume(cfg, mri);


cfg =[]
ft_sourceplot(cfg,mri_anon)

cfg = [];
cfg.filename ='Y:\Epilepsy\_OP\_Electrode Localization\Defaced MRI\MRI_Pre_Defaced_Patient_28'
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_anon);


strMRIPreDefacedFileName = 'Y:\Epilepsy\_OP\_Electrode Localization\Defaced MRI\MRI_Pre_Defaced_Patient_28.nii'
ft_write_mri(strMRIPreDefacedFileName, mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');


mri_defaced = ft_read_mri(strMRIPreDefacedFileName)
cfg =[]
ft_sourceplot(cfg,mri_defaced)


%Post
strMRIPostPath ='Z:\NCH_FL_Forschungsprojekte\Epilepsy\_Master Students\Adrian Steiner\2 Imaging Data Nifti\28 WC\MR_post\';
MRI_Post = [strMRIPostPath, 'sT1W_3D_TFE_32chSHC.nii'];
mri = ft_read_mri(MRI_Post);


cfg = [];
mri_anon2 = ft_defacevolume(cfg, mri);

cfg =[]
ft_sourceplot(cfg,mri_anon2)


cfg = [];
cfg.filename ='Y:\Epilepsy\_OP\_Electrode Localization\Defaced MRI\MRI_Post_new_Defaced_Patient_28'
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_anon2);

strMRIPostDefacedFileName = 'Y:\Epilepsy\_OP\_Electrode Localization\Defaced MRI\MRI_Post_Defaced_Patient_28.nii'
ft_write_mri(strMRIPreDefacedFileName, mri_anon.anatomy, 'transform', mri_anon.transform, 'dataformat', 'nifti');
% 