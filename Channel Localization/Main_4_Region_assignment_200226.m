%%
close all
clear
clc

%% Paths
% Folder that contains the folder 'Sternberg Task'
strPaths.Main = 'F:/Vasileios/Task Analysis/';
strPaths.Code = [strPaths.Main,'Code\'];
strPaths.FusedImgPath = [strPaths.Main,'4 Imaging Data ACPC CT Aligned\'];

%Toolboxes
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip              = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\';
%Freesurfer
strPaths.Toolboxes.FreeSurfer              = [strPaths.Main,'Programs and utilities\freesurfer\'];


cd(strPaths.Main)
% Add all subfolders to path
addpath(genpath(strPaths.Code))
addpath(strPaths.Toolboxes.FieldTrip)
addpath(genpath(strPaths.Toolboxes.FreeSurfer))
addpath(genpath(strPaths.FusedImgPath));
ft_defaults
%% Normalization and coordinate transformation

strMRI_acpc_file = [strPaths.Main, 'Data\Imaging Data\45 SS Imaging\Nifti Aligned\MRI Post\MR_acpc_200824.nii']%'Data\Imaging Data\49 TK Imaging\Nifti aligned\MRI Post\MR_acpc_210614.nii']
strPath_elec_placement = [strPaths.Main, '\Data\Imaging Data\45 SS Imaging\Elec_placement\'];
elec_placement = load([strPath_elec_placement, 'SS_elec_acpc_f_post.mat'])
elec_acpc_f = elec_placement.elec_acpc_f;
fsmri_acpc = ft_read_mri(strMRI_acpc_file)

cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.coordsys =  'ctf';
fsmri_mni = ft_volumenormalise(cfg,fsmri_acpc);

elec_acpc_fr = elec_acpc_f;
elec_mni_frv = elec_acpc_fr;

elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, ...
    elec_acpc_fr.elecpos, 'individual2sn');

elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, ...
    elec_acpc_fr.chanpos, 'individual2sn');

elec_mni_frv.coordsys = 'mni';
strPath_elec_placement = [strPaths.Main,'\Data\Imaging Data\45 SS Imaging\Volume Based Registration\']
mkdir(strPath_elec_placement)
save([strPath_elec_placement,'45 SS elec_mni_frv_post.mat'], 'elec_mni_frv');


%% Read atlas and Region assignment

atlas = ft_read_atlas([strPaths.Toolboxes.FieldTrip,'template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii']);
for ch  = 1:length(elec_mni_frv.chanpos)
    if ~isnan(elec_mni_frv.chanpos(ch,1))
        cfg                     = [];
        cfg.roi                 = elec_mni_frv.chanpos( ...
            match_str(elec_mni_frv.label, char(elec_mni_frv.label{ch})),:);
        cfg.atlas               = atlas;
        cfg.inputcoord          = 'mni';
        cfg.output              = 'label';
        labels                  = ft_volumelookup(cfg, atlas);
        [~, indx]               = max(labels.count);
        elec_mni_frv.loc{ch,2}  = labels.name(indx)
        elec_mni_frv.loc{ch,1}  = elec_mni_frv.label{ch}
    end
end
strPath_AssignedRegions = [strPaths.Main, 'Data\Imaging Data\45 SS Imaging\Anatomical Labeling\'];
mkdir(strPath_AssignedRegions)
save([strPath_AssignedRegions, '45_SS_elec_mni_frv_post_labeled.mat'],'elec_mni_frv');