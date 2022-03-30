%%
close all
clear
clc

%% Paths
% Folder that contains the folder 'Sternberg Task'
strPaths.Main = 'F:\Vasileios\Task Analysis\';
strPaths.Code = [strPaths.Main,'Code\'];
strPaths.FusedImgPath = [strPaths.Main,'\Data\Imaging Data\45 SS Imaging\Nifti Fused\'];
%Toolboxes
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip              = ['F:\Vasileios\Toolboxes\','\fieldtrip-20210529\'];
%Freesurfer
% strPaths.Toolboxes.FreeSurfer              = [strPaths.Main,'Programs and utilities\freesurfer\'];


cd(strPaths.Main)
% Add all subfolders to path
addpath(genpath(strPaths.Code))
addpath(genpath(strPaths.Toolboxes.FieldTrip))
% addpath(genpath(strPaths.Toolboxes.FreeSurfer))
addpath(genpath(strPaths.FusedImgPath));

ft_defaults


%% Load merged Images
strCT_merged_file = [strPaths.FusedImgPath '\SS_CT_acpc_MR_post.nii']
strMRI_acpc_file = [strPaths.Main, 'Data\Imaging Data\45 SS Imaging\Nifti Aligned\MRI Post\MR_acpc_200824.nii'];%'Data\Imaging Data\49 TK Imaging\Nifti aligned\CT Post\CT_acpc_210614.nii']

ct_acpc_f = ft_read_mri(strCT_merged_file);
fsmri_acpc = ft_read_mri(strMRI_acpc_file);

%% Load the header file

strPath_TaskMicroData = 'F:\Vasileios\Task Analysis\Data\Imaging Data\45 SS Imaging\45_SS_header'
load(strPath_TaskMicroData);
% hdr = ft_read_header(strPath_TaskMicroData);
hdr.label = strrep(hdr.label,'u','');

%% Place the electrode contacts on the MR/CT images
cfg = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, fsmri_acpc, fsmri_acpc);
strPath_elec_placement = [strPaths.Main, '\Data\Imaging Data\45 SS Imaging\Elec_placement\'];
mkdir(strPath_elec_placement);
save([strPath_elec_placement, 'SS_elec_acpc_f_post.mat'], 'elec_acpc_f');
