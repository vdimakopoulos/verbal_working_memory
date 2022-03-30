
%%
close all
clear
clc

%% Paths
% Folder that contains the folder 'Sternberg Task'
Drive_Letter = 'F:\';
strPaths.Main = 'F:\Vasileios\Task Analysis\';
strPaths.Code = [strPaths.Main,'Code\'];
strPaths.NiftiDir = [strPaths.Main,'Data\Imaging Data\'];
%Toolboxes
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = [Drive_Letter,'Vasileios\Toolboxes\fieldtrip-20200315\'];

%Freesurfer
% strPaths.Toolboxes.FreeSurfer              = [strPaths.Main,'Programs and utilities\freesurfer\'];
strPaths.Toolboxes.SPM = 'F:\Vasileios\Toolboxes\spm12';
addpath(strPaths.Toolboxes.SPM)

cd(strPaths.Main)
% Add all subfolders to path
addpath(genpath(strPaths.Code))
addpath(strPaths.Toolboxes.FieldTrip)
% addpath(genpath(strPaths.Toolboxes.FreeSurfer))
addpath(genpath(strPaths.NiftiDir));

ft_defaults

%% Preprocessing of the anatomical MRI
strMRIFilePath = [strPaths.NiftiDir,'45 SS Imaging\Nifti\MRI Post\t1_mprage_tra_0_9.nii'];
mri = ft_read_mri(strMRIFilePath); % we used the dcm series
ft_determine_coordsys(mri);
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);

%% Saving aligned MR nifti file
strOutputPath = [strPaths.Main 'Data\Imaging Data\'];
strPatientPath = [strOutputPath '45 SS Imaging\Nifti aligned\MRI Post\'];
mkdir(strPatientPath)


cfg = [];
cfg.filename = [strPatientPath 'MR_acpc_200824'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);


%% Preprocessing of the anatomical CT
strCTFilePath = [strPaths.NiftiDir,'45 SS Imaging\Nifti\CT Post\hr_Head_Spi_nativ_0_80_Hr60_A1.nii'];
ct = ft_read_mri(strCTFilePath)
cfg = [];
cfg.method = 'interactive';
cfg.coordsys =  'ctf';
ct_ctf = ft_volumerealign(cfg, ct);
ct_acpc = ft_convert_coordsys(ct_ctf,'acpc');

%% Saving aligned CT nifti file
strOutputPath = [strPaths.Main 'Data\Imaging Data\'];
strPatientPath = [strOutputPath '45 SS Imaging\Nifti aligned\CT Post\'];
mkdir(strPatientPath)


cfg = [];
cfg.filename = [strPatientPath 'CT_acpc_200824'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc);


%% Fusion of the CT with the MRI
strMR_acpcFileName = [strPaths.Main 'Data\Imaging Data\45 SS Imaging\Nifti aligned\MRI Post\MR_acpc_200824.nii'];
fsmri_acpc = ft_read_mri(strMR_acpcFileName);
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);


strFusionImagePath = [strPaths.Main, 'Data\Imaging Data\45 SS Imaging\Nifti Fused\']
mkdir(strFusionImagePath)
cfg = [];
cfg.filename = [strFusionImagePath,'SS_CT_acpc_MR_post'] 

cfg.filetype = 'nifti'; 
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);
print([cfg.filename],'-dpng');
