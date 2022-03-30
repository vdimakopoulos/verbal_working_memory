%%
% Script for transforming DICOM images to nifti

%%
close all
clear
clc

%% Paths
% Folder that contains the folder 'Sternberg Task'
strPaths.Main = 'F:\Vasileios\';
strPaths.Project = [strPaths.Main, 'Task Analysis'];
strPaths.Code = [strPaths.Project,'Code\'];
strPaths.NiftiCode = [strPaths.Main, 'Toolboxes\dicm2nii-master\'];
%Toolboxes
%  strPaths.Toolboxes                        = [strPaths.Main,'Programs and utilities\'];
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip              = [strPaths.Main,'Toolboxes\fieldtrip-20200315\'];
%Freesurfer
strPaths.Toolboxes.FreeSurfer              = [strPaths.Main,'Programs and utilities\freesurfer\'];

% Results
strPaths.Results = [strPaths.Main,'Data\Imaging Data\'];

cd(strPaths.Main)
% Add all subfolders to path
addpath(genpath(strPaths.Code))
addpath(genpath(strPaths.NiftiCode))
addpath(strPaths.Toolboxes.FieldTrip)
addpath(genpath(strPaths.Toolboxes.FreeSurfer))
addpath(genpath(strPaths.Results));

ft_defaults


%% Convert DICOM images to Nifti

%Input paths for DICOM images
strPaths.MainInputImagesPath = 'F:\Vasileios\Task Analysis\Data\Imaging Data\';
strPaths.CT_PostOP_input = [strPaths.MainInputImagesPath, '45 SS Imaging\CT Post\'];
strPaths.MRI_PostOp_input = [strPaths.MainInputImagesPath, '45 SS Imaging\MR Post\'];

strPaths.CT_PostOp.output = [strPaths.MainInputImagesPath, '45 SS Imaging\Nifti\CT Post\']
strPaths.MRI_PostOp.output = [strPaths.MainInputImagesPath, '45 SS Imaging\Nifti\MRI Post\']
mkdir(strPaths.CT_PostOp.output)
mkdir(strPaths.MRI_PostOp.output)

varargoutCT = dicm2nii(strPaths.CT_PostOP_input, strPaths.CT_PostOp.output, '.nii')
varargoutMRI = dicm2nii(strPaths.MRI_PostOp_input, strPaths.MRI_PostOp.output, '.nii')