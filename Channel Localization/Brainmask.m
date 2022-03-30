%% Paths
strPaths.Main = 'F:\Vasileios\';
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
strPaths.GeneralFunctions = 'F:\Vasileios\Task Analysis\Code\';
strPaths.Data = [strPaths.Main, 'Task Analysis\Data\'];
strPaths.ExtractedData = [strPaths.Main, 'Task Analysis\Extracted_Data Sternberg\'];
strPaths.Statistics = [strPaths.Main, 'Task Analysis\Code\Statistics\'];
strPaths.ChanLoc = [strPaths.Main, 'Task Analysis\Code\Channel Localization\'];
strPaths.EEGLAB_Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\EEGLAB Subfunctions\'];
strPaths.Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\']
% Results
strPaths.Results = [strPaths.Main,'Task Analysis\Analysis Results\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20191126\';
% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = 'F:\Vasileios\Toolboxes\eeglab14_1_1b\';
strPaths.Toolboxes.iElectrodes          = 'F:\Vasileios\Toolboxes\iElectrodes';
% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(strPaths.Main)
addpath(strPaths.Project)
addpath(strPaths.GeneralFunctions)
addpath(strPaths.Data)
addpath(strPaths.ExtractedData)
addpath(strPaths.Statistics)
addpath(strPaths.ChanLoc)
addpath(strPaths.Subfunctions)
addpath(strPaths.EEGLAB_Subfunctions)
addpath(strPaths.Results)
addpath(strPaths.Toolboxes.FieldTrip)
addpath(strPaths.Toolboxes.iElectrodes)

% Remove EEGLAB from path

rmpath(genpath(strPaths.Toolboxes.EEGLAB))

Figures_Path = 'F:\Vasileios\Task Analysis\Analysis Results\Statistics\PLV_significance\'

strPlotColors = {'b','g','r','c','k'};
ft_defaults


%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Convert DICOM images to Nifti

%Input paths for DICOM images
strPaths.MainInputImagesPath = 'F:\Vasileios\Task Analysis\Data\42 DS Imaging\';
strPaths.MRI_PostOp_input = [strPaths.MainInputImagesPath, 'DICOM\MR\PostOp MR'];
strPaths.MRI_PostOp_input = [strPaths.MainInputImagesPath, 'DICOM\MR\'];

strPaths.CT_PostOp_input = [strPaths.MainInputImagesPath, 'DICOM\CT\'];


strPaths.MRI_PostOp_output = [strPaths.MainInputImagesPath, 'NIFTI\MR\'];
strPaths.CT_PostOp_output = [strPaths.MainInputImagesPath, 'NIFTI\CT\'];

mkdir(strPaths.MRI_PostOp_output);
mkdir(strPaths.CT_PostOp_output);


varargoutMRI = dicm2nii(strPaths.MRI_PostOp_input, strPaths.MRI_PostOp_output, '.nii')
varargoutCT = dicm2nii(strPaths.CT_PostOp_input,strPaths.CT_PostOp_output, '.nii')

%% Read NIFTI - Post MR and Post CT
MRI_nifti_filename = [strPaths.MRI_PostOp_output,'t1_mprage_tra_p2_iso_0_9_MPR_Sag.nii'];
mri_nifti = ft_read_mri(MRI_nifti_filename) 


CT_nifti_filename = [strPaths.CT_PostOp_output,'CT_Sch_del.nii'];
CT_nifti = ft_read_mri(CT_nifti_filename) 

%% MR and CT alignment

ft_determine_coordsys(mri_nifti);
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
mri_acpc = ft_volumerealign(cfg, mri_nifti);

cfg = [];
cfg.parameter = 'anatomy';
cfg.filetype = 'nifti';
cfg.filename = [strPaths.MainInputImagesPath,'NIFTI aligned\MRI Post\MR_Nifti_acpc_aligned2']
ft_volumewrite(cfg,mri_acpc)


% CT
ft_determine_coordsys(CT_nifti);
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'ctf';
cfg.viewresult = 'yes';
ct = ft_volumerealign(cfg, CT_nifti);
ct_acpc = ft_convert_coordsys(ct,'acpc');

cfg = [];
cfg.parameter = 'anatomy';
cfg.filetype = 'nifti';
cfg.filename = [strPaths.MainInputImagesPath,'NIFTI aligned\CT Post\CT_Nifti_acpc_aligned']
ft_volumewrite(cfg,ct_acpc)

%% Fuse MR with CT

cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, mri_acpc);


cfg = [];
cfg.parameter = 'anatomy';
cfg.filetype = 'nifti';
cfg.filename = [strPaths.MainInputImagesPath,'NIFTI aligned\Fused CT MR\MR_CT_Nifti_acpc_aligned2']
ft_volumewrite(cfg,ct_acpc_f)

%% Electrode placement

strPath_TaskMicroData = 'F:/Vasileios/Task Analysis/Data/Macro Data/Macro_Data_Sessions_Patient_42_Session_01_Part_01.mat'
hdr = dataMacro.hdr;
for i = 9:72
hdr.label{i} = strrep(hdr.label{i},'_GL_','');
end

cfg = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg,ct_acpc_f, mri_acpc);

cfg = [];
cfg.channel = hdr.label;
elec_acpc_f2 = ft_electrodeplacement(cfg,ct_acpc_f, mri_acpc);

strPath_elec_placement =    'F:/Vasileios/Task Analysis/Data/Electrodes positions/';
mkdir(strPath_elec_placement);
save([strPath_elec_placement, 'DS_grid_elec_acpc_f_post_2nd_attempt.mat'], 'elec_acpc_f2');

%% Brai segmentation
mri_file =  'F:\MRI Post\MR_acpc_200225.nii';
mri = ft_read_mri(mri_file)
mri.coordsys = 'acpc'
cfg.output = {'brain'};
segmented  = ft_volumesegment(cfg, mri_nifti)


%Prepare headmodel
cfg          = [];
cfg.method   = 'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, segmented);

% Plot headmodel and grid electrodes
figure;
ft_plot_headmodel(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5;
material dull;
camlight
hold on;
labels = elec_acpc_f.label;
plot3(elec_acpc_f2.elecpos(:,1),elec_acpc_f2.elecpos(:,2),elec_acpc_f2.elecpos(:,3),'m*');
text(elec_acpc_f2.elecpos(:,1),elec_acpc_f2.elecpos(:,2),elec_acpc_f2.elecpos(:,3),labels);

cfg = [];
cfg.parameter = 'brain';
cfg.filetype = 'nifti';
cfg.filename = ['F:\brainmask']
ft_volumewrite(cfg,segmented)

% 
% figure;
% for i =1:3
% subplot(2,2,i);
% ft_plot_ortho(mri_acpc.anatomy,'transform',mri_acpc.transform)
% hold on;
% labels = elec_acpc_f.label;
% plot3(elec_acpc_f.elecpos(:,1),elec_acpc_f.elecpos(:,2),elec_acpc_f.elecpos(:,3),'m*');
% text(elec_acpc_f.elecpos(:,1),elec_acpc_f.elecpos(:,2),elec_acpc_f.elecpos(:,3),labels,'color','r')
% end


cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
mri_test = ft_volumerealign(cfg, mri_acpc);
  