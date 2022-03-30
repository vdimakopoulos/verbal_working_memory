%% Paths
strPaths.Main = 'F:\Vasileios\';
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
strPaths.GeneralFunctions = 'F:\Vasileios\Task Analysis\Code\';
strPaths.Data = [strPaths.Main, 'Task Analysis\Data\'];
% Results
strPaths.Results = [strPaths.Main,'Task Analysis\Analysis Results\'];
strPaths.FigureData = [strPaths.Data ,'Analysis Data\Figure Data\']
% FieldTrip toolbox
% strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20191126\';
strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\';

% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = 'F:\Vasileios\Toolboxes\eeglab14_1_1b\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(strPaths.Main)
addpath(genpath(strPaths.Project))
addpath(genpath(strPaths.GeneralFunctions))
addpath(strPaths.Data)
% addpath(strPaths.Subfunctions)
addpath(strPaths.Results)
addpath(strPaths.Toolboxes.FieldTrip)
% Remove EEGLAB from path

% rmpath(genpath(strPaths.Toolboxes.EEGLAB))


% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults
OmitScalp =1;
%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))


%% Import anatomical locations table and MR image

strPaths.AnatomicalLocationTable = [strPaths.Data,'Sternberg Task\FirstContactHipp_electrodes_AnatomicalLocations_13subjects.xlsx'];
% strPaths.AnatomicalLocationTable = [strPaths.Data,'Sternberg Task\FirstContactHipp_electrodes_AnatomicalLocations_KJsubject.xlsx'];

strPaths.MRI_Image = [strPaths.FigureData,'Task_figure\Images\Screenshot empty brain-2.jpg'];

strPaths.VarData = [strPaths.FigureData,'Anatomical_Locations.mat']


Anatomical_Locations = readtable(strPaths.AnatomicalLocationTable);
ImgSagittal = imread(strPaths.MRI_Image);
load(strPaths.VarData);

%% Plot Anatomical locations on sagittal slice
fig = figure;


im = imagesc(pXAxSagittal,pYAxSagittal,ImgSagittal)

set(gca,'YDir','Normal')


hold on;

nMarkerSize = 0.05;
s = [];
Temp = Anatomical_Locations;
markerColor = 'r';
for iLoc =1:size(Temp)%[1:40 42:size(Temp,1)]
%    if ~isnan(Temp.Used(iLoc) )
%        markerColor = 'b'
%    else
%        markerColor = 'r'
%    end
    s = scatter((Temp.y(iLoc)),...
        (Temp.z(iLoc)),'Marker','o','MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor,'LineWidth',nMarkerSize);

end

a = 2;
set(fig,'Position',[400,50,size(ImgSagittal,2)*a,size(ImgSagittal,1)*a])
set(gca,'Units','pixels')
a = 1.5; % 1.5
set(gca,'Position',[100,130,size(ImgSagittal,2)*a,size(ImgSagittal,1)*a])
set(gca,'Units','normalized')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'XTick',[])
set(gca,'YTick',[])
im_nan = ones(size(ImgSagittal,1),size(ImgSagittal,2));
im_nan(ImgSagittal(:,:,1)==1) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==2) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==3) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==4) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==5) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==6) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==7) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==8) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==9) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)
im_nan(ImgSagittal(:,:,1)==10) = 0;% all 'non-color [255 here, NaN in your case] should be transparent (0)


im.AlphaData = im_nan;
%% Save Figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_Task_electrode_locations','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_Task_electrode_locations','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_Task_electrode_locations','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_Task_electrode_locations.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_Task_electrode_locations.tif');


% 
% 
% 
% print(gcf,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\ReactionTime_SetSize','-dpdf','-r1000');
% print(gcf,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\ReactionTime_SetSize','-dpng','-r1000');
% print(gcf,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\ReactionTime_SetSize','-depsc','-r1000');
% saveas(gcf,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\ReactionTime_SetSize.fig');
% saveas(gcf,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\ReactionTime_SetSize.tif');
