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
addpath(strPaths.ExtractedData)
addpath(strPaths.Statistics)
addpath(strPaths.ChanLoc)
addpath(strPaths.Subfunctions)
addpath(strPaths.EEGLAB_Subfunctions)
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

%% Figure Data

strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
strFilename = [strPaths.FigureData,'42_DS_elec_mni_frv_post_labeled'];
strFilename_Img = [strPaths.FigureData,'Anatomical_Locations'];
load(strFilename);
load(strFilename_Img);


%%
% flip MR image horizontally;
%ImgBrain_flipped = flipdim(ImgBrain,2);

fig = figure;
imagesc(elec_mni_frv.chanpos(:,1),elec_mni_frv.chanpos(:,2),I)%ImgBrain)
set(gca,'YDir','Normal')

hold on
s = [];
Temp = elec_mni_frv;
for i = 11:12 %% Hippocampal contacts
   
s = scatter((Temp.elecpos(i,1)),...
    (Temp.elecpos(i,2)),'r','Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r');
s.SizeData = 20;
end
%     ll = legend('Hipp');
%     ll.Location = 'northwest'
% 
% 
% xlabel('MNI y-coordinate (mm)')
% ylabel('MNI z-coordinate (mm)')

% set(gca,'XTick',[])
% set(gca,'YTick',[])
set(gca,'FontSize',14)
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_electrode_locations','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_electrode_locations','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_electrode_locations','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Task\vWM_electrode_locations.fig');

