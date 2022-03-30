clc;
close all;
clear;

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
%% Figure
%% Make Figure
% fig = figure;
% set(fig,'Units', 'centimeters')
% % set(fig,'Position',[8.440208333333334,3.01625,20.21416666666667,22.11916666666667]);%[4.8948,2.8581,19,20]);
% h=gcf;
% 
% ha = tight_subplot(1,3,[.1 .08],[.12 .04],[.08, .05])
% set(fig,'Position',[11.482916666666668,14.73729166666667,21.29895833333333,11.112500000000002]);
% 

%% Figure Var Data
strPaths.FigureVarData = [strPaths.FigureData,'Mean TFR Granger\TFR_Granger_mean_Fig2.mat'];
load(strPaths.FigureVarData);

clim = TFR_Granger_mean_Fig.clim*100;
freqAxis = TFR_Granger_mean_Fig.freqAxis;
timeAxis = TFR_Granger_mean_Fig.timeAxis;
nSubjects = TFR_Granger_mean_Fig.nSubjects;
MeanTFR_Granger = TFR_Granger_mean_Fig.meanTFR;
cmap = TFR_Granger_mean_Fig.cmap;
%% Figure
fig = figure;
% axes(ha(1));
set(fig,'Position',[400,50,759,640.5])
set(gcf,'Color','white')
set(gca,'Units','pixels')

contourf(timeAxis,freqAxis,(MeanTFR_Granger./nSubjects)*100,100,'LineColor','none');
colormap(cmap);
ylim([4 100]);
ylim([4 30])
cb = colorbar;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
set(cb,'Ticks', [-2 -1 0 1 2],'TickLabels', [-2 -1 0 1 2]);
hold on;
line([-5 -5],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)
line([-3 -3],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)
line([-0 0],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)

%% Create textboxes

% Create textbox
annotation(fig,'textbox',...
    [0.757575757575758 0.932645590918779 0.271409741974631 0.0671350492061256],...
    'String',{'Hipp -> EEG'},...
    'FontSize',20,...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.757575757575758 0.0676963325268956 0.271409741974631 0.0671350492061255],...
    'String',{'EEG -> Hipp'},...
    'FontSize',20,...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.160737812911728,0.130147542519094,0.491436085878468,0.067135049206126],...
    'String',{'fix    encoding   maintenance'},...
    'FontSize',20,...
    'EdgeColor','none');

% annotation(fig,'textbox',...
%     [0.00131752305665188 0.932645591384081 0.0757575738924766 0.0827478512388762],...
%     'String',{'c'},...
%     'FontSize',26,...
%     'EdgeColor','none');

% Create textbox
h = text(3.151471064640926,7.581781624042329,0.237154126167283,'\DeltaGranger (%)')
set(h,'Rotation',90,'FontSize',20);

%% Axes properties
set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
set(cb,'Ticks',[-0.02, -0.01, 0, 0.01, 0.02]*100,'TickLabels',[-0.02, -0.01, 0, 0.01, 0.02]*100);
set(gca,'FontSize',24,'box','off');
set(gca,'box','off','TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.04 0.025])


%% Save Figure
set(gcf, 'InvertHardcopy', 'off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Mean TFR Granger\MeanTFR_Granger','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Mean TFR Granger\MeanTFR_Granger','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Mean TFR Granger\MeanTFR_Granger','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Mean TFR Granger\MeanTFR_Granger.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Mean TFR Granger\MeanTFR_Granger.tif');
