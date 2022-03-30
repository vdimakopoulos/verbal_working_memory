%% Close all figures, clear variables and command window
close all
clear
clc


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
addpath(strPaths.Project)
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

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))


%% Figure Data

strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
strFilename = [strPaths.FigureData,'PLV_Granger_Figure_S_2_9\PLV_Granger_var'];
load(strFilename);

%% Figure Axes

fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[4.8948,2.8581,19,20]);
h=gcf;


ha = tight_subplot(3,3,[0.05 0.05],[0.05 0.03],[0.05 0.05])
set(ha(1:size(ha,1)),'Fontsize',10)


%% Panel Labels
if 1
    % Create textbox
    annotation('textbox',...
        [0.00556973684210528 0.967345000749095 0.0459539464138971 0.0383645825842396],...
        'String',{'a'},...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    annotation('textbox',...
        [0.31889210526316 0.967345000749095 0.0459539464138972 0.0383645825842396],...
        'String','b',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    % Create textbox
    annotation('textbox',...
        [0.633607017543868 0.967345000749095 0.0459539464138972 0.0383645825842396],...
        'String','c',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    
    % Create textbox
    annotation('textbox',...
        [0.00556973684210529 0.647199167415776 0.0459539464138971 0.0383645825842396],...
        'String','d',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    % Create textbox
    annotation('textbox',...
        [0.323069736842107 0.647199167415776 0.0459539464138971 0.0383645825842396],...
        'String','e',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    % Create textbox
    annotation('textbox',...
        [0.637784649122815 0.647199167415776 0.0403837711413048 0.0383645825842396],...
        'String','f',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    % Create textbox
    annotation('textbox',...
        [0.00556973684210529 0.323084584082456 0.0459539464138971 0.0383645825842396],...
        'String','g',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    % Create textbox
    annotation('textbox',...
        [0.321677192982458 0.323084584082456 0.0459539464138972 0.0383645825842396],...
        'String','h',...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
 
end

%% Global Data
freqAxis_fix = PLV_fig_var.global_vars.freqAxis_fix;
freqAxis = PLV_fig_var.global_vars.freqAxis;
ylimit = PLV_fig_var.global_vars.ylim;
xlimit = PLV_fig_var.global_vars.xlim;


%% Panel a
axes(ha(1));
PLV_Maint_Subj{1} = PLV_fig_var.P37.maint;
PLV_Fix_Subj{1} = PLV_fig_var.P37.fixation;
f = semilogx(freqAxis_fix,PLV_Fix_Subj{1},'k','LineWidth',3);
hold on;
m = semilogx(freqAxis_fix,PLV_Maint_Subj{1},'r','LineWidth',3);

semilogx([5 8],ones(length([5 8]),1)'*0.05,'k','LineWidth',3,'HandleVisibility','off');
semilogx([10 14],ones(length([10 14]),1)'*0.05,'r','LineWidth',3,'HandleVisibility','off');
semilogx([16 21],ones(length([16 21]),1)'*0.05,'r','LineWidth',3,'HandleVisibility','off');

legend1 = legend(ha(1),[f m],'Fixation','Maintenance');
set(legend1,...
    'Position',[0.053376227701882,0.863727123589372,0.18802228063594,0.056216929718931],...
    'FontSize',11,...
    'EdgeColor','none');
strXlabel(1) = xlabel('Frequency (Hz)');
strYlabel(1) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S2');


set(ha(1),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'XMinorTick','on',...
    'YTick',[0 0.5 1]);

set(ha(1),'FontSize',10)
set(ha(1),'YTickLabelMode','auto')
set(ha(1),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(1),'Box','off');
set(strYlabel(1),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(1),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(1),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.049198596122935,0.913997956922706,0.236272456508644,0.056216929718931],...
    'String',{'Hipp - Parietal Left'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Panel b
axes(ha(2));
PLV_Maint_Subj{2} = PLV_fig_var.P44.Maint
PLV_Fix_Subj{2} = PLV_fig_var.P44.Fixation;
semilogx(freqAxis_fix,PLV_Fix_Subj{2},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{2},'r','LineWidth',3);
semilogx(PLV_fig_var.P44.significance_bar,ones(length(PLV_fig_var.P44.significance_bar),1)*0.05,'r','LineWidth',3);
semilogx([30 45],ones(length([30 45]),1)*0.05,'r','LineWidth',3);

strXlabel(2) = xlabel('Frequency (Hz)');
strYlabel(2) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S3');


set(ha(2),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(2),'FontSize',10)
set(ha(2),'YTickLabelMode','auto')
set(ha(2),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(2),'Box','off');
set(strYlabel(2),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(2),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(2),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.361128420684341,0.913997956922706,0.236272456508646,0.056216929718931],...
    'String',{'Hipp - Parietal Scalp'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Panel c
axes(ha(3));
PLV_Maint_Subj{3} = PLV_fig_var.P38.Maint
PLV_Fix_Subj{3} = PLV_fig_var.P38.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{3},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{3},'r','LineWidth',3);
semilogx([4 9],ones(length([4 9]),1)*0.05,'r','LineWidth',3);
semilogx([16 21],ones(length([16 21]),1)*0.05,'k','LineWidth',3);

strXlabel(3) = xlabel('Frequency (Hz)');
strYlabel(3) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S4');


set(ha(3),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(3),'FontSize',10)
set(ha(3),'YTickLabelMode','auto')
set(ha(3),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(3),'Box','off');
set(strYlabel(3),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(3),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(3),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.680020964543996,0.913997956922706,0.236272456508648,0.056216929718931],...
    'String',{'Hipp - Parietal Scalp'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Panel d
axes(ha(4));
PLV_Maint_Subj{4} = PLV_fig_var.P12.Maint
PLV_Fix_Subj{4} = PLV_fig_var.P12.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{4},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{4},'r','LineWidth',3);
indMaxBand = Difference_Bar(PLV_Maint_Subj{4} -PLV_Fix_Subj{4},PLV_fig_var.P12.percentile{95},freqAxis,[0 1],0.05,'r');
indMaxBand2 = Difference_Bar(PLV_Fix_Subj{4} - PLV_Maint_Subj{4},PLV_fig_var.P12.percentile{95},freqAxis,[0 1],0.05,'k');
strXlabel(4) = xlabel('Frequency (Hz)');
strYlabel(4) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S5','Position',[20.000042209233563,0.998792287232219,0]);

set(ha(4),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(4),'FontSize',10)
set(ha(4),'YTickLabelMode','auto')
set(ha(4),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(4),'Box','off');
set(strYlabel(4),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(4),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(4),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.049198596122935,0.593852123589384,0.261338245982328,0.056216929718931],...
    'String',{'Frontal Depth - Sensory Cortex'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Panel e
axes(ha(5));
PLV_Maint_Subj{5} = PLV_fig_var.P20.Maint
PLV_Fix_Subj{5} = PLV_fig_var.P20.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{5},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{5},'r','LineWidth',3);
semilogx([6 11],ones(length([6 11]),1)*0.05,'r','LineWidth',3)
strXlabel(5) = xlabel('Frequency (Hz)');
strYlabel(5) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S6','Position',[20.000042209233563,0.998792287232219,0]);

set(ha(5),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(5),'FontSize',10)
set(ha(5),'YTickLabelMode','auto')
set(ha(5),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(5),'Box','off');
set(strYlabel(5),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(5),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(5),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.3639,0.5939,0.2613,0.0562],...
    'String',{'Temporal - Parietal Left'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Panel f
axes(ha(6));
PLV_Maint_Subj{6} = PLV_fig_var.P26.Maint
PLV_Fix_Subj{6} = PLV_fig_var.P26.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{6},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{6},'r','LineWidth',3);
semilogx(PLV_fig_var.P26.significance_bars(1,:),ones(length(PLV_fig_var.P26.significance_bars(1,:)),1)*0.05,'r','LineWidth',3);
semilogx(PLV_fig_var.P26.significance_bars(2,:),ones(length(PLV_fig_var.P26.significance_bars(2,:)),1)*0.05,'r','LineWidth',3);
semilogx(PLV_fig_var.P26.significance_bars_fix(1,:),ones(length(PLV_fig_var.P26.significance_bars_fix(1,:)),1)*0.05,'k','LineWidth',3);
semilogx(PLV_fig_var.P26.significance_bars_fix(2,:),ones(length(PLV_fig_var.P26.significance_bars_fix(2,:)),1)*0.05,'k','LineWidth',3);

strXlabel(6) = xlabel('Frequency (Hz)');
strYlabel(6) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S7','Position',[20.000042209233563,0.998792287232219,0]);

set(ha(6),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(6),'FontSize',10)
set(ha(6),'YTickLabelMode','auto')
set(ha(6),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(6),'Box','off');
set(strYlabel(6),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(6),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(6),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

%Create textbox
annotation('textbox',...
    [0.678614912280705,0.593900000000001,0.261300000000001,0.0562],...
    'String',{'Parietal Depth - Sensory Cortex'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Panel g
axes(ha(7));
PLV_Maint_Subj{7} = PLV_fig_var.P40.Maint;
PLV_Fix_Subj{7} = PLV_fig_var.P40.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{7},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{7},'r','LineWidth',3);
semilogx(PLV_fig_var.P40.significance_bar,ones(length(PLV_fig_var.P40.significance_bar),1)*0.05,'r','LineWidth',3);
strXlabel(7) = xlabel('Frequency (Hz)');
strYlabel(7) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S8','Position',[20.000042209233563,0.998792287232219,0]);


set(ha(7),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(7),'FontSize',10)
set(ha(7),'YTickLabelMode','auto')
set(ha(7),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(7),'Box','off');
set(strYlabel(7),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(7),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(7),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.049198596122935,0.263122956922722,0.261338245982328,0.056216929718931],...
    'String',{'Hipp - Parietal Scalp'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Panel h
axes(ha(8));
PLV_Maint_Subj{8} = PLV_fig_var.P10.Maint;
PLV_Fix_Subj{8} = PLV_fig_var.P10.Fix;
semilogx(freqAxis_fix,PLV_Fix_Subj{8},'k','LineWidth',3);
hold on;
semilogx(freqAxis,PLV_Maint_Subj{8},'r','LineWidth',3);
semilogx(PLV_fig_var.P10.significance_bar(1,:),ones(length(PLV_fig_var.P10.significance_bar(1,:)),1)*0.05,'r','LineWidth',3);
strXlabel(8) = xlabel('Frequency (Hz)');
strYlabel(8) = ylabel('PLV');
ylim([ylimit(1) ylimit(2)]);
xlim([xlimit(1) xlimit(2)]);
title('S9','Position',[20.000042209233563,0.998792287232219,0]);

set(ha(8),'XTick',...
    [10 100],'XTickLabel',{'10','100'},'YTick',[0 0.5 1]);

set(ha(8),'FontSize',10)
set(ha(8),'YTickLabelMode','auto')
set(ha(8),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
set(ha(8),'Box','off');
set(strYlabel(8),'Position',[3.297861279709002,0.252427661303178,-1],'FontSize',11);
set(strXlabel(8),'Position',[30.997466460963274,-0.011909383349121,-1],'FontSize',11);
set(ha(8),'XMinorTick','on');
xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};
for k = 1:length(xTick)
    text(xTick(k),-0.02,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
end

% Create textbox
annotation('textbox',...
    [0.36252096454399,0.263122956922722,0.261338245982331,0.056216929718931],...
    'String',{'Hipp - Frontal'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Save Fig in Different Formats

set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');
set(ha(1:size(ha)),'box','off');
set(ha(end),'XColor','none', 'YColor','none')
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV Across Subjects\PLV_S2_9','-dpdf','-r2000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV Across Subjects\PLV_S2_9','-dpng','-r2000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV Across Subjects\PLV_S2_9','-depsc','-r2000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV Across Subjects\PLV_S2_9.fig');