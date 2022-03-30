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
strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
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
strPaths.FigureVars = [strPaths.FigureData,'Main Granger Figure\MainGranger_Fig_Data.mat'];
load(strPaths.FigureVars)

strPaths.FigureVars2 = [strPaths.FigureData,'Main Granger Figure\TFR_PSD_Hipp.mat'];
load(strPaths.FigureVars2)

ElectrodeImagesPath = [strPaths.FigureData,'Main Granger Figure\'];
% strPaths.ImageP38 = [ElectrodeImagesPath,'P38_electrodes.jpg'];
% strPaths.ImageP37 = [ElectrodeImagesPath,'P37_electrodes.jpg'];

strPaths.ImageP38 = [ElectrodeImagesPath,'USZ5_TLLS.jpg'];
strPaths.ImageP37 = [ElectrodeImagesPath,'USZ4_TOL2.jpg'];

%% Read Images
% I1 = imread(strPaths.ImageP42);
I2 = imread(strPaths.ImageP38);
I3 = imread(strPaths.ImageP37);



%% Make Figure
fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[8.440208333333334,12.170833333333334,20.21416666666667,12.964583333333337]);%[4.8948,2.8581,19,20]);
h=gcf;

ha = tight_subplot(2,5,[.1 .08],[.12 .04],[.08, .05])
% ha = tight_subplot(2,5,[.01 .08],[.1 .04],[.03 .03])%[0.05 0.05],[.1 .05],[.1 .1])%[0.05 0.03],[0.05 0.05])
set(ha(1:size(ha,1)),'Fontsize',20)


%% Plot the electrode locations

axes(ha(1));
imshow(I2);
set(gca,...ok
    'Position',[0.02,0.52,0.18,0.37]);

Im2Handle = imhandles(gca);
set(Im2Handle,'XData',[1,1250],'YData',[1,1250]);

axes(ha(6));
imshow(I3);
set(gca,... % ok
    'Position',[0.02,0.12,0.2,0.37]);

Im3Handle = imhandles(gca);
set(Im3Handle,'XData',[1,1100],'YData',[1,1100]);


%% Plot the TFR PSD for all patients

freq = MainGranger_Fig_Data.PSD.P38.TFR_baselinedPSD.freq;
time = MainGranger_Fig_Data.PSD.P38.TFR_baselinedPSD.time;
TFR_PSD_eeg = squeeze(MainGranger_Fig_Data.PSD.P38.TFR_baselinedPSD.powspctrm(15,:,:));
clim = [0 1];

axes(ha(2));
contourf(time,freq,(TFR_PSD_eeg)./2.5,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca, 'jet')
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4],'Position',[0.349912739965097,0.626530612244898,0.012041884816754,0.273469387755102]);
line([-5 -5],get(ha(2),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([-3 -3],get(ha(2),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([0 0],get(ha(2),'YLim'),'Color',[1 1 1],'LineStyle','--')
strYlab  = ylabel('Frequency (Hz)');%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')
set(gca,...
    'Position',[0.226406806282723,0.6265,0.1204,0.2735]);
set(strYlab,'Position',[-8.708405629448269,20.000030697615575,1]);
set(strXlab,'Position',[-2.0000,3.036,1]);
freq = MainGranger_Fig_Data.PSD.P37.freqAxis;
time = MainGranger_Fig_Data.PSD.P37.timeAxis;
TFR_PSD_ecog = MainGranger_Fig_Data.PSD.P37.Spectra;
clim = [0 1];

axes(ha(7));
contourf(time,freq,(TFR_PSD_ecog),100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca, 'jet')
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4],'Position',[0.347294938917976,0.197959183673469,0.012041884816754,0.287755102040816]);
line([-5 -5],get(ha(7),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([-3 -3],get(ha(7),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([0 0],get(ha(7),'YLim'),'Color',[1 1 1],'LineStyle','--')
strYlab = ylabel('Frequency (Hz)');%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')
set(gca,...
    'Position',[0.225130890052356,0.197959183673469,0.120418848167539,0.287755102040816]);
set(strYlab,'Position',[-8.894492802412614,20.000030697615575,1]);
set(strXlab,'Position',[-1.999996185302734,3.067989684526177,1]);

%% Plot the TFR PSDs Hipp
map = TFR_PSD.P38;
time = TFR_PSD.timeAxis;
freq = TFR_PSD.freqAxis;
clim  = TFR_PSD.clim;

axes(ha(3));
contourf(time,freq,map,100,'LineColor','none')
%Axes properties
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca, 'jet')
 % to achieve correct colorbar for negative red peaks use flipud(jet) and -TFR_psd_base
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4],'Position',[0.567190226876091,0.624489795918367,0.01151832460733,0.279591836734694]);

line([-5 -5],get(ha(3),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([-3 -3],get(ha(3),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([0 0],get(ha(3),'YLim'),'Color',[1 1 1],'LineStyle','--')
strYlab = ylabel('Frequency (Hz)');%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')
set(gca,...
    'Position',[0.449389179755673,0.624489795918363,0.114746945898778,0.279591836734698]);
set(strYlab,'Position',[-8.746060566468673,20.000030697615575,1])
set(strXlab,'Position',[-1.999996185302734,3.168666308303659,1])


axes(ha(8));

map = TFR_PSD.P37.map;
time = TFR_PSD.timeAxis;
freq = TFR_PSD.freqAxis;
clim  = TFR_PSD.clim;

contourf(time,freq,map,100,'LineColor','none')
%Axes properties
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca, 'jet')
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4],'Position',[0.56457242582897,0.197959183673469,0.011125654450262,0.291836734693877]);

line([-5 -5],get(ha(8),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([-3 -3],get(ha(8),'YLim'),'Color',[1 1 1],'LineStyle','--')
line([0 0],get(ha(8),'YLim'),'Color',[1 1 1],'LineStyle','--')
strYlab = ylabel('Frequency (Hz)');%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')
set(gca,...
    'Position',[0.450837696335079,0.197959183673469,0.110680628272251,0.29204081632653]);
set(strYlab,'Position',[-9.042509771795833,20.000030697615575,1]);
set(strXlab,'Position',[-2.0000,3.0081,1]);

%% Plot the Granger spectra

light_blue = [0.30,0.75,0.93];
orange       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',orange};

yshift_1= -0.02;
yshift_2 = -0.03;
gdata.Enc = MainGranger_Fig_Data.GrangerSpectra.P38.Enc;
gdata.Maint = MainGranger_Fig_Data.GrangerSpectra.P38.Maint;
significance_bar_enc = [9:16];
significance_bar_enc2 = [20:22]
significance_bar_maint = [9:16];


axes(ha(4));
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(1,:)*100,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(2,:)*100,'Color',Colors{2},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(1,:)*100,'Color',Colors{3},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(2,:)*100,'Color',Colors{4},'LineWidth',3);
xlim([4 30])
plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',3);
plot(significance_bar_enc2,yshift_1*100*ones(length(freq(significance_bar_enc2)),1)','Color',Colors{1},'LineWidth',3);
plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',3);
ylim([-0.05 0.2]*100)


strYlab = ylabel('Granger (%)','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')

set(ha(4),'FontSize',10);
set(ha(4),'YTick',[0 0.1 0.2]*100,'YTickLabel', [0 0.1 0.2]*100);
set(ha(4),'XTick',[4 10:10:30],'XTickLabel', [4 10:10:30]);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
set(gca,...
    'Position',[0.653141361256544,0.620408163265306,0.11910994764398,0.285714285714285]);
set(strYlab,'Position',[2.389606497897715,7.321440130472183,-1]);
set(strXlab,'Position',[10.954461674932862,-6.714285750474248,-1]);


gdata.Enc = MainGranger_Fig_Data.GrangerSpectra.P37.Enc;
gdata.Maint = MainGranger_Fig_Data.GrangerSpectra.P38.Maint;
significance_bar_enc = [10:18];
significance_bar_maint = [9:18];


axes(ha(9));
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(1,:)*100,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(2,:)*100,'Color',Colors{2},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(1,:)*100,'Color',Colors{3},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(2,:)*100,'Color',Colors{4},'LineWidth',3);
xlim([4 30])
plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',3);
plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',3);
ylim([-0.05 0.2]*100)


strYlab = ylabel('Granger (%)');%,'VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)');%,'VerticalAlignment','cap')

set(ha(9),'FontSize',10);
set(ha(9),'YTick',[0 0.1 0.2]*100,'YTickLabel', [0 0.1 0.2]*100);
set(ha(9),'XTick',[4 10:10:30],'XTickLabel', [4 10:10:30]);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
set(gca,...
    'Position',[0.653141361256544,0.2,0.119109947643981,0.289999999999999])
set(strYlab,'Position',[2.387490913527623,6.971843242645273,-1]);
set(strXlab,'Position',[10.954461674932862,-7.929577304863596,-1]);

%% Plot the TFR Granger

freq = MainGranger_Fig_Data.TFR_granger.P38.freq;
time = MainGranger_Fig_Data.TFR_granger.P38.time;
clim = MainGranger_Fig_Data.TFR_granger.P38.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P38.spectra;
cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap

% cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap
axes(ha(5));
contourf(time,freq,TFR_Grng*100,100,'LineColor','none')
Xlab = xlabel('Time (s)');
Ylab = ylabel('Frequency (Hz)');
ylim([5 30]);
cbh = colorbar
colormap(gca,cmap);
set(cbh, 'Ticks',[-20  20],'TickLabels', [-20  20], 'Position',[0.952006980802793,0.622448979591837,0.01086387434555,0.285714285714286]);
set(ha(5),'clim',clim,'yscale','log')
set(ha(5),'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30]);
set(ha(5),'FontSize',10);
set(ha(5),'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
set(ha(5),'box','off');
set(ha(5),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
line([-5 -5],get(ha(5),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([-3 -3],get(ha(5),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([0 0],get(ha(5),'YLim'),'Color',[0 0 0],'LineStyle','--')
set(gca,...
    'Position',[0.842931937172775,0.622448979591837,0.108638743455498,0.285714285714285]);

set(Ylab,'Position',[-8.628112491354885,12.247459177864897,1]);
set(Xlab,'Position',[-1.999996185302734,4.201246546230278,1]);



freq = MainGranger_Fig_Data.TFR_granger.P37.freq;
time = MainGranger_Fig_Data.TFR_granger.P37.time;
clim = MainGranger_Fig_Data.TFR_granger.P37.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P37.spectra;
% cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap
axes(ha(10));
contourf(time,freq,TFR_Grng*100,100,'LineColor','none')
Xlab = xlabel('Time (s)');
Ylab = ylabel('Frequency (Hz)');
ylim([5 30]);
cbh = colorbar
colormap(gca,cmap);
set(cbh, 'Ticks',[-20  20],'TickLabels', [-20  20],'Position',[0.955933682373473,0.2,0.01151832460733,0.285714285714286]);
set(ha(10),'clim',clim,'yscale','log')
set(ha(10),'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30]);
set(ha(10),'FontSize',10);
set(ha(10),'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
set(ha(10),'box','off');
set(ha(10),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
line([-5 -5],get(ha(10),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([-3 -3],get(ha(10),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([0 0],get(ha(10),'YLim'),'Color',[0 0 0],'LineStyle','--')
set(gca,...
    'Position',[0.837696335078534,0.2,0.115183246073298,0.285714285714285])
 
set(Ylab,'Position',[-8.296969630501486,12.247459177864897,1]);
set(Xlab,'Position',[-2.090905276211815,4.095073697658404,1]);


%% Axes properties

for i =1:length(ha)
    if i == 2||i == 3||i == 5||i == 7|| i == 8 ||i == 10
        set(ha(i),'box','off','FontSize',10.5,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025])
    else
        set(ha(i),'box','off','FontSize',10.5,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025])
    end
end



%% TextBoxes


textBox_FontSize = 16;

% Create textbox
annotation(fig,'textbox',...
    [0.0026178010471204 0.947491358871824 0.0392670149266408 0.0322966501116753],...
    'String',{'a'},...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.876963350785341 0.632151060215466 0.145287958115183 0.0370813389642958],...
    'Color',[1 1 1],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.00785340314136404 0.520000000000001 0.0392670149266414 0.0429556683917689],...
    'String','f',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.503926701570682 0.636232692868523 0.0772251289786472 0.0370813389642958],...
    'Color',[1 1 1],...
    'String','Hipp',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.0353403141361261 0.879382775868233 0.0890052333396144 0.0370813389642958],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.0379581151832465 0.448979591836737 0.0890052333396145 0.0438100771409022],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.287958115183248 0.636232692868521 0.077225128978647 0.037081338964296],...
    'Color',[1 1 1],...
    'String','EEG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.2749 0.2115 0.0916 0.0371],'Color',[1 1 1],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.503926701570682 0.205620447970563 0.077225128978647 0.037081338964296],...
    'Color',[1 1 1],...
    'String','Hipp',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.876963350785341 0.211742896950158 0.145287958115183 0.037081338964296],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.157068062827228 0.945450542545293 0.039267014926641 0.032296650111675],...
    'String',{'b',''},...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.382198952879591 0.949532175198355 0.0392670149266411 0.0322966501116756],...
    'String','c',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.592931937172797 0.953613807851415 0.039267014926641 0.032296650111675],...
    'String','d',...
    'FontSize',textBox_FontSize,...
     'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.778795811518335 0.953867689308302 0.0392670149266412 0.032296650111675],...
    'String','e',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.153141361256552 0.534105947290921 0.039267014926641 0.032296650111675],...
    'String','g',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.380890052356043 0.535092179104218 0.0392670149266411 0.032296650111675],...
    'String','h',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.598167539267017 0.530771312001399 0.039267014926641 0.032296650111675],...
    'String','i',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.781413612565449 0.530278196094747 0.0392670149266412 0.0322966501116748],...
    'String','j',...
    'FontSize',textBox_FontSize,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');


h5 = text(3.927560546181418,6.953785608651588,0.237200021743774,'\DeltaGranger (%)')
set(h5,'Rotation',90,'FontSize',12);

h6 = text(3.745742364363236,94.6452522693705,0.237200021743774,'\DeltaGranger (%)')
set(h6,'Rotation',90,'FontSize',12);


%% XTicks
for i = 2:size(ha,1)
    set(ha(i),'XTickLabel',[]);
end


% Create textbox
annotation(fig,'textbox',...
    [0.636125654450264 0.572469388789056 0.192408376963349 0.0551020397823682],...
    'String','4       10    20  30',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.219895287958116 0.14797959287069 0.134816753926701 0.0551020397823682],...
    'String',' -5   -3      0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.443717277486912 0.14797959287069 0.12696335078534 0.0551020397823683],...
    'String','-5    -3    0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.638743455497385 0.15002040919722 0.192408376963349 0.0551020397823682],...
    'String','4       10    20  30',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.825916230366494 0.150020409197221 0.134816753926701 0.0551020397823682],...
    'String',' -5   -3      0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');


% Create textbox
annotation(fig,'textbox',...
    [0.222513089005236 0.576551021442121 0.134816753926702 0.0551020397823684],...
    'String','-5    -3     0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.441099476439791 0.574510205115589 0.134816753926702 0.0551020397823684],...
    'String','-5    -3     0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.832460732984295 0.574510205115588 0.134816753926701 0.0551020397823682],...
    'String',' -5   -3     0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Save Figure
set(gcf,'color','white')
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Main Granger_Fig','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Main Granger_Fig','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Main Granger_Fig','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Main Granger_Fig.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Main Granger_Fig.tiff');
