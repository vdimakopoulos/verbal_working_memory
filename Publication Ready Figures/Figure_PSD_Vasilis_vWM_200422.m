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

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
%% Figure Data

strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
strFilename = [strPaths.FigureData,'FigurePSD_Data'];
load(strFilename);
nGridChans =64;


%% Figure Axes
fig = figure('units','normalized','outerposition',[0 0 1 1],'PaperSize',[20 11])
fig.Units = 'centimeters';
% fig = figure;
% set(fig,'PaperPositionMode','auto');
% set(fig,'PaperOrientation','landscape');
% set(fig,'PaperUnits','normalized');
% 
% 
% 
% fig.Position =  [40 40 1200 800];

ha = tight_subplot(3,2,[0.05 0.05],[0.05 0.03],[0.05 0.05])
% print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PSD','-dpdf');



%% Panel a
axes(ha(1));
Panel_a_dat = FigurePSD_Data.PSD_Spectra.EEG_P3;
freq = Panel_a_dat.freq;
freqFix = Panel_a_dat.freqFix;
Maint = Panel_a_dat.Maint;
Enc = Panel_a_dat.Enc;
Fix = Panel_a_dat.Fix;

semilogx(freq,Maint,'k','LineWidth',3);
hold on;
semilogx(freq,Enc,'g','LineWidth',3);
hold on;
semilogx(freqFix,Fix,'r','LineWidth',3);
ylabel('Power 10*log10(\muV^2/Hz)')
ylim([-30 30]);
xlim([4 100])

set(ha(1),'FontSize',13);
set(ha(1),'XTick',[10 100],'XTickLabel', [10 100]);

annotation('textbox',...
    [0.417941176470588 0.936446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'EEG P3'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.965565350540506 0.0259978985413909 0.0435663617188686],...
    'String',{'a'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
%% Panel b
axes(ha(2))
Panel_b_dat = FigurePSD_Data.ScalpTopo;
datavector = Panel_b_dat.datavector;
freqBand = Panel_b_dat.freqBand;
freqAxis = Panel_b_dat.freqAxis;
fr_maint_sclp_recons = Panel_b_dat.fr_maint_sclp_recons;
Scalp_topoplot(strPaths.Toolboxes,freqBand,freqAxis,datavector,fr_maint_sclp_recons,4,1,1)
colorbar('Ticks',[0 20 40],'TickLabels',[0 20 40],'FontSize',13)

annotation('textbox',...
    [0.484033613445378  0.965565350540506 0.0259978985413909 0.0435663617188686],...
    'String',{'b'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

annotation('textbox',...
    [0.544033613445378 0.936446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'EEG P3 [6 9] Hz'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     

%% Panel c
axes(ha(3));
Panel_c_dat = FigurePSD_Data.PSD_Spectra.ECoG;
freq = Panel_c_dat.freq;
freqFix = Panel_c_dat.freqFix;
Maint = Panel_c_dat.Maint;
Enc = Panel_c_dat.Enc;
Fix = Panel_c_dat.Fix;
Diff1 = Panel_c_dat.SignificanceMaint_Enc;
Diff2 = Panel_c_dat.SignificanceEnc_Maint;
Percentile = Panel_c_dat.Percentile;
ylimit = [-15 40];
yShift = -13.5;

semilogx(freq,Maint,'k','LineWidth',3);
hold on;
semilogx(freq,Enc,'g','LineWidth',3);
hold on;
semilogx(freqFix,Fix,'r','LineWidth',3);
ylim([ylimit(1) ylimit(2)])
xlim([4,100]);
color = 'm'
color2 = [0.75, 0.75, 0];
indMaxBand = Difference_Bar(Diff1,Percentile,freq,[ylimit(1) ylimit(2)],yShift,color);
indMaxBand2 = Difference_Bar(Diff2,Percentile,freq,[ylimit(1) ylimit(2)],yShift,color2);
ylabel('Power 10*log10(\muV^2/Hz)')


set(ha(3),'FontSize',13);
set(ha(3),'XTick',[10 100],'XTickLabel', [10 100]);


annotation('textbox',...
    [0.417941176470588 0.606446809023954  0.0861344537815125 0.0273556225869914],...
    'String',{'ECoG H3'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.62794 0.0259978985413909 0.0435663617188686],...
    'String',{'c'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');


%% Panel d 
axes(ha(4));
Panel_d_dat = FigurePSD_Data.TFR_PSD.ECoG;
time = Panel_d_dat.time;
freq = Panel_d_dat.freq;
TFR_psd = Panel_d_dat.TFR_baselinedPSD;
clim = Panel_d_dat.clim;

contourf(time,freq,TFR_psd,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',13)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar
% title('TFR Power - ECoG H3','FontSize',13)  ;

annotation('textbox',...
    [0.484033613445378 0.62794 0.0259978985413909 0.0435663617188686],...
    'String',{'d'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

annotation('textbox',...
    [0.855441176470585 0.600367781879698 0.0569852925758778 0.0344478209571877],...
    'Color',[1 1 1],...
    'String',{'ECoG H3'},...
    'FontWeight','bold',...
    'FontSize',13,...
    'EdgeColor','none');

%% Panel e
axes(ha(5));
Panel_e_dat = FigurePSD_Data.PSD_Spectra.Hipp;
freq = Panel_e_dat.freq;
freqFix = Panel_e_dat.freqFix;
Maint = Panel_e_dat.Maint;
Enc = Panel_e_dat.Enc;
Fix = Panel_e_dat.Fix;

semilogx(freq,Maint,'k','LineWidth',3);
hold on;
ylabel('Power 10*log10(\muV^2/Hz)');
% xlabel('Frequency (Hz)','VerticalAlignment','baseline');
semilogx(freq,Enc,'g','LineWidth',3);
hold on;
semilogx(freqFix,Fix,'r','LineWidth',3);
xlim([4 100])
ylim([-20 20])

set(ha(5),'FontSize',13);
set(ha(5),'XTick',[10 100],'XTickLabel', [10 100]);


annotation('textbox',...
    [0.417941176470588 0.286446809023954  0.0861344537815125 0.0273556225869914],...
    'String',{'iEEG Hipp'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.30794 0.0259978985413909 0.0435663617188686],...
    'String',{'e'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');


%% Panel f
axes(ha(6));
Panel_f_dat = FigurePSD_Data.TFR_PSD.iEEG_Hipp;
time = Panel_f_dat.time;
freq = Panel_f_dat.freq;
TFR_psd =Panel_f_dat.TFRbaselinedPSD;
clim = Panel_f_dat.clim;

contourf(time,freq,TFR_psd,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',13)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar
% title('TFR Power - iEEG Hipp','FontSize',13);

annotation('textbox',...
    [0.484033613445378 0.30794 0.0259978985413909 0.0435663617188686],...
    'String',{'f'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

annotation('textbox',...
    [0.863319327731085 0.278179332031678 0.0569852925758778 0.0344478209571877],...
    'Color',[1 1 1],...
    'String',{'iEEG Hipp'},...
    'FontWeight','bold',...
    'FontSize',13,...
    'EdgeColor','none');

%% Save Fig in Different formats
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');
set(ha(1:size(ha)),'box','off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PSD','-dpdf');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PSD','-dpng');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PSD','-depsc');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PSD.fig');