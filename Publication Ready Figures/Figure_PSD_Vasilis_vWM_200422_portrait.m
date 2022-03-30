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
OmitScalp =1;
%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))
%% Figure Data

strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
strFilename = [strPaths.FigureData,'FigurePSD_Data'];
load(strFilename);
nGridChans =64;


%% Figure Axes
fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[4.8948,2.8581,15,20]);
h=gcf;


ha = tight_subplot(2,2,[0.05 0.05],[0.05 0.03],[0.05 0.05])

%%
if ~OmitScalp
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
    ylabel('Power 10*log10(\muV^2/Hz)')
    ylim([-30 30]);
    xlim([4 100])
    
    set(ha(1),'FontSize',10);
    set(ha(1),'XTick',[10 100],'XTickLabel', [10 100]);
    
    annotation('textbox',...
        [0.30514512383901 0.926760765693373 0.105833330594265 0.0383645825842396],...
        'String',{'EEG P3'},...
        'FontSize',11,...
        'FitBoxToText','on',...
        'EdgeColor','none');
    
    annotation('textbox',...
        [0.0084033613445378 0.965565350540506 0.0259978985413909 0.0435663617188686],...
        'String',{'a'},...
        'FontWeight','bold',...
        'FontSize',11,...
        'EdgeColor','none');
    
    set(ha(1),'Parent',fig,...
        'Position',[0.0708881578947369 0.535 0.332949122807018 0.435]);
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
end

%% Panel c
axes(ha(1));
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
yShift = -9;

semilogx(freqFix,Fix,'k','LineWidth',3);
hold on;
semilogx(freq,Enc,'g','LineWidth',3);
hold on;
semilogx(freq,Maint,'r','LineWidth',3);

ylim([ylimit(1) ylimit(2)])
xlim([4,100]);
color = 'm'
color2 = [0.75, 0.75, 0];
indMaxBand = Difference_Bar(Diff1,Percentile,freq,[ylimit(1) ylimit(2)],yShift,color);
indMaxBand2 = Difference_Bar(Diff2,Percentile,freq,[ylimit(1) ylimit(2)],yShift,color2);
ylabel('10*log_{10}PSD (\muV^2/Hz)')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
Lgd = legend('fixation','encoding','maintenance','box', 'off');

set(ha(1),'FontSize',10);
set(ha(1),'XTick',[10 100],'XTickLabel', [10 100]);
set(ha(1),'YTick',[-10:10:40],'YTickLabel', [-10:10:40]);
ylim([-10 40])
strAnnot = annotation('textbox',...
      [0.293834613767534,0.930289921337018,0.151675481014151,0.038359787610788],...
    'String',{'ECoG H3'},...
    'FontSize',11,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
set(ha(1),'Parent',fig,...
    'Position',[0.0864197530864199 0.734126984126984 0.358024691358025 0.232363315696648]);
axis square

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'-10','0','10','20','30','40'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end


ttLetter(1) = text(2,2,'a');
ttLetter(1).VerticalAlignment = 'bottom';
ttLetter(1).HorizontalAlignment = 'center';
ttLetter(1).FontWeight = 'bold';
ttLetter(1).Position(2) = ha(1).Position(4)+42;
ttLetter(1).Position(1) =ha(1).YLabel.Position(1)-0.3;

% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-3;
Lgd.Location = 'southwest';
% Lgd.Position = [0.269248158973721,0.868385159432256,0.148148146413621,0.070105818213609];
%% Panel d 
axes(ha(2));
Panel_d_dat = FigurePSD_Data.TFR_PSD.ECoG;
time = Panel_d_dat.time;
freq = Panel_d_dat.freq;
TFR_psd = Panel_d_dat.TFR_baselinedPSD;
clim = [0 1]%Panel_d_dat.clim;

contourf(time,freq,(TFR_psd-0.6)./3.5,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4]);
line([-5 -5],get(ha(2),'YLim'),'Color',[1 1 1])
line([-3 -3],get(ha(2),'YLim'),'Color',[1 1 1])
line([0 0],get(ha(2),'YLim'),'Color',[1 1 1])

ylabel('Frequency (Hz)','VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')

annotation('textbox',...
    [0.790732843137252,0.932264639236529,0.155222218121919,0.03836458258424],...
    'Color',[1 1 1],...
    'String',{'ECoG H3'},...
    'FontWeight','bold',...
    'FitBoxToText','on',...
    'FontSize',11,...
    'EdgeColor','none');

set(ha(2),'Parent',fig,...
    'Position',[0.505761574074074 0.736772486772487 0.414987870370371 0.231904596560846])

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'-5','-3','0'};
yTickLabel = {'4','10','20','40','100'};

for k = 1:length(xTick)
    text(xTick(k),3.5,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','cap')
end

ttLetter(2) = text(2,2,'b');
ttLetter(2).VerticalAlignment = 'bottom';
ttLetter(2).HorizontalAlignment = 'center';
ttLetter(2).FontWeight = 'bold';
ttLetter(2).Position(1) =ha(2).YLabel.Position(1)-0.03;
ttLetter(2).Position(2) = ha(2).Position(4)+120
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = 3.5-0.85;
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
%% Panel e
axes(ha(3));
Panel_e_dat = FigurePSD_Data.PSD_Spectra.Hipp;
freq = Panel_e_dat.freq;
freqFix = Panel_e_dat.freqFix;
Maint = Panel_e_dat.Maint;
Enc = Panel_e_dat.Enc;
Fix = Panel_e_dat.Fix;

semilogx(freqFix,Fix,'k','LineWidth',3);
hold on;
semilogx(freq,Enc,'g','LineWidth',3);
hold on;
semilogx(freq,Maint,'r','LineWidth',3);
ylabel('10*log_{10}PSD (\muV^2/Hz)');
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
% xlabel('Frequency (Hz)','VerticalAlignment','baseline');
indMaxBand =[41 52];
color = 'm';
yShift = -19;
semilogx(freq(indMaxBand),yShift*ones(length(freq(indMaxBand)),1)','Color',color,'LineWidth',3)
xlim([4 100])
ylim([-20 20])


set(ha(3),'FontSize',10);
set(ha(3),'XTick',[10 100],'XTickLabel', [10 100]);
set(ha(3),'YTick',[-20:10:20],'YTickLabel',[-20:10:20]);


annotation('textbox',...
    [0.343805733696988,0.538667015693373,0.095249997687009,0.03836458258424],...
    'String',{'Hipp'},...
    'FontSize',11,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
set(ha(3),'Parent',fig,...
    'Position',[0.0864197530864199 0.3441 0.358024691358025 0.232363315696648]);

axis square

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'-20','-10','0','10','20'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end

ttLetter(3) = text(2,2,'c');
ttLetter(3).VerticalAlignment = 'bottom';
ttLetter(3).HorizontalAlignment = 'center';
ttLetter(3).FontWeight = 'bold';
ttLetter(3).Position(2) = ha(3).Position(4)+22;
ttLetter(3).Position(1) =ttLetter(1).Position(1);

yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);
strXlab.Position(1) = 20;
strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-3;
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
%% Panel f
axes(ha(4));
Panel_f_dat = FigurePSD_Data.TFR_PSD.iEEG_Hipp;
time = Panel_f_dat.time;
freq = Panel_f_dat.freq;
TFR_psd =Panel_f_dat.TFRbaselinedPSD;
clim = [-1 0]% Panel_f_dat.clim;

contourf(time,freq,(TFR_psd-3.5)./3.5,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar('Ticks',[-1:1:4],'TickLabels',[-1:1:4]);
ylabel('Frequency (Hz)','VerticalAlignment','middle')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')
% title('TFR Power - iEEG Hipp','FontSize',13);
line([-5 -5],get(ha(4),'YLim'),'Color',[1 1 1])
line([-3 -3],get(ha(4),'YLim'),'Color',[1 1 1])
line([0 0],get(ha(4),'YLim'),'Color',[1 1 1])
annotation('textbox',...
    [0.844816176470586,0.542687754066322,0.100541664195971,0.03836458258424],...
    'Color',[1 1 1],...
    'String',{'Hipp'},...
    'FontWeight','bold',...
    'FitBoxToText','on',...
    'FontSize',11,...
    'EdgeColor','none');

set(ha(4),'Parent',fig,...
    'Position',[0.505761574074074  0.3441  0.414987870370371 0.231904596560846])



xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'-5','-3','0'};
yTickLabel = {'4','10','20','40','100'};

for k = 1:length(xTick)
    text(xTick(k),3.5,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','cap')
end

ttLetter(4) = text(2,2,'d');
ttLetter(4).VerticalAlignment = 'bottom';
ttLetter(4).HorizontalAlignment = 'center';
ttLetter(4).FontWeight = 'bold';
ttLetter(4).Position(1) =ha(4).YLabel.Position(1)-0.03;
ttLetter(4).Position(2) = ha(4).Position(4)+125
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = 3.5-0.85;
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);



% Create textbox
% Create textbox
annotation('textbox',...
    [0.557388333333334 0.729220000551965 0.119944441395501 0.0317499994480361],...
    'Color',[1 1 1],...
    'String',{'encoding'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.502707777777778 0.729220000551965 0.0634999986332324 0.0317499994480361],...
    'Color',[1 1 1],...
    'String',{'fix.'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.506235555555557 0.336313750551966 0.0634999986332324 0.0317499994480362],...
    'Color',[1 1 1],...
    'String',{'fix.'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.557388333333335 0.336313750551966 0.119944441395501 0.0317499994480362],...
    'Color',[1 1 1],...
    'String',{'encoding'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.666749444444446 0.729220000551966 0.151694440449278 0.0317499994480361],...
    'Color',[1 1 1],...
    'String',{'maintenance'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.668513333333336 0.336313750551966 0.151694440449278 0.0317499994480362],...
    'Color',[1 1 1],...
    'String',{'maintenance'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');


%% Save Fig in Different formats
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');
set(ha(1:size(ha)),'box','off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PSD\vWM grid PSD','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PSD\vWM grid PSD','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PSD\vWM grid PSD','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PSD\vWM grid PSD.fig');

