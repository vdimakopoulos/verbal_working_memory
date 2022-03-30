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
strFilename = [strPaths.FigureData,'FigureMaintData'];
load(strFilename);

nGridChans =64;

%% Figure Axes
% fig = figure;
% fig = figure('units','normalized','outerposition',[0 0 1 1],'PaperSize',[20 11])
% fig.Units = 'centimeters';
% % fig = figure;
% % set(fig,'PaperPositionMode','auto');
% % set(fig,'PaperOrientation','landscape');
% % set(fig,'PaperUnits','normalized');
% % 
% % 
% % 
% % fig.Position =  [40 40 1200 800];
fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[4.8948,2.8581,19,20]);
h=gcf;


ha = tight_subplot(4,3,[0.05 0.05],[0.05 0.03],[0.05 0.05])
set(ha(1:size(ha,1)),'Fontsize',20)



%% Panel a
axes(ha(1));
Panel_a_Dat = FigureMaint_Data.PLV_Spectra_Scalp_Hipp;
freqHi = 30:5:100;
iSS = 4;

nChannelPairs = Panel_a_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_a_Dat.strChannelNameList(nChannelPairs(:,1))','_AHL2-_AHL3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_a_Dat.strChannelNameList(nChannelPairs(:,2))','_P3')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

Vars = Panel_a_Dat.SignificanceBars.Scalp_Hipp;
Maint_Fix_Diff = abs(Panel_a_Dat.PLV_PairsMaint{nPair,iSS})-abs(Panel_a_Dat.PLV_PairsFix{nPair,iSS});
semilogx(Panel_a_Dat.freqAxis,abs(Panel_a_Dat.PLV_PairsFix{nPair,iSS}),'k','LineWidth',3)
hold on;
semilogx(Panel_a_Dat.freqAxis,abs(Panel_a_Dat.PLV_PairsMaint{nPair,iSS}),'r','LineWidth',3)
hold on
indMaxBandMaintFix = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_a_Dat.freqAxis,[0 1], 0.02,'r');
hold on;
semilogx(freqHi,abs(Panel_a_Dat.PLV_Pairs_FixHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_a_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'r','LineWidth',3);
ylim([0,1])
xlim([4 100])
Lgd = legend('fixation','maintenance','box','off')

set(ha(1),'Fontsize',10)
set(ha(1),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(1),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);
strYlab = ylabel('PLV','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
strAnnot = annotation('textbox',...
    [0.12941176470588 0.936446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'Hipp - EEG P3'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

set(ha(1),'Parent',fig,...
    'Position',[0.05 0.7775 0.227533039647577 0.1925]);


xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end

ttLetter(1) = text(2,2,'a');
ttLetter(1).VerticalAlignment = 'bottom';
ttLetter(1).HorizontalAlignment = 'center';
ttLetter(1).FontWeight = 'bold';
ttLetter(1).Position(2) = ha(1).Position(4)+0.8;
ttLetter(1).Position(1) =ha(1).YLabel.Position(1)-0.4;

axis square
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-0.07;
strYlab.Position(1)= 2.1373;
Lgd.Position(2) = strAnnot.Position(2)-0.06;
Lgd.Position(1) = 0.100401309379775;

set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel b
axes(ha(2));

Panel_b_Dat = FigureMaint_Data.TFR_PLV;

clim = [0.2 0.6];
iSS_To_Plot = 4;

contourf(Panel_b_Dat.Time,Panel_b_Dat.FreqTFR,abs(Panel_b_Dat.PLV{iSS_To_Plot}),100,'LineColor','none')
%Axes properties
set(gca,'ytick', [5,11,23,47,97],'YTickLabel',[4,10,20,40,100]) %background color of grid

set(gca,'clim',clim,'yscale','log')
set(gca,'color',[0.01 0.01 0.56]);
set(gcf,'color','w');
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar
line([-5 -5],get(ha(2),'YLim'),'Color',[1 1 1])
line([-3 -3],get(ha(2),'YLim'),'Color',[1 1 1])
line([0 0],get(ha(2),'YLim'),'Color',[1 1 1])
ylabel('Frequency (Hz)','VerticalAlignment','bottom')
strXlab = xlabel('Time (s)','VerticalAlignment','cap')

annotation('textbox',...
    [0.487307883679589,0.928083682281187,0.161535083320188,0.035718749329758],...
    'Color',[1 1 1],...
    'String',' Hipp - EEG P3',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.717839352883604 0.712 0.0259978985413909 0.0435663617188686],...
    'String',{'f'},...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.455361403508776 0.776845000551966 0.119758768775745 0.0317499994480361],...
    'Color',[1 1 1],...
    'String',{'maintenance'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.374593859649126 0.776845000551966 0.0946929800490801 0.0317499994480361],...
    'Color',[1 1 1],...
    'String',{'encoding'},...
    'FontWeight','bold',...
    'FontSize',8,...
    'EdgeColor','none');

% Create textbox
annotation('textbox',...
    [0.338387719298249 0.776845000551966 0.0362061403508738 0.0317499994480361],...
    'Color',[1 1 1],...
    'String','fix.',...
    'FontWeight','bold',...
    'FontSize',8,...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(ha(2),'Parent',fig,...
    'Position',[0.343483796296296 0.7775 0.289571759259259 0.1925]);


xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'-5','-3','0'};
yTickLabel = {'4','10','20','40','100'};

for k = 1:length(xTick)
    text(xTick(k),3.5,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end
% yTick(1)-0.05*(yTick(end)-yTick(1))
ttLetter(2) = text(2,2,'b');
ttLetter(2).VerticalAlignment = 'bottom';
ttLetter(2).HorizontalAlignment = 'center';
ttLetter(2).FontWeight = 'bold';
ttLetter(2).Position(1) =ha(2).YLabel.Position(1)-0.2;
ttLetter(2).Position(2) = ha(2).Position(4)+100
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = 3.5-0.7;
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel c
axes(ha(3));
Panel_c_Dat = FigureMaint_Data.ScalpTopos.Sclp_Hipp.Freq_6_7;
FreqBand = Panel_c_Dat.freqBand;
FreqAxis = Panel_c_Dat.freqAxis;
dataBipolar = Panel_c_Dat.dataBipolarSclp;
datavector = Panel_c_Dat.datavector;
nPairs = Panel_c_Dat.nChannelPairs;
Scalp_topoplot(strPaths.Toolboxes,FreqBand,FreqAxis,datavector,dataBipolar,4,nPairs,1);


annotation('textbox',...
    [0.717839352883604 ttLetter(1).Position(2)-0.04 0.0259978985413909 0.0435663617188685],...
    'String',{'c'},...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);

annotation('textbox',...
    [0.810065052336529,0.933459283730825,0.193563591137594,0.035718749329758],...
    'String',{'Hipp - EEG [6 7] Hz'},...
    'FitBoxToText','on',...
    'EdgeColor','none');

set(ha(3),'Parent',fig,...
    'Position',[0.787543421052632 0.7775 0.138345452002634 0.1925]);
% 
% ttLetter(3) = text(2,2,'c');
% ttLetter(3).VerticalAlignment = 'bottom';
% ttLetter(3).HorizontalAlignment = 'center';
% ttLetter(3).FontWeight = 'bold';
% ttLetter(3).Position(2) =ha(3).Position(1)+100

%% Panel d
axes(ha(4));
axis square

iSS_to_plot =4;
Panel_d_Dat = FigureMaint_Data.PLV_Spectra_Scalp_Grid;
freqHi = 30:5:100;
iSS_to_plot = 4;


nChannelPairs = Panel_d_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_d_Dat.strChannelNameList(nChannelPairs(:,1))','_P3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_d_Dat.strChannelNameList(nChannelPairs(:,2))','_GL_D6')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

for i = 1
    for j = nPair
        s = nChannelPairs(find(nChannelPairs== nChannelPairs((j+(i-1)*nGridChans),1)),1);
        s_loc = s(1);
        semilogx(Panel_d_Dat.freqAxis,abs(Panel_d_Dat.PLV_PairsFix{(j+(s_loc-1)*nGridChans),...
            iSS_to_plot}),'k','LineWidth',3);
        hold on;
        semilogx(Panel_d_Dat.freqAxis,abs(Panel_d_Dat.PLV_PairsMaint{(j+(s_loc-1)*nGridChans),... 
            iSS_to_plot}),'r','LineWidth',3);
        semilogx(freqHi,abs(Panel_d_Dat.PLV_Pairs_FixHi{nPair,iSS_to_plot}),'k','LineWidth',3);
        semilogx(freqHi,abs(Panel_d_Dat.PLV_Pairs_MaintHi{nPair,iSS_to_plot}),'r','LineWidth',3);

        ylim([0 1])
        xlim([4 100]);
        hold on;
        Maint_Fix_Diff = abs(Panel_d_Dat.PLV_PairsMaint{(j+(s_loc-1)*nGridChans),iSS_to_plot}) ...
            - abs(Panel_d_Dat.PLV_PairsFix{(j+(s_loc-1)*nGridChans),iSS_to_plot});
        Vars = Panel_d_Dat.SignificanceBars.Scalp_Grid;
        indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_d_Dat.freqAxis,[0 1],0.02,'r');
    end
  
end
strYlab = ylabel('PLV','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
     
set(ha(4),'Fontsize',10)
set(ha(4),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(4),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.09111176470588 0.696446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'ECoG D6 - EEG P3'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');
set(ha(4), 'Position',[0.05 0.535 0.227533039647577 0.1925]);

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end

ttLetter(3) = text(2,2,'d');
ttLetter(3).VerticalAlignment = 'bottom';
ttLetter(3).HorizontalAlignment = 'center';
ttLetter(3).FontWeight = 'bold';
ttLetter(3).Position(2) = ha(3).Position(4)+0.8;
ttLetter(3).Position(1) =ttLetter(1).Position(1);%ha(3).YLabel.Position(1);

axis square
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-0.053;
strYlab.Position(1)= 2.1373;
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel e

axes(ha(5));
Panel_e_Dat = FigureMaint_Data.PLV_maps.ScalpGrid.Freq6_9;
imagesc(Panel_e_Dat.PLV_In_Band_To_Plot,[0 0.7])
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);
colormap jet
set(ha(5),'FontSize',10)
% title('ECoG - EEG P3 [6 9] Hz','FontSize',13);

annotation('textbox',...
    [0.413130252100642 0.696446809023954 0.108193277311122 0.0314083073699607],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','ECoG - EEG P3 [6 9] Hz',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

set(ha(5),'Position',[0.343483796296296 0.535 0.289571759259259 0.1925]);

% yTick(1)-0.05*(yTick(end)-yTick(1))
ttLetter(4) = text(2,2,'e');
ttLetter(4).VerticalAlignment = 'bottom';
ttLetter(4).HorizontalAlignment = 'center';
ttLetter(4).FontWeight = 'bold';
% ttLetter(4).Position(1) = ttLetter(2).Position(1);
ttLetter(4).Position(2) = ttLetter(3).Position(2)+7.5;
ttLetter(4).Position(1) = ttLetter(4).Position(1)-3;

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
yTickLabel = {'1','2','3','4','5','6','7','8'};
xTickLabel = {'A','B','C','D','E','F','G','H'};

for k = 1:length(xTick)
    text(xTick(k),0.1,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel f
axes(ha(6));
Panel_f_Dat = FigureMaint_Data.ScalpTopos.Sclp_Grid.Freq_6_9;
datavector = Panel_f_Dat.datavector;
FreqBand = Panel_f_Dat.freqBand;
freqAxis = Panel_f_Dat.freqAxis;
dataBipolar = Panel_f_Dat.dataBipolarSclp;
nPairs = Panel_f_Dat.nChannelPairs;

Scalp_topoplot(strPaths.Toolboxes,FreqBand,freqAxis,datavector,dataBipolar,4,nPairs,1)



annotation('textbox',...
    [0.717839352883604 0.712 0.0259978985413909 0.0435663617188686],...
    'String',{'f'},...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);

annotation('textbox',...
    [0.765887420757581,0.692136367064156,0.235339905682036,0.035718749329758],...
    'String',{'ECoG D6 - EEG [6 9] Hz'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

set(ha(6),'Position',[0.787543421052632  0.535 0.138345452002634 0.1925]);

%% Panel g
axes(ha(7))
Panel_g_Dat = FigureMaint_Data.PLV_Spectra_Grid_Hipp;
strPlotColors = {'k','g','b','r','c'};
iSS = 4;
freqHi =[30:5:100];

nChannelPairs = Panel_g_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,1))','_AHL2-_AHL3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,2))','_GL_H2')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

semilogx(Panel_g_Dat.freqAxis,abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS}),strPlotColors{1},'LineWidth',3);
hold on;
semilogx(Panel_g_Dat.freqAxis,abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS}),strPlotColors{4},'LineWidth',3);
hold on;
semilogx(freqHi,abs(Panel_g_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_g_Dat.PLV_Pairs_FixHi{nPair,iSS}),'r','LineWidth',3)

Vars = Panel_g_Dat.SignificanceBars.GridHipp;
Maint_Fix_Diff = abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS});
indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'r');
Fix_Maint_Diff =abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'k');


ylim([0,1]);
xlim([4 100])
strYlab = ylabel('PLV','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')

set(ha(7),'Fontsize',10)
set(ha(7),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(7),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.117941176470588 0.4564 0.0861344537815125 0.0273556225869914],...
    'String',{'Hipp - ECoG H2'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
% annotation('textbox',...
%     [0.0010033613445378 0.4808 0.0259978985413909 0.0435663617188686],...
%     'String',{'g'},...
%     'FontWeight','bold',...
%     'FontSize',11,...
%     'EdgeColor','none');
    

set(ha(7), 'Position',[0.05 0.2950 0.227533039647577 0.1925]);


xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end

ttLetter(5) = text(2,2,'g');
ttLetter(5).VerticalAlignment = 'bottom';
ttLetter(5).HorizontalAlignment = 'center';
ttLetter(5).FontWeight = 'bold';
ttLetter(5).Position(2) = ha(7).Position(4)+0.8;
ttLetter(5).Position(1) =ttLetter(1).Position(1);%ha(3).YLabel.Position(1);

axis square
% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-0.053;
strYlab.Position(1)= 2.1373;
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel h

axes(ha(8))
Panel_h_Dat = FigureMaint_Data.PLV_maps.GridHipp.Freq9_10.PLV_In_Band_To_Plot;

imagesc(Panel_h_Dat,[0 0.7])
ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);
colormap jet

set(ha(8),'FontSize',10)

annotation('textbox',...
    [0.4313130252100642 0.4564 0.120273109243694 0.031408307369961],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','Hipp - ECoG [9 10] Hz',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');


set(ha(8), 'Position',[0.343483796296296 0.2950 0.289571759259259 0.1925]);



% yTick(1)-0.05*(yTick(end)-yTick(1))
ttLetter(6) = text(2,2,'h');
ttLetter(6).VerticalAlignment = 'bottom';
ttLetter(6).HorizontalAlignment = 'center';
ttLetter(6).FontWeight = 'bold';
% ttLetter(4).Position(1) = ttLetter(2).Position(1);
ttLetter(6).Position(2) = ttLetter(3).Position(2)+7.5;
ttLetter(6).Position(1) = ttLetter(4).Position(1);

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
yTickLabel = {'1','2','3','4','5','6','7','8'};
xTickLabel = {'A','B','C','D','E','F','G','H'};

for k = 1:length(xTick)
    text(xTick(k),0.1,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel i
axes(ha(9));
Panel_i_Dat = FigureMaint_Data.ScalpTopos.Grid_Hipp.Freq9_10;
FreqBand_Sclp = Panel_i_Dat.freqBand;
freqAxis_sclp = Panel_i_Dat.freqAxis;
datavector = Panel_i_Dat.datavector;
dataBipolar = Panel_i_Dat.dataBipolarSclp;
nPairs = Panel_i_Dat.nChannelPairs;
Scalp_topoplot(strPaths.Toolboxes,FreqBand_Sclp,freqAxis_sclp,datavector,dataBipolar,4,nPairs,1);


annotation('textbox',...
    [0.717839352883604 0.473445833333333 0.0259978985413909 0.0435663617188686],...
    'String',{'i'},...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);

annotation('textbox',...
    [0.753354526020739,0.452089558040203,0.245087712409073,0.035718749329758],...
    'String',{'ECoG H2 - EEG [9 10] Hz'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');
set(ha(9),'Position',[0.787543421052632 0.2950 0.138345452002634 0.1925]);



%% Panel j
axes(ha(10))
axis square

Panel_j_Dat = FigureMaint_Data.PLV_Spectra_Grid_Hipp;
nChannelPairs = Panel_j_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,1))','_AHL2-_AHL3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,2))','_GL_C2')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

strPlotColors = {'k','g','b','r','c'};
iSS = 4;
freqHi = [30:5:100];
semilogx(freqAxis,abs(Panel_j_Dat.PLV_PairsFix{nPair,iSS}),strPlotColors{1},'LineWidth',3);
hold on;
semilogx(freqAxis,abs(Panel_j_Dat.PLV_PairsMaint{nPair,iSS}),strPlotColors{4},'LineWidth',3);
hold on;
semilogx(freqHi,abs(Panel_j_Dat.PLV_Pairs_FixHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_j_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'r','LineWidth',3);

Vars = Panel_j_Dat.SignificanceBars.GridHipp2;
Maint_Fix_Diff = abs(Panel_j_Dat.PLV_PairsMaint{nPair,iSS})-abs(Panel_j_Dat.PLV_PairsFix{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'k');
indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'r');
Fix_Maint_Diff =abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'r');
strYlab = ylabel('PLV','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
ylim([0,1]);
xlim([4 100])

set(ha(10),'Fontsize',10)
set(ha(10),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(10),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.117941176470588 0.2164 0.0861344537815125 0.0273556225869914],...
    'String',{'Hipp - ECoG C2'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

    
set(ha(10), 'Position',[0.05 0.0550 0.227533039647577 0.1925]);

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','100'};
yTickLabel = {'0','0.5','1'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end

ttLetter(7) = text(2,2,'j');
ttLetter(7).VerticalAlignment = 'bottom';
ttLetter(7).HorizontalAlignment = 'center';
ttLetter(7).FontWeight = 'bold';
ttLetter(7).Position(2) = ha(10).Position(4)+0.8;
ttLetter(7).Position(1) =ttLetter(1).Position(1);%ha(3).YLabel.Position(1);

% yTickLabel= cellstr(char(num2str([0 0.5 1]')));
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-0.07;
strXlab.Position(1) = strXlab.Position(1)+10;
strYlab.Position(1)= 2.1373;

axis square
set(gca,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel k
axes(ha(11));
Panel_k_Dat = FigureMaint_Data.PLV_maps.GridHipp.Freq16_29.PLV_In_Band_To_Plot;

imagesc(Panel_k_Dat,[0 0.7])
ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',10);
colormap jet

% annotation('textbox',...
%     [0.284033613445378 0.2418 0.0259978985413909 0.0435663617188686],...
%     'String',{'k'},...
%     'FontWeight','bold',...
%     'FontSize',11,...
%     'EdgeColor','none');
set(ha(11),'Fontsize',10)

% title('sEEG Hipp - ECoG [16 29] Hz','FontSize',13);

annotation('textbox',...
    [0.419970175438597,0.212089558040203,0.23255481804574,0.035718749329758],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','Hipp - ECoG [16 29] Hz',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');

set(ha(11), 'Position',[0.343483796296296 0.0550 0.289571759259259 0.1925]);




% yTick(1)-0.05*(yTick(end)-yTick(1))
ttLetter(8) = text(2,2,'k');
ttLetter(8).VerticalAlignment = 'bottom';
ttLetter(8).HorizontalAlignment = 'center';
ttLetter(8).FontWeight = 'bold';
% ttLetter(4).Position(1) = ttLetter(2).Position(1);
ttLetter(8).Position(2) = ttLetter(3).Position(2)+7.5;
ttLetter(8).Position(1) = ttLetter(4).Position(1);

xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
yTickLabel = {'1','2','3','4','5','6','7','8'};
xTickLabel = {'A','B','C','D','E','F','G','H'};

for k = 1:length(xTick)
    text(xTick(k),0.1,xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','top')
end
yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);
set(gca,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);

%% Panel l 
axes(ha(12));
figure;
Panel_l_Dat = FigureMaint_Data.Granger;
Grngr_Enc = Panel_l_Dat.Grng_Encoding{1};
Grngr_Maint = Panel_l_Dat.Grng_Maint{1};
indMaxBandEnc = Panel_l_Dat.significance.Enc_CortexHipp;
indMaxBandMaint = Panel_l_Dat.significance.Maint_HippCortex;
yshift_1= 0.01;
yshift_2 = 0.005;
freq = Panel_l_Dat.freq;
dark_Green = [0.02 0.47 0.04];
orange       = [1 0.45 0.45];
Colors     = {dark_Green,'g','r',orange};


hold on;
semilogx(freq,Grngr_Enc.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
semilogx(freq,Grngr_Maint.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
semilogx(freq,Grngr_Maint.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
semilogx(freq,Grngr_Enc.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
plot(indMaxBandEnc,yshift_1*ones(length(freq(indMaxBandEnc)),1)','Color',Colors{1},'LineWidth',3);
plot(indMaxBandMaint,yshift_2*ones(length(freq(indMaxBandMaint)),1)','Color',Colors{3},'LineWidth',3);

xlim([4 30]);
ylim([0 0.2]);
strYlab = ylabel('Granger causality','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')
% 
% lgd = legend('ECoG C2- sEEG Hipp enc','sEEG Hipp - ECoG C2 enc','sEEG Hipp - ECoG C2 maint','ECoG C2 - sEEG Hipp maint');
% lgd.Location = 'northwest';
% lgd.FontSize =10;

set(ha(12),'FontSize',10);
set(ha(12),'YTick',[0 0.1 0.2],'YTickLabel', [0 0.1 0.2]);
set(ha(12),'XTick',[10:10:30],'XTickLabel', [10:10:30]);

annotation('textbox',...
    [0.717839352883604 0.2418 0.0259978985413909 0.0435663617188686],...
    'String',{'l'},...
    'FontWeight','bold',...
    'FontSize',10,...
    'EdgeColor','none');

annotation('textbox',...
    [0.7861,0.208,0.1977,0.0357],...
    'String',{'  Hipp - ECoG C2'},...
    'FontSize',10,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     

set(ha(12),'Fontsize',10)
set(ha(12),'Position',[0.770076315789474 0.055 0.208881578947369 0.1925]);



xTick = get(gca,'xtick');
yTick = get(gca,'ytick');
set(gca,'xticklabel',[])
xTickLabel = {'10','20','30'};
yTickLabel = {'0','0.1','0.2'};

for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.05*(yTick(end)-yTick(1)),xTickLabel{k},'HorizontalAlignment','center','VerticalAlignment','middle')
end

yTickLabel =  cellfun(@(x) strrep(x,' ','  '),yTickLabel,'UniformOutput',false);
set(gca,'YTickLabel',yTickLabel);

strXlab.Position(2) = yTick(1)-0.05*(yTick(end)-yTick(1))-0.003;
% strYlab.Position(1) = ;
% strXlab.Position(1) = strXlab.Position(1)-10;

ha(12).YLabel.Position = [-0.3813,0.100000095367432,-1.000000000000014];
axis square
set(gca,'XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);

%% Save Fig in Different formats
set(ha(1:size(ha)),'box','off');
set(gcf,'color','white')
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV\vWM grid PLV','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV\vWM grid PLV','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV\vWM grid PLV','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\PLV\vWM grid PLV.fig');