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
addpath(strPaths.GeneralFunctions)
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
fig = figure('units','normalized','outerposition',[0 0 1 1],'PaperSize',[20 11])
fig.Units = 'centimeters';
% fig = figure;
% set(fig,'PaperPositionMode','auto');
set(fig,'PaperOrientation','portrait');
% set(fig,'PaperUnits','normalized');
% 
% 
% 
% fig.Position =  [40 40 1200 800];



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
semilogx(Panel_a_Dat.freqAxis,abs(Panel_a_Dat.PLV_PairsMaint{nPair,iSS}),'k','LineWidth',3)
hold on
semilogx(Panel_a_Dat.freqAxis,abs(Panel_a_Dat.PLV_PairsFix{nPair,iSS}),'r','LineWidth',3)
hold on;
indMaxBandMaintFix = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_a_Dat.freqAxis,[0 1], 0.05,'k');
hold on;
semilogx(freqHi,abs(Panel_a_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_a_Dat.PLV_Pairs_FixHi{nPair,iSS}),'r','LineWidth',3);
ylim([0,1])
xlim([4 100])

set(ha(1),'Fontsize',13)
set(ha(1),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(1),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.207941176470588 0.936446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'iEEG Hipp - EEG P3'},...
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
set(gca,'FontSize',13)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar
% title('TFR iEEG Hipp - EEG P3','FontSize',13);

annotation('textbox',...
    [0.324033613445378 0.965565350540506 0.0259978985413909 0.0435663617188686],...
    'String',{'b'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
% set(gca, 'InvertHardcopy', 'off')

annotation('textbox',...
    [0.513130252100642 0.935433637798021 0.108193277311121 0.0314083073699607],...
    'Color',[1 1 1],...
    'String','iEEG Hipp - EEG P3',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

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
    [0.654033613445378 0.965565350540506 0.0259978985413909 0.0435663617188686],...
    'String',{'c'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);

% annotation('textbox',...
%     [0.64915966386535 0.931380953015048 0.0887605017129364 0.0314083073699607],...
%     'String',{'iEEG Hipp - EEG [6 7] Hz'},...
%     'FitBoxtoText','off',...
%     'FontSize',12,...
%     'EdgeColor','none');
% 

annotation('textbox',...
    [0.667016806722492 0.936446809023954 0.108193277311121 0.0314083073699607],...
    'String',{'iEEG Hipp - EEG [6 7] Hz'},...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Panel d
axes(ha(4));
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
        semilogx(Panel_d_Dat.freqAxis,abs(Panel_d_Dat.PLV_PairsMaint{(j+(s_loc-1)*nGridChans),... 
            iSS_to_plot}),'k','LineWidth',3);
        hold on;      
        semilogx(Panel_d_Dat.freqAxis,abs(Panel_d_Dat.PLV_PairsFix{(j+(s_loc-1)*nGridChans),...
            iSS_to_plot}),'r','LineWidth',3);
        semilogx(freqHi,abs(Panel_d_Dat.PLV_Pairs_MaintHi{nPair,iSS_to_plot}),'k','LineWidth',3);
        semilogx(freqHi,abs(Panel_d_Dat.PLV_Pairs_FixHi{nPair,iSS_to_plot}),'r','LineWidth',3);
        ylim([0 1])
        xlim([4 100]);
        hold on;
        Maint_Fix_Diff = abs(Panel_d_Dat.PLV_PairsMaint{(j+(s_loc-1)*nGridChans),iSS_to_plot}) ...
            - abs(Panel_d_Dat.PLV_PairsFix{(j+(s_loc-1)*nGridChans),iSS_to_plot});
        Vars = Panel_d_Dat.SignificanceBars.Scalp_Grid;
        indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_d_Dat.freqAxis,[0 1],0.05,'k');
    end
  
end

     
set(ha(4),'Fontsize',13)
set(ha(4),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(4),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.207941176470588 0.696446809023954 0.0861344537815125 0.0273556225869914],...
    'String',{'ECoG D6- EEG P3'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.71794 0.0259978985413909 0.0435663617188686],...
    'String',{'d'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

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

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);
colormap jet
set(ha(5),'FontSize',13)
% title('ECoG - EEG P3 [6 9] Hz','FontSize',13);


annotation('textbox',...
     [0.324033613445378 0.71794 0.0259978985413909 0.0435663617188686],...
    'String',{'e'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

annotation('textbox',...
    [0.498949579831733 0.695312057250918 0.108193277311122 0.0314083073699607],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','ECoG - EEG P3 [6 9] Hz',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
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
    [0.654033613445378 0.71794 0.0259978985413909 0.0435663617188686],...
    'String',{'f'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);

annotation('textbox',...
    [0.661764705882334 0.698351571049507 0.105567223901133 0.0314083073699607],...
    'String',{'ECoG D6 - EEG [6 9] Hz'},...
    'FontSize',12,...
    'EdgeColor','none');
%% Panel g
axes(ha(7))
Panel_g_Dat = FigureMaint_Data.PLV_Spectra_Grid_Hipp;
strPlotColors = {'r','g','b','k','c'};
iSS = 4;
freqHi =[30:5:100];

nChannelPairs = Panel_g_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,1))','_AHL2-_AHL3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,2))','_GL_H2')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

semilogx(Panel_g_Dat.freqAxis,abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS}),strPlotColors{4},'LineWidth',3);
hold on;
semilogx(Panel_g_Dat.freqAxis,abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS}),strPlotColors{1},'LineWidth',3);
hold on;
semilogx(freqHi,abs(Panel_g_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_g_Dat.PLV_Pairs_FixHi{nPair,iSS}),'r','LineWidth',3)
Vars = Panel_g_Dat.SignificanceBars.GridHipp;
Maint_Fix_Diff = abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS});
indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'k');
Fix_Maint_Diff =abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'r');


ylim([0,1]);
xlim([4 100])


set(ha(7),'Fontsize',13)
set(ha(7),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(7),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.207941176470588 0.4564 0.0861344537815125 0.0273556225869914],...
    'String',{'iEEG Hipp - ECoG H2'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.4718 0.0259978985413909 0.0435663617188686],...
    'String',{'g'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
    
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

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);
colormap jet
set(ha(8),'FontSize',13)
% title('iEEG Hipp - ECoG [9 10] Hz','FontSize',13);
annotation('textbox',...
    [0.324033613445378 0.4718 0.0259978985413909 0.0435663617188686],...
    'String',{'h'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

annotation('textbox',...
    [0.483193277310924 0.453164134251936 0.120273109243694 0.031408307369961],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','iEEG Hipp - ECoG [9 10] Hz',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');
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
    [0.654033613445378 0.4718 0.0259978985413909 0.0435663617188686],...
    'String',{'i'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);

annotation('textbox',...
    [0.662815126050401 0.45316413437272 0.109768904448182 0.0314083073699608],...
    'String',{'ECoG H2 - EEG [9 10] Hz'},...
    'FontSize',12,...
    'EdgeColor','none');


%% Panel j
axes(ha(10))
Panel_j_Dat = FigureMaint_Data.PLV_Spectra_Grid_Hipp;
nChannelPairs = Panel_j_Dat.nChannelPairs;
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,1))','_AHL2-_AHL3')));
ind2 = find(~cellfun(@isempty,strfind(Panel_g_Dat.dataBipolar.label(nChannelPairs(:,2))','_GL_C2')));
nPair = intersect(nPairs_ToRun,ind);
nPair = intersect(nPair,ind2);

strPlotColors = {'r','g','b','k','c'};
iSS = 4;
freqHi = [30:5:100];

semilogx(freqAxis,abs(Panel_j_Dat.PLV_PairsMaint{nPair,iSS}),strPlotColors{4},'LineWidth',3);
hold on;
semilogx(freqAxis,abs(Panel_j_Dat.PLV_PairsFix{nPair,iSS}),strPlotColors{1},'LineWidth',3);
hold on;
semilogx(freqHi,abs(Panel_j_Dat.PLV_Pairs_MaintHi{nPair,iSS}),'k','LineWidth',3);
semilogx(freqHi,abs(Panel_j_Dat.PLV_Pairs_FixHi{nPair,iSS}),'r','LineWidth',3);
Vars = Panel_j_Dat.SignificanceBars.GridHipp2;
Maint_Fix_Diff = abs(Panel_j_Dat.PLV_PairsMaint{nPair,iSS})-abs(Panel_j_Dat.PLV_PairsFix{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},Panel_g_Dat.freqAxis,[0 1],0.03,'r');indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'k');
Fix_Maint_Diff =abs(Panel_g_Dat.PLV_PairsFix{nPair,iSS})-abs(Panel_g_Dat.PLV_PairsMaint{nPair,iSS});
indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'r');

ylim([0,1]);
xlim([4 100])



set(ha(10),'Fontsize',13)
set(ha(10),'XTick',[10 100],'XTicklabel',[10 100]);
set(ha(10),'YTick',[0 0.5 1],'YTicklabel',[0 0.5 1]);

annotation('textbox',...
    [0.207941176470588 0.2164 0.0861344537815125 0.0273556225869914],...
    'String',{'iEEG Hipp - ECoG C2'},...
    'FontSize',14,...
    'FitBoxToText','on',...
    'EdgeColor','none');
     
annotation('textbox',...
    [0.0084033613445378 0.2418 0.0259978985413909 0.0435663617188686],...
    'String',{'j'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
    

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

colorbar('Ticks',[0 0.2 0.5 0.7],'FontSize',13);
colormap jet

annotation('textbox',...
    [0.324033613445378 0.2418 0.0259978985413909 0.0435663617188686],...
    'String',{'k'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
set(ha(11),'Fontsize',13)

% title('iEEG Hipp - ECoG [16 29] Hz','FontSize',13);

annotation('textbox',...
    [0.476365546218487 0.212029382478887 0.1234243697479 0.031408307369961],...
    'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
    'String','iEEG Hipp - ECoG [16 29] Hz',...
    'FontWeight','bold',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'EdgeColor','none');

%% Panel l 
axes(ha(12));
Panel_l_Dat = FigureMaint_Data.Granger;
Grngr_Enc = Panel_l_Dat.Grng_Encoding{1};
Grngr_Maint = Panel_l_Dat.Grng_Maint{1};
freq = Panel_l_Dat.freq;
dark_Green = [0.02 0.47 0.04];
gray       = [0.5 0.5 0.5];
Colors     = {dark_Green,'g','k',gray};


semilogx(freq,Grngr_Enc.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(freq,Grngr_Enc.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
semilogx(freq,Grngr_Maint.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
semilogx(freq,Grngr_Maint.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
xlim([4 30]);
ylim([0 0.2]);

lgd = legend('ECoG C2- iEEG Hipp enc','iEEG Hipp - ECoG C2 enc','iEEG Hipp - ECoG C2 maint','ECoG C2 - iEEG Hipp maint');
set(ha(12),'FontSize',13);
set(ha(12),'YTick',[0 0.1 0.2],'YTickLabel', [0 0.1 0.2]);
set(ha(12),'XTick',[0 10 20 30],'XTickLabel', [0 10 20 30]);

lgd.Location = 'northwest';
lgd.FontSize =10;
annotation('textbox',...
    [0.654033613445378 0.2418 0.0259978985413909 0.0435663617188686],...
    'String',{'l'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'EdgeColor','none');
set(ha(12),'Fontsize',13)


%% Save Fig in Different formats
set(ha(1:size(ha)),'box','off');
set(gcf,'color','white')
set(gcf, 'InvertHardcopy', 'off')

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PLV','-dpdf');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PLV','-dpng');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PLV','-depsc');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\vWM grid PLV.fig');