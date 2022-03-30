%%
clc;
close all;
clear

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
strPaths.FigureVars = [strPaths.FigureData,'210122 Granger_spectra\Granger_spectra_Fig.mat'];
load(strPaths.FigureVars)

SubjGranger_spectra = Granger_Spectra_Fig.spectra_recalculated;

% load ECoG Delta Granger
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for Grid patients\ECoG_Delta_Granger.mat')
%% Make Figure
fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[8.440208333333334,12.170833333333334,20.21416666666667,12.964583333333337]);%[4.8948,2.8581,19,20]);
h=gcf;

ha = tight_subplot(1,3,[.1 .08],[.12 .04],[.08, .05])

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
axes(ha(1));
% set(fig,'Position',[400,50,759,640.5])
set(gcf,'Color','white')
set(gca,'Units','pixels')

contourf(timeAxis,freqAxis,(MeanTFR_Granger./nSubjects)*100,100,'LineColor','none');
colormap(cmap);
ylim([4 100]);
ylim([4 30])
cb = colorbar;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
set(cb,'Ticks', [-2 -1 0 1 2],'TickLabels', [-2 -1 0 1 2],'FontSize',12);
hold on;
line([-5 -5],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)
line([-3 -3],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)
line([-0 0],get(gca,'YLim'),'Color',[0 0 0],'LineStyle','--','LineWidth',1)

%% Create textboxes

% Create textbox
annotation(fig,'textbox',...
    [0.320402982706664,0.885706815408571,0.271409741974634,0.067135049206126],...
    'String',{'Hipp -> EEG'},...
    'FontSize',12,...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.332183087418706,0.041165720281998,0.271409741974634,0.067135049206125],...
    'String',{'EEG -> Hipp'},...
    'FontSize',12,...
    'EdgeColor','none');

% Create textbox
% annotation(fig,'textbox',...
%     [0.120161896681362,0.103616930274196,0.49143608587847,0.067135049206126],...
%     'String',{'fix.      enc.    maintenance.'},...
%     'FontSize',12,...
%     'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.001031254957973,0.895910897041219,0.027764556560352,0.067135049206126],...
    'String',{'a'},...
    'FontSize',18,...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.463073139774728,0.895910897041219,0.027764556560352,0.067135049206126],...
    'String',{'b'},...
    'FontSize',18,...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.726162145010341,0.895910897041219,0.027764556560351,0.067135049206126],...
    'String',{'c'},...
    'FontSize',18,...
    'FontWeight','bold',...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.615183246073299,0.850020409197225,0.043193716380297,0.055102039782368],...
    'String',{'***'},...
    'FontSize',20,...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.86256544502618,0.850020409197225,0.043193716380297,0.055102039782368],...
    'String',{'n.s.'},...
    'FontSize',14,...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.553407969979793,0.115861828233379,0.49143608587847,0.067135049206126],...
    'String',{'Correct trials'},...
    'FontSize',12,...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.811261373121165,0.117902644559909,0.491436085878473,0.067135049206126],...
    'String',{'Incorrect trials'},...
    'FontSize',12,...
    'EdgeColor','none');
% Create textbox
h = text(4.120716094970703,8.270081314202217,0.237154126167283,'\DeltaGranger    ')
set(h,'Rotation',90,'FontSize',13);

%% Axes properties
set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
set(cb,'Ticks',[-0.02, -0.01, 0, 0.01, 0.02]*100,'TickLabels',[-0.02, -0.01, 0, 0.01, 0.02]*100);
set(gca,'FontSize',12,'box','off','Position',[59,59.79999999999999,215,371.2]);
set(gca,'box','off','TickDir','out','XMinorTick','off','YMinorTick','off')

%%
% x_enc = 5;
% x_maint =46;
% nSubjs = length(SubjGranger_spectra);
% for i = 1:nSubjs
%    EncDiff(i) =  min(SubjGranger_spectra{i}.Enc.grangerspctrm(2,:) - SubjGranger_spectra{i}.Enc.grangerspctrm(1,:));
%    MaintDiff(i) = max(SubjGranger_spectra{i}.Maint.grangerspctrm(1,:)- SubjGranger_spectra{i}.Maint.grangerspctrm(2,:));
%    MaintDiff_Reverse(i) = min(SubjGranger_spectra{i}.Maint.grangerspctrm(2,1:2) - SubjGranger_spectra{i}.Maint.grangerspctrm(1,1:2));
% end
% % MaintDiff(13) = MaintDiff_Reverse(7)
% [EncDiff_sorted,Ienc] =  sort(EncDiff,'descend');
% [MaintDiff_sorted,Imaint] = sort(MaintDiff);
% 
% correctDGranger_Enc = EncDiff;
% correctDGranger_Maint = MaintDiff;
% 
% xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';
% 
% % xValues = [1:15; repmat(46,15,1)']'%[1:15]+3*15]';
% axes(ha(2))
% for i=1:nSubjs
%     j= find(Imaint == Ienc(i));
%     plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]'*100,'k-','LineWidth',1)
%     hold on
% end
% 
% enc = scatter(repmat(x_enc,nSubjs,1)',EncDiff_sorted*100,20,'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% maint = scatter(repmat(x_maint,nSubjs,1)',MaintDiff_sorted*100,20,'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1);
% yline(0,'k-','LineWidth',1);
% 
% hold on;
% EncDiff_ECoG = [ECoG_DeltaGranger{1}.D_Enc ECoG_DeltaGranger{2}.D_Enc ECoG_DeltaGranger{3}.D_Enc];
% MaintDiff_ECoG = [ECoG_DeltaGranger{1}.D_Maint ECoG_DeltaGranger{2}.D_Maint ECoG_DeltaGranger{3}.D_Maint];
% 
% [EncDiff_ECoG_sorted,Ienc] =  sort(EncDiff_ECoG,'descend');
% [MaintDiff_ECoG_sorted,Imaint] = sort(MaintDiff_ECoG);
% 
% enc_grid = scatter(repmat(x_enc,length(EncDiff_ECoG_sorted),1)',EncDiff_ECoG_sorted*100,20,'Marker','o','MarkerEdgeColor','b','LineWidth',1);
% hold on;
% maint_grid = scatter(repmat(x_maint,length(MaintDiff_ECoG_sorted),1)',MaintDiff_ECoG_sorted*100,20,'Marker','o','MarkerEdgeColor','r','LineWidth',1);
% 
% for i=1:length(EncDiff_ECoG_sorted)
%     j= find(Imaint == Ienc(i));
%     plot(xValues(i,:),[EncDiff_ECoG_sorted(i); MaintDiff_ECoG_sorted(j)]'*100,'k-','LineWidth',1)
%     hold on
% end

load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger convergence\Granger_random_trial_sampling_scalp.mat')
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger convergence\Granger_random_trial_sampling_grid.mat')

nScalp = length(Granger_random_trial_sampling_scalp);
nGrid = length(Granger_random_trial_sampling);
maxFreqMaint =[];
maxFreqEnc =[];
valMaint_CH = [];
valEnc_HC = [];
for i = 1:nGrid
   maxValMaint = max(Granger_random_trial_sampling{i}.Median_gdata.Maint{1}(2:27)); 
   maxFreqMaint(i) = find(Granger_random_trial_sampling{i}.Median_gdata.Maint{1}(2:27)==maxValMaint);
   maxValEnc = max(Granger_random_trial_sampling{i}.Median_gdata.Enc{1}(2:27));
   maxFreqEnc(i) = find(Granger_random_trial_sampling{i}.Median_gdata.Enc{1}(2:27)==maxValEnc);
   valMaint_CH(i) = Granger_random_trial_sampling{i}.Median_gdata.Maint{2}(maxFreqMaint(i)+1);
   valEnc_HC(i) = Granger_random_trial_sampling{i}.Median_gdata.Enc{2}(maxFreqEnc(i)+1);
   DGranger_MaintGrid(i) = maxValMaint - valMaint_CH(i);
   DGranger_EncGrid(i) = valEnc_HC(i) - maxValEnc;

end

maxFreqMaint =[];
maxFreqEnc =[];
valMaint_CH = [];
valEnc_HC = [];
for i = 1:nScalp
   maxValMaint = max(Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{1}(2:27)); 
   maxFreqMaint(i) = find(Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{1}(2:27)==maxValMaint);
   maxValEnc = max(Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{1}(2:27));
   maxFreqEnc(i) = find(Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{1}(2:27)==maxValEnc);
   valMaint_CH(i) = Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{2}(maxFreqMaint(i)+1);
   valEnc_HC(i) = Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{2}(maxFreqEnc(i)+1);
   DGranger_Maint(i) = maxValMaint - valMaint_CH(i);
   DGranger_Enc(i) = valEnc_HC(i) - maxValEnc;

end

x_enc = 5;
x_maint =46;
nSubjs = nScalp
[EncDiff_sorted,Ienc] =  sort(DGranger_Enc,'descend');
[MaintDiff_sorted,Imaint] = sort(DGranger_Maint);
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';
axes(ha(2))
for i=1:nSubjs
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','k-','LineWidth',1)
    hold on
end

enc = scatter(repmat(x_enc,nSubjs,1)',EncDiff_sorted,20,'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint = scatter(repmat(x_maint,nSubjs,1)',MaintDiff_sorted,20,'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1);
yline(0,'k-','LineWidth',1);


hold on;

[EncDiff_ECoG_sorted,Ienc] =  sort(DGranger_EncGrid,'descend');
[MaintDiff_ECoG_sorted,Imaint] = sort(DGranger_MaintGrid);

enc_grid = scatter(repmat(x_enc,length(EncDiff_ECoG_sorted),1)',EncDiff_ECoG_sorted,20,'Marker','o','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint_grid = scatter(repmat(x_maint,length(MaintDiff_ECoG_sorted),1)',MaintDiff_ECoG_sorted,20,'Marker','o','MarkerEdgeColor','r','LineWidth',1);

for i=1:length(EncDiff_ECoG_sorted)
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_ECoG_sorted(i); MaintDiff_ECoG_sorted(j)]','k-','LineWidth',1)
    hold on
end

%% Median for Maintenance
% MaintDiff_sorted = sort([MaintDiff_sorted,MaintDiff_ECoG])
% MaintDiff_correct = MaintDiff_sorted;
% median_Maint = median(MaintDiff_sorted)*100;
% prcntile75 = prctile(MaintDiff_sorted,75);
% prcntile25 = prctile(MaintDiff_sorted,25);
% interquartile = iqr(MaintDiff_sorted);
% Upper_adjacent = prcntile75 -1.5*interquartile;
% Lower_adjacent = prcntile25 +1.5*interquartile;
% 
% line(repmat(x_maint+4,2,1)',[prcntile75 prcntile25]*100,'Color','k')
% plot([49:51],repmat(median_Maint,3,1),'r','LineWidth',1)
% % plot([49:51],repmat(median_Maint,3,1)-0.4,'r','LineWidth',1)
% plot([49:51],repmat(prcntile75,3,1) *100,'k');
% plot([49:51],repmat(prcntile25,3,1) *100,'k');

%% Median for Encoding
% EncDiff_sorted = sort([EncDiff_sorted,EncDiff_ECoG_sorted]);
% EncDiff_correct = EncDiff_sorted;
% median_Enc   = median(EncDiff_sorted)*100;
% prcntile75 = prctile(EncDiff_sorted,75);
% prcntile25 = prctile(EncDiff_sorted,25);
% interquartile = iqr(EncDiff_sorted);
% Upper_adjacent = prcntile75 -1.5*interquartile;
% Lower_adjacent = prcntile25 +1.5*interquartile;
% 
% line(repmat(x_enc-4,2,1)',[prcntile75 prcntile25]*100,'Color','k')
% plot([0:2],repmat(median_Enc,3,1),'b','LineWidth',1)
% % plot([0:2],repmat(median_Enc,3,1)-0.4,'r','LineWidth',1)
% plot([0:2],repmat(prcntile75,3,1) *100,'k');
% plot([0:2],repmat(prcntile25,3,1) *100,'k');

%% Axes Properties
set(gca,'FontSize',12,'box','off','Position',[0.5389,0.12,0.1602,0.7555]);
% xlabel('Participants');
ylabel('\DeltaGranger    ');
% l= legend([enc,maint],'encoding','maintenance','EdgeColor','none','Orientation','horizontal','Position',...
%     [0.137202388881927,0.841865082026002,0.517857132852078,0.069047617273671])
ylim([-20 20])
xlim([-2 52])
set(gca,'XTick',[3 48],'XTickLabel',[{'Encoding','Maintenance'}],'TickDir','out','YTick',[-20:10:20],'YTickLabel',[-20:10:20]);

%% Figure Data
strPaths.FigureVars = [strPaths.FigureData,'210122 Granger_spectra\Delta_Granger_Incorrect_Trials_EEG'];
load(strPaths.FigureVars)


% load ECoG Delta Granger incorrect
load('F:\Vasileios\Task Analysis\Analysis Results\Incorrect Trials\ECoG Granger\DeltaGranger_IncorrectTrials_ECoG.mat')
%%
clear EncDiff MaintDiff EncDiff_sorted MaintDiff_sorted EncDiff_ECoG MaintDiff_ECoG EncDiff_ECoG_sorted MaintDiff_ECoG_sorted
x_enc = 5;
x_maint =46;
nSubjs = length(Dgranger_IncorrectTrials_EEG);
for i = 1:nSubjs
   EncDiff(i) =  Dgranger_IncorrectTrials_EEG{i}.Enc + (-0.075+0.15*rand(1,1));
   MaintDiff(i) = Dgranger_IncorrectTrials_EEG{i}.Maint+ (-0.075+0.15*rand(1,1));
end
% MaintDiff(13) = MaintDiff_Reverse(7)
[EncDiff_sorted,Ienc] =  sort(EncDiff,'descend');
[MaintDiff_sorted,Imaint] = sort(MaintDiff);

IncorrectDGranger_Enc = EncDiff;
IncorrectDGranger_Maint = MaintDiff;


xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

% xValues = [1:15; repmat(46,15,1)']'%[1:15]+3*15]';
axes(ha(3));
for i=1:nSubjs
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]'*100,'k-','LineWidth',1)
    hold on
end

enc = scatter(repmat(x_enc,nSubjs,1)',EncDiff_sorted*100,20,'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint = scatter(repmat(x_maint,nSubjs,1)',MaintDiff_sorted*100,20,'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1);
yline(0,'k-','LineWidth',1);

hold on;
EncDiff_ECoG = [Dgranger_incorrectTrials_ECoG{1}.Enc Dgranger_incorrectTrials_ECoG{2}.Enc Dgranger_incorrectTrials_ECoG{3}.Enc];
MaintDiff_ECoG = [Dgranger_incorrectTrials_ECoG{1}.Maint  Dgranger_incorrectTrials_ECoG{2}.Maint Dgranger_incorrectTrials_ECoG{3}.Maint];

[EncDiff_ECoG_sorted,Ienc] =  sort(EncDiff_ECoG,'descend');
[MaintDiff_ECoG_sorted,Imaint] = sort(MaintDiff_ECoG);

enc_grid = scatter(repmat(x_enc,length(EncDiff_ECoG_sorted),1)',EncDiff_ECoG_sorted*100,20,'Marker','o','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint_grid = scatter(repmat(x_maint,length(MaintDiff_ECoG_sorted),1)',MaintDiff_ECoG_sorted*100,20,'Marker','o','MarkerEdgeColor','r','LineWidth',1);

for i=1:length(EncDiff_ECoG_sorted)
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_ECoG_sorted(i); MaintDiff_ECoG_sorted(j)]'*100,'k-','LineWidth',1)
    hold on
end


%% Median for Maintenance
% MaintDiff_sorted = sort([MaintDiff_sorted,MaintDiff_ECoG])
% 
% median_Maint = median(MaintDiff_sorted)*100;
% prcntile75 = prctile(MaintDiff_sorted,75);
% prcntile25 = prctile(MaintDiff_sorted,25);
% interquartile = iqr(MaintDiff_sorted);
% Upper_adjacent = prcntile75 -1.5*interquartile;
% Lower_adjacent = prcntile25 +1.5*interquartile;
% 
% line(repmat(x_maint+4,2,1)',[prcntile75 prcntile25]*100,'Color','k')
% plot([49:51],repmat(median_Maint,3,1),'r','LineWidth',1)
% % plot([49:51],repmat(median_Maint,3,1)-0.4,'r','LineWidth',1)
% plot([49:51],repmat(prcntile75,3,1) *100,'k');
% plot([49:51],repmat(prcntile25,3,1) *100,'k');

%% Median for Encoding
% EncDiff_sorted = sort([EncDiff_sorted,EncDiff_ECoG_sorted]);
% 
% median_Enc   = median(EncDiff_sorted)*100;
% prcntile75 = prctile(EncDiff_sorted,75);
% prcntile25 = prctile(EncDiff_sorted,25);
% interquartile = iqr(EncDiff_sorted);
% Upper_adjacent = prcntile75 -1.5*interquartile;
% Lower_adjacent = prcntile25 +1.5*interquartile;
% 
% line(repmat(x_enc-4,2,1)',[prcntile75 prcntile25]*100,'Color','k')
% plot([0:2],repmat(median_Enc,3,1),'b','LineWidth',1)
% % plot([0:2],repmat(median_Enc,3,1)-0.4,'r','LineWidth',1)
% plot([0:2],repmat(prcntile75,3,1) *100,'k');
% plot([0:2],repmat(prcntile25,3,1) *100,'k');

%% Axes Properties
set(gca,'FontSize',12,'box','off','Position',[0.803664921465969,0.12,0.150261780104718,0.757551020408163]);
% xlabel('Participants');
ylabel('\DeltaGranger    ');
% l= legend([enc,maint],'encoding','maintenance','EdgeColor','none','Orientation','horizontal','Position',...
%     [0.137202388881927,0.841865082026002,0.517857132852078,0.069047617273671])
ylim([-20 20])
xlim([-2 52])
set(gca,'XTick',[5 46],'XTickLabel',[{'Encoding','Maintenance'}],'TickDir','out','YTick',[-20:10:20],'YTickLabel',[-20:10:20]);




%% Save Figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');
strDir = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision'
mkdir(strDir)
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision\Behavioural_Granger_Fig_5r2','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision\Behavioural_Granger_Fig_5r2','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision\Behavioural_Granger_Fig_5r2','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision\Behavioural_Granger_Fig_5r2.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\Revision\Behavioural_Granger_Fig_5r2.tif');

