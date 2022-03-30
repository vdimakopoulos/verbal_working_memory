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
strPaths.FigureVars = [strPaths.FigureData,'210122 Granger_spectra\Delta_Granger_Incorrect_Trials_EEG'];
load(strPaths.FigureVars)


% load ECoG Delta Granger incorrect
load('F:\Vasileios\Task Analysis\Analysis Results\Incorrect Trials\ECoG Granger\DeltaGranger_IncorrectTrials_ECoG.mat')
%%
x_enc = 5;
x_maint =46;
nSubjs = length(Dgranger_IncorrectTrials_EEG);
for i = 1:nSubjs
   EncDiff(i) =  Dgranger_IncorrectTrials_EEG{i}.Enc;
   MaintDiff(i) = Dgranger_IncorrectTrials_EEG{i}.Maint;
end
% MaintDiff(13) = MaintDiff_Reverse(7)
[EncDiff_sorted,Ienc] =  sort(EncDiff,'descend');
[MaintDiff_sorted,Imaint] = sort(MaintDiff);


xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

% xValues = [1:15; repmat(46,15,1)']'%[1:15]+3*15]';
fig = figure;
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
MaintDiff_sorted = sort([MaintDiff_sorted,MaintDiff_ECoG])

median_Maint = median(MaintDiff_sorted)*100;
prcntile75 = prctile(MaintDiff_sorted,75);
prcntile25 = prctile(MaintDiff_sorted,25);
interquartile = iqr(MaintDiff_sorted);
Upper_adjacent = prcntile75 -1.5*interquartile;
Lower_adjacent = prcntile25 +1.5*interquartile;

line(repmat(x_maint+2,2,1)',[prcntile75 prcntile25]*100,'Color','k')
plot([47:49],repmat(median_Maint,3,1),'r','LineWidth',1)
plot([47:49],repmat(median_Maint,3,1)-0.4,'r','LineWidth',1)
plot([47:49],repmat(prcntile75,3,1) *100,'k');
plot([47:49],repmat(prcntile25,3,1) *100,'k');

%% Median for Encoding
EncDiff_sorted = sort([EncDiff_sorted,EncDiff_ECoG_sorted]);

median_Enc   = median(EncDiff_sorted)*100;
prcntile75 = prctile(EncDiff_sorted,75);
prcntile25 = prctile(EncDiff_sorted,25);
interquartile = iqr(EncDiff_sorted);
Upper_adjacent = prcntile75 -1.5*interquartile;
Lower_adjacent = prcntile25 +1.5*interquartile;

line(repmat(x_enc-2,2,1)',[prcntile75 prcntile25]*100,'Color','k')
plot([2:4],repmat(median_Enc,3,1),'r','LineWidth',1)
plot([2:4],repmat(median_Enc,3,1)-0.4,'r','LineWidth',1)
plot([2:4],repmat(prcntile75,3,1) *100,'k');
plot([2:4],repmat(prcntile25,3,1) *100,'k');

%% Axes Properties
set(gca,'FontSize',16,'box','off');
xlabel('Participants');
ylabel('\DeltaGranger (%)');
% l= legend([enc,maint],'encoding','maintenance','EdgeColor','none','Orientation','horizontal','Position',...
%     [0.137202388881927,0.841865082026002,0.517857132852078,0.069047617273671])
ylim([-20 20])
set(gca,'XTick',[5 46],'XTickLabel',[{'Encoding','Maintenance'}],'TickDir','out');


%% Statistics
% 
% n1 = length(find(MaintDiff_sorted>0));
% N1 = nSubjs;
% 
% n2 = length(find(DGrangerMaint_Correct>0));
% N2 = nSubjs;
% 
% x1 = [repmat('a',N1,1); repmat('b',N2,1)];
% x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
% [tbl,chi2stat,pval] = crosstab(x1,x2)
% 
% 
% 
% n1 = length(find(EncDiff_sorted<0));
% N1 = nSubjs;
% 
% n2 = length(find(DGrangerEnc_Correct<0));
% N2 = nSubjs;
% 
% x1 = [repmat('a',N1,1); repmat('b',N2,1)];
% x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
% [tbl,chi2stat,pval] = crosstab(x1,x2)

%% Save Figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off');

print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\DeltaGrangerIncorrectTrials','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\DeltaGrangerIncorrectTrials','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\DeltaGrangerIncorrectTrials','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\DeltaGrangerIncorrectTrials.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\BehaviourGranger\DeltaGrangerIncorrectTrials.tif');

