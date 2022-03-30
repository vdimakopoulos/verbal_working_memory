%% Paths
Drive_Letter = 'F:\';
strPaths.Main = [Drive_Letter,'Vasileios\'];
strPaths.Project = [strPaths.Main, 'Task Analysis\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = [Drive_Letter,'Vasileios\Toolboxes\fieldtrip-20200315\'];

% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = [Drive_Letter,'Vasileios\Toolboxes\eeglab14_1_1b\'];

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(strPaths.Main)
addpath(genpath(strPaths.Project))
addpath(strPaths.Toolboxes.FieldTrip)
% Remove EEGLAB from path
rmpath(genpath(strPaths.Toolboxes.EEGLAB))

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

% open meeg
setenv('PATH','F:\Vasileios\Toolboxes\OpenMEEG\bin\')
setenv('LD_LIBRARY_PATH','F:\Vasileios\Toolboxes\OpenMEEG\lib\')
setenv('TMP','F:/Vasileios/Temp/')
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')


%% create figure
fig = figure;
set(gcf,'color','white')
set(fig,'Position',[680,681,591,297]);%[680,291,610,687]);680,461,591,517
% ha = tight_subplot(2,3,[.1 .08],[.12 .04],[.08, .05])
ha = tight_subplot(1,2,[.01 .01],[.01 .01],[.01, .01])


%% load data
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceStatistics_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceStatistics_group_maint.mat')
%% mean statistical source power
nSubjects = length(statint_group)
mean_stat_enc = statint_group{1};
sum_parcel_enc = statint_group{1}.parcel;
for i = 2:length(statint_group)
    sum_parcel_enc = sum_parcel_enc+statint_group{i}.parcel;
end
mean_parcel_enc = sum_parcel_enc./nSubjects;
mean_stat_enc.parcel = mean_parcel_enc;    

mean_stat_maint = statint_group_maint{1};
sum_parcel_maint = statint_group_maint{1}.parcel;
for i = 2:length(statint_group)
    sum_parcel_maint = sum_parcel_maint+statint_group_maint{i}.parcel;
end
mean_parcel_maint = sum_parcel_maint./nSubjects;
mean_stat_maint.parcel = mean_parcel_maint;    

%% visualization
axes(ha(1))
cfg = [];
cfg.figure = gcf;
cfg.axis = gca;
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 10];
cfg.funcolormap    = colormap(bluewhitered_neg(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, mean_stat_enc);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(ha(1),flipud(bluewhitered_neg(128)))
set(ha(1),'FontSize',12)
cb = get(ha(1),'colorbar');
set(cb,'Ticks',[0 10],'TickLabels',[0 10],'Position',[0.436628664946881,0.24820606060606,0.018532079553966,0.542944107744107])
text('Parent',ha(1),'Rotation',90,'String','t-value',...
    'Position',[-60.03510047318288,-151.3664016163944,-34.57126391474049],'FontSize',12);
set(ha(1),'Position',[-0.05000005,0.15,0.536666666666667,0.77])


axes(ha(2))
cfg = [];
cfg.figure = gcf;
cfg.axis = gca;
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 10];
cfg.funcolormap    = colormap(bluewhitered_pos(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, mean_stat_maint);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(ha(2),bluewhitered_pos(128))
colormap(ha(1),flipud(bluewhitered_neg(128)))
set(ha(2),'FontSize',12)
cb = get(ha(2),'colorbar');
set(cb,'Ticks',[0 10],'TickLabels',[0 10],'Position',[0.937474688635544,0.252525252525252,0.018532079553966,0.542944107744107])
text('Parent',ha(2),'Rotation',90,'String','t-value',...
    'Position',[-60.53385258936943,-148.9460746267742,-36.237637192915834],'FontSize',12);
set(ha(2),'Position',[0.4500005,0.15,0.536666666666667,0.77])

%% Textboxes

annotation(fig,'textbox',...
    [0.008460236886633,0.81718182089112,0.070219964570604,0.124579121869823],...
    'Color',[0 0 0],...
    'String',{'a'},...
    'FontSize',16,...
     'FontWeight','bold', ...
    'EdgeColor','none');


annotation(fig,'textbox',...
    [0.510998307952621,0.81718182089112,0.073604059224403,0.124579121869823],...
    'Color',[0 0 0],...
    'String',{'b'},...
    'FontSize',16,...
     'FontWeight','bold', ...
    'EdgeColor','none');


%% Save figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\Main Beamforming sources','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\Main Beamforming sources','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\Main Beamforming sources','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\Main Beamforming sources.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming sources_Fig\Main Beamforming sources.tiff');
