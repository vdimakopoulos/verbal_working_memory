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
addpath(genpath('F:\Vasileios\Toolboxes\iELvis-iElvis_pm\'))
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

strPaths.FigureVars2 = [strPaths.FigureData,'FigureMaintData.mat'];
load(strPaths.FigureVars2)

strPaths.FigureVars3 = [strPaths.FigureData,'FigurePSD_Data.mat'];
load(strPaths.FigureVars3)

ElectrodeImagesPath = [strPaths.FigureData,'Main Granger Figure\'];
strPaths.ImageP42 = [ElectrodeImagesPath,'P42_electrodes.jpg'];

%% Read Image
% I1 = imread(strPaths.ImageP42);


%% Make Figure
fig = figure;
set(fig,'Units', 'centimeters')
set(fig,'Position',[8.440208333333334,6.958541666666668,17.885833333333338,18.176875000000006]);%[4.8948,2.8581,19,20]);
h=gcf;

ha = tight_subplot(4,3,[.1 .08],[.12 .04],[.08, .05])
% ha = tight_subplot(2,5,[.01 .08],[.1 .04],[.03 .03])%[0.05 0.05],[.1 .05],[.1 .1])%[0.05 0.03],[0.05 0.05])
set(ha(1:size(ha,1)),'Fontsize',20)


%% Plot the electrode locations

myPatient='USZ';
globalFsDir = 'F:\Vasileios\Task Analysis\Data\Freesurfer Data\'
filename=fullfile(globalFsDir,myPatient,'mri\orig.mgz');

% load MRI
mri=ft_read_mri(filename);
mri.coordsys='ras';

% segment scalp
cfg_seg=[];
cfg_seg.output='scalp';
scalp=ft_volumesegment(cfg_seg,mri);

% prepare mesh
cfg_mesh=[];
cfg_mesh.numvertices=65536;
cfg_mesh.tissue={'scalp'};
scalp_mesh=ft_prepare_mesh(cfg_mesh,scalp);


% copy coordinates after inspecting mesh
switch myPatient
    case 'USZ'%'USZpostopMR'
        nas=[7.535 81.77 -0.8851]; rpa=[77.53 6.774 -47.89]; lpa=[-82.47 17.77 -35.89];
    case 'USZ4'
        nas=[2.658 106.8 26.02]; rpa=[77.66 19.78 -19.98]; lpa=[-70.34 19.78 -24.98];
    case 'USZ5'
        nas=[2.016 95.26 -8.905]; rpa=[77.02 12.26 -28.9]; lpa=[-71.98 9.26 -27.9];
end
tenten=ft_read_sens(fullfile(strPaths.Toolboxes.FieldTrip,'template\electrode\standard_1020.elc'));

% align electrodes to fiducials
cfg_align2fid=[];
cfg_align2fid.method='fiducial';
cfg_align2fid.target.elecpos(1,:)=nas;
cfg_align2fid.target.elecpos(2,:)=lpa;
cfg_align2fid.target.elecpos(3,:)=rpa;
cfg_align2fid.target.chanpos=cfg_align2fid.target.elecpos;
cfg_align2fid.target.label={'Nz','LPA','RPA'};
cfg_align2fid.elec=tenten;
cfg_align2fid.fiducial={'Nz','LPA','RPA'};
tenten_fid=ft_electroderealign(cfg_align2fid);
cfg_snap=[];
cfg_snap.method='project';
cfg_snap.elec=tenten_fid;
cfg_snap.headshape=scalp_mesh;
tenten_snap=ft_electroderealign(cfg_snap);
elec_2keep={'T3'}; % for USZpostopMR

elec_2keep_mask=ismember(tenten_snap.label,elec_2keep);
elec=tenten_snap;
elec.chanpos=elec.chanpos(elec_2keep_mask,:);
elec.chantype=elec.chantype(elec_2keep_mask);
elec.chanunit=elec.chanunit(elec_2keep_mask);
elec.elecpos=elec.elecpos(elec_2keep_mask,:);
elec.label=elec.label(elec_2keep_mask);

axes(ha(1));
clear elecNames elecColors
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColors((i-1)*8+j) = -300;
        if (i-1)*8+j == 18 %% contact C2
            elecColors((i-1)*8+j) = 60;
        end
        
        if (i-1)*8+j == 59  %% contact H3
            elecColors((i-1)*8+j) = -200;
        end
        
    end
end

cfg = [];
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.elecColors = elecColors';
cfg.elecNames = elecNames'
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg. elecCbar = 'n';
cfg.title = [];
cfg.axis = gca;
cfgOut = plotPialSurf('USZ',cfg)

set(gca,... ok
    'Position',[0.03,0.825,0.236666666666667,0.135]);
clear Axes_objs
Axes_objs = get(gca,'Children');
for i =1:length(Axes_objs)-2
    set(Axes_objs(i),'FaceColor',[0 0 0])
    if i == 47 %% Contact C2
        set(Axes_objs(i),'FaceColor',[0 1 1])
    elseif i== 6 %% Contact H3
        set(Axes_objs(i),'FaceColor',[1 0 1])
        
    end
end
hold on;
plot3(elec.elecpos(:,1),elec.elecpos(:,2),elec.elecpos(:,3),'Marker','o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',5,'Color','none');
hold off;
title('a \color{cyan}C2, \color{magenta}H3')

%% TFR PSD
clim = FigurePSD_Data.TFR_PSD.ECoG.clim;
freq = FigurePSD_Data.TFR_PSD.ECoG.freq;
time = FigurePSD_Data.TFR_PSD.ECoG.time;
TFR_PSD_ECoG = FigurePSD_Data.TFR_PSD.ECoG.TFR_baselinedPSD;

axes(ha(2));

contourf(time,freq,TFR_PSD_ECoG/4,100,'LineColor','none');
set(gca,'clim',[0 1],'yscale','log');
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]); %background color of grid;
set(gca,'FontSize',10);
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
colormap(gca,jet);
colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4])%,'Location')%,'northoutside');
line([-5 -5],get(ha(2),'YLim'),'Color',[1 1 1]);
line([-3 -3],get(ha(2),'YLim'),'Color',[1 1 1]);
line([0 0],get(ha(2),'YLim'),'Color',[1 1 1]);
ylabel('Frequency (Hz)')%,'VerticalAlignment','middle');
strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap');
title('b H3')
% set(strXlab,'Position',[-1.999996185302734,3.087865666894011,1])

% set(gca, 'Position',...
%     [0.250000000000007,0.657142857142854,0.119109947643973,0.21918367346938]);

%% TFR PSD parietal ECoG
for i = 1:2
    fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\210219_Supplementary Figure S42\Supplementary_Figure_S42_Parietal_ECoG_TFR_PSD.fig';
    H_sup1 = openfig(fig_path,'invisible');
    ax_sup1 = gca;
    figS1a = get(ax_sup1,'children');
    axes(ha(3));
    
    copyobj(figS1a,gca)
    colormap(gca,jet);
    
    set(gca,'clim',[0 1],'yscale','log');
    set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]); %background color of grid;
    set(gca,'FontSize',10);
    set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
    colorbar('Ticks',[0 1],'TickLabels',[0 1])%,'Location')%,'northoutside');
    line([-5 -5],get(ha(3),'YLim'),'Color',[1 1 1]);
    line([-3 -3],get(ha(3),'YLim'),'Color',[1 1 1]);
    line([0 0],get(ha(3),'YLim'),'Color',[1 1 1]);
    ylabel('Frequency (Hz)')%,'VerticalAlignment','middle');
    strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap');
    % set(gca,'Position',...
    %     [0.476171204188482,0.657142857142854,0.119109947643973,0.21918367346938]);
    % set(strXlab,'Position',[-1.999996185302734,3.087865666894011,1])
end
title('c C2')

%% PSD Topo Grid Enc
axes(ha(4));
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end

strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\PSD_PLV_Granger_topos.mat';
load(strVar_toLoad)
elecColors = reshape(PSD_PLV_Granger.PSD_topo_enc,1,64)
cfg = [];
cfg.elecColorScale =[2 3]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
% 
% for i =1:2
%     
%     fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\42 DS Results\S42_PSD_Topo_Grid_enc_60_80Hz.fig';
%     H_sup1 = openfig(fig_path,'invisible');
%     ax_sup1 = gca;
%     figS1a = get(ax_sup1,'children');
%     axes(ha(4));
%     xlim([1 8])
%     ylim([1 8])
%     copyobj(figS1a,ha(4));
%     colormap(gca,'bluewhitered(128)');
%     set(gca,'clim',[0 1]);
%     colorbar('Ticks',[0 1],'FontSize',10)%,'Location')%,'northoutside');
%     % set(gca,'Position',...
%     %     [0.670000000000008,0.246938775510192,0.1191,0.2192]);
% end
% ax1 = gca;
% ax1.XTick = 1:8;
% ax1.YTick = 1:8;
% ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
% ax1.YTickLabel = 1:8;
% ax1.YDir = 'normal'
% xlabel('Grid column')
% ylabel('Grid row')
% set(gca,'dataAspectRatio',[1 1 1])
% title('d [-3.5 -3] s, [60 80] Hz')

%% PSD Hipp

clim = FigurePSD_Data.TFR_PSD.iEEG_Hipp.clim
freq = FigurePSD_Data.TFR_PSD.iEEG_Hipp.freq;
time = FigurePSD_Data.TFR_PSD.iEEG_Hipp.time
TFR_PSD_Hipp = FigurePSD_Data.TFR_PSD.iEEG_Hipp.TFRbaselinedPSD;

axes(ha(5))

contourf(time,freq,(TFR_PSD_Hipp-3.5)./3.5,100,'LineColor','none')
set(gca,'clim',[-1 0],'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca,jet);
colorbar('Ticks',[-1:1:4],'TickLabels',[-1:1:4])%,'Location')%,'northoutside');
ylabel('Frequency (Hz)')%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap')
% set(strXlab,'Position',[-1.999996185302734,3.199670443946102,1])
line([-5 -5],get(ha(5),'YLim'),'Color',[1 1 1])
line([-3 -3],get(ha(5),'YLim'),'Color',[1 1 1])
line([0 0],get(ha(5),'YLim'),'Color',[1 1 1])
title('e Hipp')
% set(gca, 'Position',...
%     [0.670000000000008,0.657142857142854,0.119109947643973,0.21918367346938]);

%% PSD Topo Grid Maint

axes(ha(6));
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end

strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
load(strVar_toLoad)
PSD_maint_topo = reshape(PLV_PSD_Granger_Topologies.PSD_maint,1,64);
ind = find(PSD_maint_topo>(max(PSD_maint_topo) - max(PSD_maint_topo)*0.2))
PSD_maint_topo(setdiff([1:length(PSD_maint_topo)],ind)) = 0;
elecColors = PSD_maint_topo;
cfg = [];
cfg.elecColorScale =[0 1]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
% 
% for i =1:2
%     
%     fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\42 DS Results\S42_PSD_Topo_Grid_maint.fig';
%     H_sup1 = openfig(fig_path,'invisible');
%     ax_sup1 = gca;
%     figS1a = get(ax_sup1,'children');
%     axes(ha(6));
%     xlim([1 8])
%     ylim([1 8])
%     copyobj(figS1a,ha(6));
%     colormap(gca,'bluewhitered(128)');
%     set(gca,'clim',[0 1]);
%     colorbar('Ticks',[0 1],'FontSize',10)%,'Location')%,'northoutside');
%     % set(gca,'Position',...
%     %     [0.670000000000008,0.246938775510192,0.1191,0.2192]);
% end
% ax1 = gca;
% ax1.XTick = 1:8;
% ax1.YTick = 1:8;
% ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
% ax1.YTickLabel = 1:8;
% ax1.YDir = 'normal'
% xlabel('Grid column')
% ylabel('Grid row')
% title('f [-2 0] s, [11 14] Hz')


%% PLV Hipp ECoG
for i =1:2
fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\210219_Supplementary Figure S42\PLV_Hipp_ECoG_C2.fig';
H_sup1 = openfig(fig_path,'invisible');
ax_sup1 = gca;
figS1a = get(ax_sup1,'children');
axes(ha(7));

copyobj(figS1a,gca)
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off',...
    'TickLength',[0.02 0.025],'YTick',[0 1],'YTickLabel',[0 1],'XTick',[4 10 20 100],'XTickLabel',[4 10 20 100]);

end
strXlab = xlabel('Frequency (Hz)');
% set(strXlab,'Position',[10.954461674932862,-7.044256077816259,-1]);
ylabel('PLV');
% set(gca,'FontSize',10,'Position',...
%     [0.859528795811527,0.657142857142854,0.119109947643973,0.21918367346938]);
xlim([4 100])
ylim([0 1])
title('g Hipp - C2')
%% Plot the PLV heatmap


axes(ha(8));
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end

strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
load(strVar_toLoad)
PLV_maint_topo = reshape(PLV_PSD_Granger_Topologies.PLV_maint,1,64);
ind = find(PLV_maint_topo>(max(PLV_maint_topo) - max(PLV_maint_topo)*0.15))
PLV_maint_topo(setdiff([1:length(PLV_maint_topo)],ind)) = 0;
elecColors = PLV_maint_topo;
cfg = [];
cfg.elecColorScale =[0.3 0.6]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
% 
% PLVmap = FigureMaint_Data.PLV_maps.GridHipp.Freq16_29.PLV_In_Band_To_Plot;
% 
% axes(ha(8));
% 
% imagesc(PLVmap,[0 0.6])
% ax1 = gca;
% ax1.XTick = 1:8;
% ax1.YTick = 1:8;
% ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
% ax1.YTickLabel = 1:8;
% ax1.YDir = 'normal'
% 
% colorbar('Ticks',[0 0.6],'FontSize',10)%,'Location')%,'northoutside');
% colormap(gca,'parula')
% % set(gca, 'Position',...
% %     [0.04869109947644,0.246938775510192,0.119109947643973,0.219183673469379]);
% xlabel('Grid column')
% ylabel('Grid row')
% title('h [-2 0] s, [16 29] Hz')
%% Granger spectra
Granger_enc = MainGranger_Fig_Data.GrangerSpectra.P42.Grng_Encoding{1}.grangerspctrm;
Granger_maint = MainGranger_Fig_Data.GrangerSpectra.P42.Grng_Maint{1}.grangerspctrm;
freq = MainGranger_Fig_Data.GrangerSpectra.P42.freq;
significance_bar_enc = MainGranger_Fig_Data.GrangerSpectra.P42.significance.Enc_CortexHipp;
significance_bar_maint = MainGranger_Fig_Data.GrangerSpectra.P42.significance.Maint_HippCortex;
light_blue = [0.30,0.75,0.93];
orange       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',orange};

yshift_1= -0.02;
yshift_2 = -0.03;
axes(ha(9));

semilogx(freq,Granger_enc(1,:)*100,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(freq,Granger_enc(2,:)*100,'Color',Colors{2},'LineWidth',3);
semilogx(freq,Granger_maint(1,:)*100,'Color',Colors{3},'LineWidth',3);
semilogx(freq,Granger_maint(2,:)*100,'Color',Colors{4},'LineWidth',3);
xlim([4 30])
plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',3);
plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',3);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
ylim([-0.05 0.2]*100)
strXlab = xlabel('Frequency (Hz)');
% set(strXlab,'Position',[10.954461674932862,-7.044256077816259,-1]);
ylabel('Granger (%)');
set(gca,'YTick',[0 10 20],'YTickLabel',[0 10 20],'XTick',[4 10 20 30],'XTickLabel',[4 10 20 30]);
title('i Hipp - C2')

% set(gca, 'Position',...
%     [0.250000000000007,0.246938775510192,0.119109947643973,0.219183673469379])




%% Granger heatmap encoding

axes(ha(10));
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end

strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
load(strVar_toLoad)
Granger_topo = reshape(PLV_PSD_Granger_Topologies.Granger_enc,1,64);
ind = find(Granger_topo<(min(Granger_topo) - min(Granger_topo)*0.2))
Granger_topo(setdiff([1:length(Granger_topo)],ind)) = 0;
elecColors = Granger_topo;
cfg = [];
cfg.elecColorScale =[-10 0]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_neg(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)


% for i =1:2
%     fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\210219_Supplementary Figure S42\Granger_heatmap_encoding_AHL2_3_grid_9_18Hz.fig';
%     H_sup1 = openfig(fig_path,'invisible');
%     ax_sup1 = gca;
%     figS1a = get(ax_sup1,'children');
%     axes(ha(10));
%     xlim([1 8])
%     ylim([1 8])
%     colormap(gca,'bluewhitered(256)');
%     colorbar('Ticks',[-10 0],'FontSize',10)%,'Location')%,'northoutside');
%     copyobj(figS1a,gca)
%     %     set(gca,'Position',...
%     %         [0.476171204188482,0.246938775510192,0.1191,0.2192]);
%     set(gca,'clim',[-10 0]);
% end
% ax1 = gca;
% ax1.XTick = 1:8;
% ax1.YTick = 1:8;
% ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
% ax1.YTickLabel = 1:8;
% ax1.YDir = 'normal'
% xlabel('Grid column')
% ylabel('Grid row')
% title('j [-5 -3] s, [8 18] Hz')

%% Granger heatmap maintenance

axes(ha(11));
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end

strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
load(strVar_toLoad)
Granger_topo = reshape(PLV_PSD_Granger_Topologies.Granger_maint,1,64);
ind = find(Granger_topo>(max(Granger_topo) - max(Granger_topo)*0.3))
Granger_topo(setdiff([1:length(Granger_topo)],ind)) = 0;
elecColors = Granger_topo;
cfg = [];
cfg.elecColorScale =[0 5]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)



% for i =1:2
%     
%     fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\210219_Supplementary Figure S42\Granger_heatmap_maintenance_AHL2_3_grid_9_18Hz.fig';
%     H_sup1 = openfig(fig_path,'invisible');
%     ax_sup1 = gca;
%     figS1a = get(ax_sup1,'children');
%     axes(ha(11));
%     xlim([1 8])
%     ylim([1 8])
%     copyobj(figS1a,gca);
%     colormap(gca,'bluewhitered(128)');
%     set(gca,'clim',[0 5]);
%     colorbar('Ticks',[0 5],'FontSize',10)%,'Location')%,'northoutside');
%     % set(gca,'Position',...
%     %     [0.670000000000008,0.246938775510192,0.1191,0.2192]);
% end
% ax1 = gca;
% ax1.XTick = 1:8;
% ax1.YTick = 1:8;
% ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
% ax1.YTickLabel = 1:8;
% ax1.YDir = 'normal'
% xlabel('Grid Column')
% ylabel('Grid contact')
% title('k [-2 0] s, [8 14] Hz')

%% TFR Granger
freq = MainGranger_Fig_Data.TFR_granger.P42.freq;
time = MainGranger_Fig_Data.TFR_granger.P42.time;
clim = MainGranger_Fig_Data.TFR_granger.P42.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P42.spectra;
cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap

axes(ha(12));
contourf(time,freq,TFR_Grng*100,100,'LineColor','none')
Xlab = xlabel('Time (s)');
% set(Xlab,'Position',[-1.999996185302734,4.372218008848717,1])
Ylab = ylabel('Frequency (Hz)');
% set(Ylab,'Position',[-7.76619346757089,12.247459177864897,1]);
ylim([5 30]);
cbh = colorbar%('Location')%,'northoutside')
colormap(gca,cmap);
set(cbh, 'Ticks',[-20  20],'TickLabels', [-20  20]);
set(ha(12),'clim',clim,'yscale','log')
set(ha(12),'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30]);
set(ha(12),'FontSize',10);
set(ha(12),'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
set(ha(12),'box','off');
set(ha(12),'TickDir','out','XMinorTick','off','YMinorTick','off')%,'TickLength',[0.02 0.025]);
line([-5 -5],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','-')
line([-3 -3],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','-')
line([0 0],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','-')
% set(gca,'Position',...
%     [0.859528795811527,0.246899999999999,0.119109947643973,0.2192]);
title('l Hipp - C2')



%% Axes properties

for i =1:length(ha)
    if i == 2||i == 3||i == 4||i == 6||i == 8||i == 9||i == 10
        set(ha(i),'box','off','FontSize',10,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025])
    else
        set(ha(i),'box','off','FontSize',10,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025])
    end
end

%% Colorbar labels

% Create text
text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
    'Position',[-12.241903233750968,8.449045717183134,0]);

% Create text
text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
    'Position',[3.758096694946289,8.449045730256017,0]);

% % Create text
text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
    'Position',[-28.3167,0.031,0]);
% 
% % Create text
text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
    'Position',[-12.615735080754646,0.03102318376432,0]);

text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
    'Position',[3.758096694946289,0.031023183044308,0]);
% % Create text
text('Parent',ha(9),'Rotation',90,'String','PLV',...
    'Position',[1.45359248910069 5.78674467180365 0]);

% % Create text
text('Parent',ha(9),'Rotation',90,'String','\DeltaGranger',...
    'Position',[1.38218488267912 -44.2132553281963 0]);

% % Create text
text('Parent',ha(9),'Rotation',90,'String','\DeltaGranger',...
    'Position',[0.093362334823627,-44.48207148300704,0]);

text('Parent',ha(12),'Rotation',90, 'String','\DeltaGranger',...
    'Position',[3.758625984191895,6.862013213252109,0]);
%% Text boxes
% Create textbox
annotation(fig,'textbox',...
    [0,0.945899999999999,0.04908376856078,0.069387753642336],...
    'String',{'a'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0,0.482469389214803,0.04908376856078,0.069387753642335],...
    'String',{'d'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.329842931937173 0.9459 0.04908376856078 0.0693877536423351],...
    'String',{'b'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.329842931937173 0.482469389214804 0.04908376856078 0.069387753642335],...
    'String',{'e'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.654450261780106,0.482469389214804,0.0425,0.0694],...
    'String',{'f'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.654450261780106 0.945938776969908 0.0477748680762283 0.0693877536423353],...
    'String',{'c'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
% annotation(fig,'textbox',...
%     [0.873036649214662 0.57246938903234 0.0798429299477508 0.0632653048452066],...
%     'Color',[1 1 1],...
%     'String',{'Hipp'},...
%     'FontSize',12,...
%     'EdgeColor','none');

% Create textbox
% annotation(fig,'textbox',...
%     [0.539267015706807 0.570428572705812 0.0916230343087182 0.0632653048452066],...
%     'Color',[1 1 1],...
%     'String',{'ECoG'},...
%     'FontSize',12,...
%     'EdgeColor','none');

% Create textbox
% annotation(fig,'textbox',...
%     [0.0877,0.4337,0.2029,0.0633],...
%     'Color',[1 1 1],...%[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String',{'PLV Hipp - ECoG'},...
%     'FontSize',12,...
%     'FitBoxToText','on',...
%     'EdgeColor','none');

% Create textbox
% annotation(fig,'textbox',...
%     [0.802356020942414,0.433693878828259,0.157068058536315,0.063265304845207],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String',{'Hipp - ECoG'},...
%     'FontSize',12,...
%     'FitBoxToText','on',...
%     'EdgeColor','none');

% annotation(fig,'textbox',...
%     [0.500000000000009,0.433693878828258,0.157068058536315,0.063265304845207],...
%     'Color',[0.149019607843137 0.149019607843137 0.149019607843137],...
%     'String',{'Hipp - ECoG'},...
%     'FontSize',12,...
%     'FitBoxToText','on',...
%     'EdgeColor','none');

%% Xticks
annotation(fig,'textbox',...
    [0.412303664921466 0.539816327564568 0.198952879581152 0.0551020397823684],...
    'String','-5       -3            0',...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.73691099476441 0.539816327564567 0.198952879581153 0.0551020397823684],...
    'String','-5       -3            0',...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.392670157068064 0.0704285724625266 0.298429319371727 0.0551020397823685],...
    'String','4                 10           20    30',...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.73036649214661 0.0683877561359955 0.198952879581153 0.055102039782368],...
    'String','-5       -3            0',...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.0392670157068063 0.0683877561359958 0.244764397905759 0.0551020397823686],...
    'String',{'A   B   C  D   E   F   G   H'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%% Save figure
set(gcf,'color','white')
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42.tiff');
