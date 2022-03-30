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
addpath(genpath('F:\Vasileios\Toolboxes\iELvis-iElvis_pm\'))

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
set(fig,'Position',[8.440208333333334,12.805833333333336,27.33145833333334,12.329583333333336])%[8.440208333333334,10.027708333333335,29.5275,15.107708333333337]);%[4.8948,2.8581,19,20]);
h=gcf;

ha = tight_subplot(2,5,[.1 .08],[.12 .04],[.08, .05])
% ha = tight_subplot(2,5,[.01 .08],[.1 .04],[.03 .03])%[0.05 0.05],[.1 .05],[.1 .1])%[0.05 0.03],[0.05 0.05])
set(ha(1:size(ha,1)),'Fontsize',20)


%% Plot the electrode locations
myPatient='USZ5';
globalFsDir = 'F:\Vasileios\Task Analysis\Electrodes Placement\By Piere\';
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
elec_2keep={'T5'}; % for USZpostopMR

elec_2keep_mask=ismember(tenten_snap.label,elec_2keep);
elec=tenten_snap;
elec.chanpos=elec.chanpos(elec_2keep_mask,:);
elec.chantype=elec.chantype(elec_2keep_mask);
elec.chanunit=elec.chanunit(elec_2keep_mask);
elec.elecpos=elec.elecpos(elec_2keep_mask,:);
elec.label=elec.label(elec_2keep_mask);


axes(ha(1));
% figure;
clear elecNames elecColors elecColorsEdge
GridLines = {'TLLS'}
ind = 4;

for i =1:4
        elecNames{i} = sprintf('%s%d',GridLines{1},i)
        elecColorsEdge(i,:) = [1 1 1];
        if ismember(i,ind)
            elecColorsEdge(i,:) = [1 1 0];
        end
end
strVarToLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Strip.mat';
load(strVarToLoad);
elecColors = DeltaGranger_S38_37.Granger_S38*100;
cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.view = 'l'
cfg.elecColorScale = [0 10]%[0 8];
cfg.elecCmapName = 'bluewhitered_pos(128)';
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.elecColors = elecColors';
cfg.ignoreDepthElec      ='n'
cfg.elecNames = elecNames'
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:\Vasileios\Task Analysis\Electrodes Placement\By Piere\'
cfg.elecSize = 6;
cfg. elecCbar = 'n';
cfg.opaqueness    = 0.3;
cfg.title = [];
cfg.axis = gca;
cfgOut = plotPialSurf('USZ5',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
cb = get(gca,'colorbar');
cb.Label.String = '\DeltaGranger (%)'
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment = 'bottom';
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.01;
CbPos(2) = 0.62214592274678;%CbPos(2)+0.03;
CbPos(3) = 0.009162955792191;
CbPos(4) = 0.27;
set(cb,'Position',CbPos)
set(gca, 'Position', [0.08,0.59,0.11,0.37])
set(cb.Label,'Position',[3.243,5.0000,0])
Axes_objs = get(gca,'Children');
indObj = [1];
set(Axes_objs(indObj),'LineWidth',0.8,'MarkerEdgeColor',[1 1 0]);
indWhite = setdiff([1:4],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none'); 
% 
hold on;
plot3(elec.elecpos(:,1),elec.elecpos(:,2)+9,elec.elecpos(:,3)+7,'Marker','o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',5,'Color','none', 'LineWidth',1.5);
hold off;

%%


% 
% imshow(I2);
set(gca,...ok
    'Position',[0.02,0.58,0.11,0.37]);
% 
% Im2Handle = imhandles(gca);
% set(Im2Handle,'XData',[1,1250],'YData',[1,1250]);
%%
axes(ha(6));
% figure
clear elecNames elecColors
GridLines = {'TOL'}
ind = 6
for i =1:6
    
        elecNames{i} = sprintf('%s%d',GridLines{1},i)
        elecColorsEdge(i,:) = [1 1 1];

         if ismember(i,ind)
            elecColorsEdge(i,:) = [1 1 0];
        end
end
strVarToLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Strip.mat';
load(strVarToLoad);
DeltaGranger_S38_37.Granger_S37 = DeltaGranger_S38_37.Granger_S37
elecColors = DeltaGranger_S38_37.Granger_S37*100-1.2;
cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.elecColorScale = [0 10]
cfg.view = 'l'
cfg.elecCmapName = 'bluewhitered_pos(128)'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.elecColors = elecColors';
cfg.ignoreDepthElec      ='n'
cfg.elecNames = elecNames'
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:\Vasileios\Task Analysis\Electrodes Placement\By Piere\'
cfg.elecSize = 7;
cfg. elecCbar = 'n';
% cfg.opaqueness    =0.5
cfg.title = [];
cfg.axis = gca;
cfgOut = plotPialSurf('USZ4',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
cb = get(gca,'colorbar');
cb.Label.String = '\DeltaGranger (%)'
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment = 'bottom';
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.01;
CbPos(2) = 0.211459227467811;%CbPos(2)+0.07;
CbPos(3) = 0.009162955792191;
CbPos(4) = 0.27;

set(cb,'Position',CbPos)
set(cb.Label,'Position',[3.1318519115448,5.060004351139069,0]);
set(gca,... % ok
    'Position',[0.02,0.16,0.11,0.37]);
% 
% Im3Handle = imhandles(gca);
% set(Im3Handle,'XData',[1,1100],'YData',[1,1100]);
Axes_objs = get(gca,'Children');
indObj = [1];
set(Axes_objs(indObj),'LineWidth',1.5,'MarkerEdgeColor',[1 1 0]);
indWhite = setdiff([1:6],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none'); 


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
set(strYlab,'Position',[-8.06324444278594,20.000030697615575,1]);
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
set(strYlab,'Position',[-7.934212184721424,20.000030697615575,1]);
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
set(strYlab,'Position',[-8.141018593010783,20.000030697615575,1])
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
set(strYlab,'Position',[-8.410930846866808,20.000030697615575,1]);
set(strXlab,'Position',[-2.0000,3.0081,1]);

%% Plot the Granger spectra

light_blue = [0.30,0.75,0.93];
orange       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',orange};

yshift_1= -0.005;
yshift_2 = -0.008;
gdata.Enc = MainGranger_Fig_Data.GrangerSpectra.P38.Enc;
gdata.Maint = MainGranger_Fig_Data.GrangerSpectra.P38.Maint;
significance_bar_enc = [6:8];%[9:16];
% significance_bar_enc2 = [20:22]
significance_bar_maint = [6:10]%[9:16];


axes(ha(4));
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(1,:)*100,'Color',Colors{1},'LineWidth',2);
hold on;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(2,:)*100,'Color',Colors{2},'LineWidth',2);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(1,:)*100,'Color',Colors{3},'LineWidth',2);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(2,:)*100,'Color',Colors{4},'LineWidth',2);
xlim([4 20])
plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',2);
% plot(significance_bar_enc2,yshift_1*100*ones(length(freq(significance_bar_enc2)),1)','Color',Colors{1},'LineWidth',3);
plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',2);
ylim([-0.01 0.2]*100)


strYlab = ylabel('Granger (%)','VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)','VerticalAlignment','cap')

set(ha(4),'FontSize',10);
set(ha(4),'YTick',[0 0.1 0.2]*100,'YTickLabel', [0 0.1 0.2]*100);
set(ha(4),'XTick',[4 10 20],'XTickLabel', [4 10 20]);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
set(gca,...
    'Position',[0.653141361256544,0.620408163265306,0.11910994764398,0.285714285714285]);
set(strYlab,'Position',[2.551430067086884,8.426702990343697,-1]);
set(strXlab,'Position',[10.954461674932862,-2.451128015392704,-1]);


gdata.Enc = MainGranger_Fig_Data.GrangerSpectra.P37.Enc;
gdata.Maint = MainGranger_Fig_Data.GrangerSpectra.P37.Maint;
significance_bar_enc = [6:9]%[10:16]%[10:18];
significance_bar_maint = [6:8]%[9:18];


axes(ha(9));

yshift_1= -0.003;
yshift_2 = -0.002;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(1,:)*100*2,'Color',Colors{1},'LineWidth',2);
hold on;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(2,:)*100,'Color',Colors{2},'LineWidth',2);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(1,:)*100*2,'Color',Colors{3},'LineWidth',2);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(2,:)*100,'Color',Colors{4},'LineWidth',2);
xlim([4 20])
plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',2);
plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',2);
ylim([-0.005 0.05]*100)


strYlab = ylabel('Granger (%)');%,'VerticalAlignment','bottom')
strXlab = xlabel('Frequency (Hz)');%,'VerticalAlignment','cap')

set(ha(9),'FontSize',10);
set(ha(9),'YTick',[0 0.05 0.1]*100,'YTickLabel', [0 0.05 0.1]*100);
set(ha(9),'XTick',[4 10 20],'XTickLabel', [4 10 20]);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
set(gca,...
    'Position',[0.653141361256544,0.2,0.119109947643981,0.289999999999999])
set(strYlab,'Position',[2.591273957556311,2.490361802224759,-1]);
set(strXlab,'Position',[10.954461674932862,-1.125873371075701,-1]);

%% Plot the TFR Granger

freq = MainGranger_Fig_Data.TFR_granger.P38.freq;
time = MainGranger_Fig_Data.TFR_granger.P38.time;
clim = MainGranger_Fig_Data.TFR_granger.P38.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P38.spectra;
cmap = MainGranger_Fig_Data.TFR_granger.P38.cmap

% cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap
axes(ha(5));
contourf(time,freq,TFR_Grng*100,100,'LineColor','none')
hold on;
% contour(time,freq,TFR_Grng,'k')
Xlab = xlabel('Time (s)');
Ylab = ylabel('Frequency (Hz)');
ylim([4 20]);
cbh = colorbar
colormap(gca,cmap);
set(cbh, 'Ticks',[-10  10],'TickLabels', [-10  10], 'Position',[0.952006980802793,0.622448979591837,0.01086387434555,0.285714285714286]);
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

set(Ylab,'Position',[-8.270969697407313,8.62261504662936,1]);
set(Xlab,'Position',[-2.214281899588449,3.37894189560413,1]);



freq = MainGranger_Fig_Data.TFR_granger.P37.freq;
time = MainGranger_Fig_Data.TFR_granger.P37.time;
clim = MainGranger_Fig_Data.TFR_granger.P37.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P37.spectra;
% cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap
axes(ha(10));
contourf(time,freq,TFR_Grng,100,'LineColor','none')
hold on;
% contour(time,freq,TFR_Grng,'k')
Xlab = xlabel('Time (s)');
Ylab = ylabel('Frequency (Hz)');
ylim([5 20]);
cbh = colorbar
colormap(gca,cmap);
set(cbh, 'Ticks',[-5  5],'TickLabels', [-5  5],'Position',[0.955933682373473,0.2,0.01151832460733,0.285714285714286]);
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
 
set(Ylab,'Position',[-7.893608307638088,9.942844955240146,1]);
set(Xlab,'Position',[-2.09090518951416,4.18133783467852,1]);


%% Axes properties

for i =1:length(ha)
    if ismember(i,[2 3 5 7 8 10])
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
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.00785340314136404 0.520000000000001 0.0392670149266414 0.0429556683917689],...
    'String','f',...
    'FontWeight','bold',...
    'FontSize',16,...
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
    [0.157068062827228 0.945450542545293 0.039267014926641 0.032296650111675],...
    'String',{'b',''},...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.382198952879591 0.949532175198355 0.0392670149266411 0.0322966501116756],...
    'String','c',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.592931937172797 0.953613807851415 0.039267014926641 0.032296650111675],...
    'String','d',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.778795811518335 0.953867689308302 0.0392670149266412 0.032296650111675],...
    'String','e',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.153141361256552 0.534105947290921 0.039267014926641 0.032296650111675],...
    'String','g',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.380890052356043 0.535092179104218 0.0392670149266411 0.032296650111675],...
    'String','h',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.598167539267017 0.530771312001399 0.039267014926641 0.032296650111675],...
    'String','i',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.781413612565449 0.530278196094747 0.0392670149266412 0.0322966501116748],...
    'String','j',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.645806196560623 0.574615311535837 0.192408376963349 0.0551020397823681],...
    'String','4              10         20',...
'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');
% Create textbox
annotation(fig,'textbox',...
    [0.447875855917042 0.574510205115589 0.134816753926702 0.0551020397823681],...
    'String','-5     -3         0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.643583726552564 0.15002040919722 0.192408376963349 0.055102039782368],...
    'String','4              10         20',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.221831396380188 0.14797959287069 0.134816753926701 0.055102039782368],...
    'String',' -5     -3          0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.450493656964163 0.14797959287069 0.12696335078534 0.055102039782368],...
    'String','-5    -3         0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.83850093510996 0.150020409197221 0.134816753926701 0.055102039782368],...
    'String','-5    -3         0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.227353360060415 0.576551021442121 0.134816753926702 0.0551020397823681],...
    'String','-5      -3        0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.839237112461546 0.574510205115588 0.134816753926701 0.0551020397823679],...
    'String','-5     -3        0',...
    'FontSize',10.5,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.520383623158291 0.209912293464125 0.077225128978647 0.037081338964296],...
    'Color',[1 1 1],...
    'String','Hipp',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.295229138431752 0.2115 0.0916 0.0371],...
    'Color',[1 1 1],...
    'String','ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.521351677369327 0.636232692868523 0.0772251289786472 0.0370813389642958],...
    'Color',[1 1 1],...
    'String','Hipp',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.861474483408771,0.878184549356223,0.145287958115183,0.037081338964296],...
    'String','Hipp-ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.863410591830842,0.459729613733906,0.145287958115182,0.037081338964296],...
    'String','Hipp-ECoG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.305383090981894,0.636232692868521,0.077225128978647,0.037081338964296],...
    'Color',[1 1 1],...
    'String','EEG',...
    'FontSize',12,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');


h2 = text(-23.346694160910214,51.82388849199689,0.237200021743774,' Rel. PSD')
set(h2,'Rotation',90,'FontSize',12);

h3 = text(-38.40551769032198,51.82388849199689,0.237200021743774,' Rel. PSD')
set(h3,'Rotation',90,'FontSize',12);

h5 = text(3.476835250854492,46.69413634654633,0.237200021743774,'\DeltaGranger (%)')
set(h5,'Rotation',90,'FontSize',12);

h6 = text(3.591426368521049,6.030939630001639,0.237200021743774,'\DeltaGranger (%)')
set(h6,'Rotation',90,'FontSize',12);


h7= text(-23.346694160910214,7.0781,0.2372,' Rel. PSD')
set(h7,'Rotation',90,'FontSize',12);


h8 = text(-38.40551769032198,7.0781,0.2372,' Rel. PSD')
set(h8,'Rotation',90,'FontSize',12);



%% XTicks
for i = 2:size(ha,1)
    set(ha(i),'XTickLabel',[]);
end


%% Save Figure
set(gcf,'color','white')
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_r2_1','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_r2_1','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_r2_1','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_r2_1.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_r2_1.tiff');
