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

% strPaths.FigureVars2 = [strPaths.FigureData,'FigureMaintData.mat'];
% load(strPaths.FigureVars2)
%
% strPaths.FigureVars3 = [strPaths.FigureData,'FigurePSD_Data.mat'];
% load(strPaths.FigureVars3)
%
% ElectrodeImagesPath = [strPaths.FigureData,'Main Granger Figure\'];
% strPaths.ImageP42 = [ElectrodeImagesPath,'P42_electrodes.jpg'];

strPaths.FigureDataReref = [strPaths.FigureData 'Main Granger Figure\Reref_with_2_WhiteMatter_contacts\'];
cd(strPaths.FigureDataReref)
files = dir('*.mat')
for i = 1:length(files)
    load(files(i).name)
end
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
PSD_zscore.enc([3 4],4) = [0 0];
PSD_zscore.enc([2 4],6) = [0 0];

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
GridLines = {'A','B','C','D','E','F','G','H'}


% strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\PSD_PLV_Granger_topos.mat';
% load(strVar_toLoad)
elecColors = reshape(PSD_zscore.enc,1,64)
ind = find(elecColors>=2);
% elecColors(59) = elecColors(59)+2.5;
% ind = [ind 59]
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColorsEdge((i-1)*8+j,:) = [1 1 1];
        if ismember((i-1)*8+j,ind)
            elecColorsEdge((i-1)*8+j,:) = [1 1 0]
        end
    end
end
% figure;
cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.elecColorScale =[2 3]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName = 'custom_blue_white_map';%'bluewhitered_neg(128)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 6;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
set(gca,... ok
    'Position',[0.01,0.825,0.236666666666667,0.135]);
firstPos =get(gca,'Position');
cb = get(gca,'colorbar');
cb.Label.String = 'z-score'
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment = 'bottom';
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.05;
CbPos(3) = 0.01508875739645;
set(cb,'Position',CbPos)


Axes_objs = get(gca,'Children');
indObj = [4 6 7 10 15 47 55]%[2 6 7 46 47 55 ];
set(Axes_objs(indObj),'LineWidth',1.5);
indWhite = setdiff([1:64],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none');
%
% hold on;
% plot3(elec.elecpos(:,1),elec.elecpos(:,2),elec.elecpos(:,3),'Marker','o','MarkerFaceColor','white','MarkerEdgeColor','black','MarkerSize',5,'Color','none', 'LineWidth',1.5);
% hold off;

% title('a \color{cyan}C2, \color{magenta}H3')
%% TFR PSD
clim = FigurePSD_Data.TFR_PSD.ECoG.clim;
freq = FigurePSD_Data.TFR_PSD.ECoG.freq;
time = FigurePSD_Data.TFR_PSD.ECoG.time;
TFR_PSD_ECoG = FigurePSD_Data.TFR_PSD.ECoG.TFR_baselinedPSD;

axes(ha(2));

contourf(time,freq,TFR_PSD_ECoG,100,'LineColor','none');
set(gca,'clim',[0 1],'yscale','log');
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]); %background color of grid;
set(gca,'FontSize',10);
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
colormap(gca,jet);
cb = colorbar('Ticks',[-1 0 1 2 3 4],'TickLabels',[-1 0 1 2 3 4])%,'Location')%,'northoutside');
cb.Label.String = 'Rel. PSD';
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment ='bottom'
line([-5 -5],get(ha(2),'YLim'),'Color',[1 1 1]);
line([-3 -3],get(ha(2),'YLim'),'Color',[1 1 1]);
line([0 0],get(ha(2),'YLim'),'Color',[1 1 1]);
ylabel('Frequency (Hz)')%,'VerticalAlignment','middle');
strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap');
% title('b H3')

ax2Pos = get(gca,'Position');
newPos = ax2Pos;
newPos(1) = ax2Pos(1)+0.03;
set(gca,'Position',newPos);
ax2Pos = newPos;
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
    cb = colorbar('Ticks',[0 1],'TickLabels',[0 1])%,'Location')%,'northoutside');
    line([-5 -5],get(ha(3),'YLim'),'Color',[1 1 1]);
    line([-3 -3],get(ha(3),'YLim'),'Color',[1 1 1]);
    line([0 0],get(ha(3),'YLim'),'Color',[1 1 1]);
    ylabel('Frequency (Hz)')%,'VerticalAlignment','middle');
    strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap');
end
% title('c C2')
ax3Pos = get(gca,'Position');
newPos = ax3Pos;
newPos(1) = ax3Pos(1)+0.029;
set(gca,'Position',newPos);
ax3Pos = newPos;
cb.Label.String = 'Rel. PSD';
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment ='bottom'
cb.Label.Position = [2.84606053612449,0.500000476837158,0]
%% PSD Topo Grid maint
axes(ha(4));
GridLines = {'A','B','C','D','E','F','G','H'}
PSD_zscore.maint(8,7) =0;
%
% strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
% load(strVar_toLoad)
PSD_maint_topo = reshape(PSD_zscore.maint,1,64);
ind = find(PSD_maint_topo>(max(PSD_maint_topo) - max(PSD_maint_topo)*0.2))
elecColors = PSD_maint_topo;

for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColorsEdge((i-1)*8+j,:) = [1 1 1];
        if ismember((i-1)*8+j,ind)
            elecColorsEdge((i-1)*8+j,:) = [1 1 0]
        end
    end
end

cfg = [];
cfg.elecColorScale =[2 3]
cfg.elecColorsEdge = elecColorsEdge;
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(512)';
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 6;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)

colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
oldPos = get(gca,'Position');
newPos = [firstPos(1) oldPos(2) firstPos(3) oldPos(4)];
set(gca, 'Position', newPos);
cb = get(gca,'colorbar');
cb.Label.String = 'z-score'
cb.Label.FontSize =11;
cb.Label.VerticalAlignment = 'bottom';
% cb.Label.Position = [2.846060536124492,0.500000476837158,0]
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.05;
CbPos(3) = 0.01508875739645;
set(cb,'Position',CbPos)

Axes_objs = get(gca,'Children');
indObj = [6:7 47 55];
set(Axes_objs(indObj),'LineWidth',1.5);
indWhite = setdiff([1:64],indObj);
indWhite = setdiff([1:64],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none');
%% PSD Hipp

clim = FigurePSD_Data.TFR_PSD.iEEG_Hipp.clim
freq = FigurePSD_Data.TFR_PSD.iEEG_Hipp.freq;
time = FigurePSD_Data.TFR_PSD.iEEG_Hipp.time
TFR_PSD_Hipp = FigurePSD_Data.TFR_PSD.iEEG_Hipp.TFRbaselinedPSD;
% var = TFR_PSD_Hipp(:,1:5);
% [ind1,ind2] = find(~isnan(TFR_PSD_Hipp(:,1:5))==1);
% TFR_PSD_Hipp = (TFR_PSD_Hipp-3.5)./3.5;
% TFR_PSD_Hipp(:,1:5) = var;
% TFR_PSD_Hipp(ind1,ind2) = 0;
axes(ha(5))

contourf(time,freq,TFR_PSD_Hipp,100,'LineColor','none')
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'FontSize',10)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap(gca,jet);
cb = colorbar('Ticks',[-1:1:4],'TickLabels',[-1:1:4])%,'Location')%,'northoutside');
cb.Label.String = 'Rel. PSD';
cb.Label.FontSize = 11;
cb.Label.VerticalAlignment ='bottom'

ylabel('Frequency (Hz)')%,'VerticalAlignment','middle')
strXlab = xlabel('Time (s)')%,'VerticalAlignment','cap')
% set(strXlab,'Position',[-1.999996185302734,3.199670443946102,1])
line([-5 -5],get(ha(5),'YLim'),'Color',[1 1 1])
line([-3 -3],get(ha(5),'YLim'),'Color',[1 1 1])
line([0 0],get(ha(5),'YLim'),'Color',[1 1 1])
% title('e Hipp')

oldPos = get(gca,'Position')
newPos = [ax2Pos(1) oldPos(2) ax2Pos(3) oldPos(4)]
set(gca,'Position',newPos)



%% PLV Hipp ECoG
% for i =1:2
%     fig_path = 'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\210219_Supplementary Figure S42\PLV_Hipp_ECoG_C2.fig';
%     H_sup1 = openfig(fig_path,'invisible');
%     ax_sup1 = gca;
%     figS1a = get(ax_sup1,'children');
%     axes(ha(6));
%
%     copyobj(figS1a,gca)
%     set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off',...
%         'TickLength',[0.02 0.025],'YTick',[0 1],'YTickLabel',[0 1],'XTick',[4 10 20 100],'XTickLabel',[4 10 20 100]);
%
% end
PLV_maint = PLV_after_Reref_2WM.PLV_maint2;
PLV_fix =  PLV_after_Reref_2WM.PLV_fix2;
PLV_enc = PLV_after_Reref_2WM.PLV_enc2%*1.35;
freq =  PLV_after_Reref_2WM.freq;
prc =  PLV_after_Reref_2WM.Prc
yshift = 0.05;
axes(ha(6));

semilogx(freq,PLV_fix,'LineWidth',2,'Color','k')
hold on;
semilogx(freq,PLV_enc,'LineWidth',2,'Color','b')
semilogx(freq,PLV_maint,'LineWidth',2,'Color','r')

indMaxPLVbar = [4:6]%[16:29]; %Difference_Bar(PLV_maint-PLV_fix,prc,freq,[0 1],0.05,'r')
semilogx(freq(indMaxPLVbar),ones(1,length(indMaxPLVbar))*yshift,'r','LineWidth',2)

indMaxPLVbar = [8:10]%[16:29]; %Difference_Bar(PLV_maint-PLV_fix,prc,freq,[0 1],0.05,'r')
semilogx(freq(indMaxPLVbar),ones(1,length(indMaxPLVbar))*yshift,'r','LineWidth',2)
strXlab = xlabel('Frequency (Hz)');
ylabel('PLV');
xlim([4 20])
ylim([0 1])
% title('f Hipp - C2')
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off',...
    'TickLength',[0.02 0.025],'YTick',[0 1],'YTickLabel',[0 1],'XTick',[4 10 20 100],'XTickLabel',[4 10 20 100]);

oldPos = get(gca,'Position')
newPos = [ax3Pos(1) oldPos(2) ax3Pos(3)+0.028 oldPos(4)]
set(gca,'Position',newPos)

xlab = get(gca, 'XLabel');
set(xlab, 'Position',[20.000030697615575,-0.323010737857511,-1]);
set(gca,'YTick',[0 1],'YTickLabel',[0 0.6]);
%% Plot the PLV heatmap


axes(ha(7));

%
strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid2.mat';
load(strVar_toLoad)
% figure
PLV_topo = PLV_PSD_Granger_Topologies.PLV_maint2;
PLV_maint_topo = reshape(PLV_topo,1,64);
ind = [find(PLV_maint_topo>(max(PLV_maint_topo) - max(PLV_maint_topo)*0.3)) 19 32]
% PLV_maint_topo(setdiff([1:length(PLV_maint_topo)],ind)) = 0;
elecColors = PLV_maint_topo;

GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColorsEdge((i-1)*8+j,:) = [1 1 1];
        if ismember((i-1)*8+j,ind)
            elecColorsEdge((i-1)*8+j,:) = [1 1 0]
        end
        
        
    end
end

cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.elecColorScale =[0.4 .8]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(128)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 6;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',[0 0.6])
caxis(cfgOut.elecCbarLimits)
oldPos = get(gca,'Position');
newPos = [firstPos(1) oldPos(2) firstPos(3) oldPos(4)];
set(gca, 'Position', newPos);
cb = get(gca,'colorbar');
CbPos = cb.Position;
cb.Label.String = 'PLV'
cb.Label.FontSize =11;
cb.Label.VerticalAlignment = 'bottom';
% cb.Label.Position = [3.027878717942671,0.650000333786011,0]
CbPos(1) = CbPos(1)-0.05;
CbPos(3) = 0.01508875739645;
set(cb,'Position',CbPos)

Axes_objs = get(gca,'Children');
indObj = [46 47 54 55 63 64 ];
set(Axes_objs(indObj),'LineWidth',1.5);
indWhite = setdiff([1:64],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none');
% title('g [-2 0] s, [16 29] Hz')
%% Granger spectra
Granger_enc1 = GrangerVarsP42_2WM_reref.Enc1;
Granger_enc2 = GrangerVarsP42_2WM_reref.Enc2;
Granger_maint1 = GrangerVarsP42_2WM_reref.Maint1;
Granger_maint2 = GrangerVarsP42_2WM_reref.Maint2;

freq = GrangerVarsP42_2WM_reref.freq;
Colors     = GrangerVarsP42_2WM_reref.Colors;
prc_enc = GrangerVarsP42_2WM_reref.Prc_enc;
prc_maint = GrangerVarsP42_2WM_reref.Prc_maint;
yshift_1= -0.005;
yshift_2 = -0.007;
axes(ha(8));

semilogx(freq,Granger_enc2*100,'Color',Colors{1},'LineWidth',2);
hold on;
semilogx(freq,Granger_enc1*100,'Color',Colors{2},'LineWidth',2);
semilogx(freq,Granger_maint2*100,'Color',Colors{3},'LineWidth',2);
semilogx(freq,Granger_maint1*100,'Color',Colors{4},'LineWidth',2);
xlim([4 20])
indCortexHipp_enc = [121:161];%Difference_Bar(Granger_enc2-Granger_enc1,prc_enc-0.019,freq,[0 20],yshift_1*100, Colors{1});
indHippCortex_maint = [101:161];%Difference_Bar(Granger_maint2-Granger_maint1,prc_maint-0.0028,freq,[0 20],yshift_2*100, Colors{3});

semilogx(freq(indCortexHipp_enc), ones(1,length(indCortexHipp_enc))*yshift_1*100,'color', Colors{1},'LineWidth',2);
semilogx(freq(indHippCortex_maint), ones(1,length(indHippCortex_maint))*yshift_2*100,'color', Colors{3},'LineWidth',2);


% plot(significance_bar_enc,yshift_1*100*ones(length(freq(significance_bar_enc)),1)','Color',Colors{1},'LineWidth',3);

% plot(significance_bar_maint,yshift_2*100*ones(length(freq(significance_bar_maint)),1)','Color',Colors{3},'LineWidth',3);
set(gca,'box','off','XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
ylim([-0.01 0.1]*100)
strXlab = xlabel('Frequency (Hz)');
ylabel('Granger (%)');
set(gca,'YTick',[0 10 20],'YTickLabel',[0 10 20],'XTick',[4 10 20 30],'XTickLabel',[4 10 20 30]);
% title('i Hipp - C2')

oldPos = get(gca,'Position')
newPos = [ax2Pos(1) oldPos(2) ax2Pos(3)+0.028 oldPos(4)]
set(gca,'Position',newPos)
xlab = get(gca, 'XLabel');
set(xlab, 'Position', [10.954461674932862,-13.599354350759135,-1]);
%% Granger heatmap encoding

axes(ha(9));
GridLines = {'A','B','C','D','E','F','G','H'}
% strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
% load(strVar_toLoad)

Granger_topo = reshape(GrangerHeatmaps.Enc,1,64)*100;
ind = find(Granger_topo<(min(Granger_topo) - min(Granger_topo)*0.4))
% Granger_topo(setdiff([1:length(Granger_topo)],ind)) = 0;
elecColors = Granger_topo;

for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColorsEdge((i-1)*8+j,:) = [1 1 1];
        if ismember((i-1)*8+j,ind)
            elecColorsEdge((i-1)*8+j,:) = [1 1 0]
        end
    end
end

cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.elecColorScale = GrangerHeatmaps.Enc_clim*100
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_neg(128)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 6;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
% oldPos = get(gca,'Position');

oldPos = get(gca,'Position')
newPos = [ax3Pos(1)-0.05 oldPos(2) firstPos(3) oldPos(4)]
set(gca,'Position',newPos)

% newPos = [firstPos(1) oldPos(2) firstPos(3) oldPos(4)];
% set(gca, 'Position', newPos);
cb = get(gca,'colorbar');
cb.Label.String = '\DeltaGranger (%)'
cb.Label.FontSize =11;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.Position = [2.755151358517733,-4.999995231628418,0];
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.055;
CbPos(3) = 0.01508875739645;
set(cb,'Position',CbPos)
set(cb.Label,'Position',[4.209696726365523,-4.569887704746698,0]);
% title('j [-5 -3] s, [8 18] Hz')
Axes_objs = get(gca,'Children');
indObj = [39 40 46 47 64];
set(Axes_objs(indObj),'LineWidth',1.5);
indWhite = setdiff([1:64],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none');
%% Granger heatmap maintenance

axes(ha(10));
GridLines = {'A','B','C','D','E','F','G','H'}

%
% strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\PLV_PSD_Granger_Topologies_Grid.mat';
% load(strVar_toLoad)

Granger_topo = reshape(GrangerHeatmaps.Maint,1,64)*100;
ind = find(Granger_topo>(max(Granger_topo) - max(Granger_topo)*0.5))
elecColors = Granger_topo;

for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)
        elecColorsEdge((i-1)*8+j,:) = [1 1 1];
        if ismember((i-1)*8+j,ind)
            elecColorsEdge((i-1)*8+j,:) = [1 1 0]
        end
    end
end

cfg = [];
cfg.elecColorsEdge = elecColorsEdge;
cfg.elecColorScale = GrangerHeatmaps.Maint_clim*100;
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='n';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_pos(128)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'marker'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 6;
cfg.elecCbar = 'n'
cfg.axis = gca;
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
colormap(gca,cfgOut.elecCmapName)
colorbar('Ticks',cfgOut.elecCbarLimits,'TickLabels',cfgOut.elecCbarLimits)
caxis(cfgOut.elecCbarLimits)
oldPos = get(gca,'Position');
newPos = [firstPos(1) oldPos(2) firstPos(3) oldPos(4)];
set(gca, 'Position', newPos);
cb = get(gca,'colorbar');
cb = get(gca,'colorbar');
cb.Label.String = '\DeltaGranger (%)'
cb.Label.FontSize =11;
cb.Label.VerticalAlignment = 'bottom';
CbPos = cb.Position;
CbPos(1) = CbPos(1)-0.05;
CbPos(3) = 0.01508875739645;
set(cb,'Position',CbPos);
set(cb.Label,'Position',[4.562424291263927,4.892477886651152,0]);
Axes_objs = get(gca,'Children');
indObj = [16 18 47 55 64];
set(Axes_objs(indObj),'LineWidth',1.5);
indWhite = setdiff([1:64],indObj);
set(Axes_objs(indWhite),'MarkerEdgeColor','none');
% title('k [-2 0] s, [8 14] Hz')
%% Scalar Product Granger

% strVar_toLoad = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\Main Granger Figure\ScalarProductGranger_maint.mat';
% load(strVar_toLoad)
hfodPath = 'F:\Vasileios\HFO Analysis\Geneva HFO analysis\hfod\';
addpath(genpath(hfodPath))
custom_map     = Visual.Utility.getMyColourMap();
axes(ha(11));
ActualDistofSP = Granger_scalarProduct.ActualDist;
RandomDistofSP =Granger_scalarProduct.RandomDist;

PercentileVal = prctile(RandomDistofSP(:), 95);
PersentageOfDataThatsucceed = (1 - mean(ActualDistofSP(:) < PercentileVal))*100;
binStarts = 0:0.01:1;
%     distScalarFig = figure('units','normalized','outerposition',[0 0 1 1]);
colormap(gca, custom_map)
hold on
hist1 = histogram(RandomDistofSP(:), binStarts,'Normalization','probability','EdgeAlpha',0,'facecolor',[0.8,0.8,0.8],'facealpha', 0.85);
hist2 = histogram(ActualDistofSP(:), binStarts,'Normalization','probability','EdgeAlpha',0,'facecolor',[1,0,0],'facealpha', 0.85);
PercentileLine = line([PercentileVal, PercentileVal], [0, 0.1]);
PercentileLine.LineWidth = 2;
hold off
ylabel('Prob. Density')
xlabel('Scalar Product')
ylim([0 nanmax([hist1.Values, hist2.Values])+0.01])
xlim([0 1])
oldPos = get(gca,'Position')
newPos = [ax2Pos(1) oldPos(2) ax2Pos(3)+0.028 oldPos(4)]
set(gca,'Position',newPos)
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 0.5 1], 'YTick', [0 nanmax([hist1.Values, hist2.Values])+0.01], 'YTickLabel',[0 round(nanmax([hist1.Values, hist2.Values])+0.01,1)]);

xlab = get(gca,'XLabel')
set(xlab, 'Position',[0.500000476837158,-0.044350191434253,-1]);

%% TFR Granger
freq = MainGranger_Fig_Data.TFR_granger.P42.freq;
time = MainGranger_Fig_Data.TFR_granger.P42.time;
clim = MainGranger_Fig_Data.TFR_granger.P42.clim*100;
TFR_Grng =MainGranger_Fig_Data.TFR_granger.P42.spectra;
cmap = MainGranger_Fig_Data.TFR_granger.P42.cmap

axes(ha(12));
contourf(time,freq,TFR_Grng,100,'LineColor','none')
Xlab = xlabel('Time (s)');
% set(Xlab,'Position',[-1.999996185302734,4.372218008848717,1])
Ylab = ylabel('Frequency (Hz)');
% set(Ylab,'Position',[-7.76619346757089,12.247459177864897,1]);
ylim([5 20]);
cbh = colorbar%('Location')%,'northoutside')
cbh.Label.String = '\DeltaGranger (%)';
cbh.Label.FontSize = 11;
cbh.Label.VerticalAlignment ='bottom'
cbh.Label.Position = [4.755151336843317,0.000019073486328,0];
colormap(gca,cmap);
% colormap(gca,bluewhitered(128));

set(cbh, 'Ticks',[-10  10],'TickLabels', [-10 10]);
set(ha(12),'clim',clim,'yscale','log')
set(ha(12),'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30]);
set(ha(12),'FontSize',10);
set(ha(12),'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0]);
set(ha(12),'box','off');
set(ha(12),'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025]);
line([-5 -5],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([-3 -3],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','--')
line([0 0],get(ha(12),'YLim'),'Color',[0 0 0],'LineStyle','--')
% set(gca,'Position',...
%     [0.859528795811527,0.246899999999999,0.119109947643973,0.2192]);
% title('l Hipp - C2')


oldPos = get(gca,'Position')
newPos = [ax3Pos(1) oldPos(2) ax3Pos(3) oldPos(4)]
set(gca,'Position',newPos)


%% Axes properties

for i =1:length(ha)
    if ismember(i,[2 3 5 12])
        set(ha(i),'box','off','FontSize',10.5,'TickDir','out','XMinorTick','off','YMinorTick','off','TickLength',[0.02 0.025])
    else
        set(ha(i),'box','off','FontSize',10.5,'TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025])
    end
    Xlab = get(ha(i),'XLabel')
    if ~isempty(Xlab)
        set(Xlab,'VerticalAlignment', 'middle');
    end
    if ~ismember(i,[6 8 11])
        Xlab = get(ha(i),'XLabel')
        Xlab.Position(2) = Xlab.Position(2)-0.18;
    end
    if i==8
        Xlab = get(ha(i),'XLabel');
        set(Xlab,'Position',[10.954461697121948,-4.620859645387179,-1]);
    elseif i == 11
        Xlab = get(ha(i),'XLabel');
        set(Xlab,'Position',[0.5,-0.0341,-1]);
        
    end
    
    
    Ylab = get(ha(i),'YLabel')
    if i== 12
        Ylab.Position = [-8.38641511016297,10.000006610368823,1];
    else
        Ylab.Position(1) = Ylab.Position(1)-0.01;
    end
    
    if ~isempty(Ylab)
        set(Ylab,'VerticalAlignment', 'baseline');
    end
    
    
    
end

% %% Colorbar labels
%
% % Create text
% text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
%     'Position',[-12.241903233750968,8.449045717183134,0]);
%
% % Create text
% text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
%     'Position',[3.758096694946289,8.449045730256017,0]);
%
% % % Create text
% text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
%     'Position',[-28.3167,0.031,0]);
% %
% % % Create text
% text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
%     'Position',[-12.615735080754646,0.03102318376432,0]);
%
% text('Parent',ha(3),'Rotation',90,'String','Rel. PSD',...
%     'Position',[3.758096694946289,0.031023183044308,0]);
% % % Create text
% text('Parent',ha(9),'Rotation',90,'String','PLV',...
%     'Position',[1.45359248910069 5.78674467180365 0]);
%
% % % Create text
% text('Parent',ha(9),'Rotation',90,'String','\DeltaGranger',...
%     'Position',[1.38218488267912 -44.2132553281963 0]);
%
% % % Create text
% text('Parent',ha(9),'Rotation',90,'String','\DeltaGranger',...
%     'Position',[0.093362334823627,-44.48207148300704,0]);
%
% text('Parent',ha(12),'Rotation',90, 'String','\DeltaGranger',...
%     'Position',[3.758625984191895,6.862013213252109,0]);
% %% Text boxes
% % Create textbox
annotation(fig,'textbox',...
    [0,0.945899999999999,0.04908376856078,0.069387753642336],...
    'String',{'a'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.323925772173859,0.9459,0.04908376856078,0.069387753642335],...
    'String',{'b'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.648533102016792,0.9473943810456,0.047774868076228,0.069387753642336],...
    'String',{'c'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0,0.715366041325435,0.04908376856078,0.069387753642335],...
    'String',{'d'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',...
    [0.321005917159763,0.715366041325436,0.04908376856078,0.069387753642335],...
    'String',{'e'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.652366863905327,0.715366041325436,0.04908376856078,0.069387753642335],...
    'String',{'f'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');



% Create textbox
annotation(fig,'textbox',[0,0.483924993290508,0.04908376856078,0.069387753642335],...
    'String',{'g'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.323964497041419,0.483924993290508,0.04908376856078,0.069387753642335],...
    'String',{'h'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.647928994082841,0.482469389214817,0.04908376856078,0.069387753642335],...
    'String',{'i'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0,0.251028341179885,0.04908376856078,0.069387753642335],...
    'String',{'j'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.323964497041419,0.249572737104193,0.04908376856078,0.069387753642335],...
    'String',{'k'},...
    'FontSize',16,...
    'FitBoxToText','off',...
    'FontWeight','bold',...
    'EdgeColor','none');

% Create textbox
annotation(fig,'textbox',[0.646449704142013,0.248117133028502,0.04908376856078,0.069387753642335],...
    'String',{'l'},...
    'FontSize',16,...
    'FontWeight','bold',...
    'FitBoxToText','off',...
    'EdgeColor','none');
%
%Create textbox
annotation(fig,'textbox',...
    [0.523924223179159,0.560824556426809,0.079842929947751,0.063265304845207],...
    'Color',[1 1 1],...
    'String',{'Hipp'},...
    'FontSize',12,...
    'FontWeight','bold', ...
    'EdgeColor','none');
%
% % Create textbox
annotation(fig,'textbox',...
    [0.50820192694941,0.797502808513678,0.091623034308718,0.063265304845207],...
    'Color',[1 1 1],...
    'String',{'ECoG'},...
    'FontSize',12,...
    'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.821811394405027,0.797502808513678,0.091623034308718,0.063265304845207],...
    'Color',[1 1 1],...
    'String',{'ECoG'},...
    'FontSize',12,...
    'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.746143003190934,0.221850014278788,0.170118338566207,0.045123725435446],...
    'Color',[0 0 0],...
    'String',{'Hipp - ECoG'},...
    'FontSize',12,...
    'FitBoxToText','on', ...
    'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.060650887573965,0.952420670402097,0.164201178979239,0.042212517370824],...
    'String',{'f = [60 80] Hz'},...
    'FontSize',11,...
    'FitBoxToText','on', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.060650887573965,0.719524018291455,0.164201178979239,0.042212517370824],...
    'String',{'f = [11 14] Hz'},...
    'FontSize',11,...
    'FitBoxToText','on', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.060650887573965,0.483716158029441,0.164201178979239,0.042212517370824],...
    'String',{'f = [4 8] Hz'},...
    'FontSize',11,...
    'FitBoxToText','on', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.757396449704136,0.483716158029441,0.140532540631365,0.042212517370824],...
    'String',{'f = [4 8] Hz'},...
    'FontSize',11,...
    'FitBoxToText','on', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.060650887573965,0.246452693691741,0.140532540631365,0.042212517370824],...
    'String',{'f = [4 8] Hz'},...
    'FontSize',11,...
    'FitBoxToText','on', ...
    'EdgeColor','none');
%% Save figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42_r2_1.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Granger_Fig\Revision\Main Granger_Fig_S42_r2_1.tiff');


%%
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Granger_Fig\Revision\');
print(fig,'F:\Vasileios\Task Analysis\Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Granger_Fig\Revision\Main Granger_Fig_S42_r2_1','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Granger_Fig\Revision\Main Granger_Fig_S42_r2_1.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Granger_Fig\Revision\Main Granger_Fig_S42_r2_1.tiff');
