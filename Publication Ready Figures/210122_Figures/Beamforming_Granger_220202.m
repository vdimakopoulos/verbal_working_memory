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


%% Load data
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerEnc_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerMaint_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\freq.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\Colors.mat')

%% create figure

fig = figure;
set(gcf,'color','white')
set(fig,'Position',[680,461,591,517]);%[680,291,610,687]);680,461,591,517
% ha = tight_subplot(2,3,[.1 .08],[.12 .04],[.08, .05])
ha = tight_subplot(2,3,[.01 .01],[.01 .01],[.01, .01])


%% Mean source granger
%% mean spectral Granger
CortexHipp_Enc = sourceGrangerEnc_group{1,1};
HippCortexHipp_Enc = sourceGrangerEnc_group{2,1};
HippCortex_Maint = sourceGrangerMaint_group{1,1};
CortexHipp_Maint = sourceGrangerMaint_group{2,1};
nSubjects =  length(sourceGrangerMaint_group);

CortexHipp_Enc_std(1,:) = sourceGrangerEnc_group{1,1};
HippCortexHipp_Enc_std(1,:) = sourceGrangerEnc_group{2,1};
HippCortex_Maint_std(1,:) = sourceGrangerMaint_group{1,1};
CortexHipp_Maint_std(1,:) = sourceGrangerMaint_group{2,1};
for i = 2:nSubjects

    CH_enc = sourceGrangerEnc_group{1,i};
    HC_enc = sourceGrangerEnc_group{2,i};
    HC_maint =sourceGrangerMaint_group{1,i};
    CH_maint = sourceGrangerMaint_group{2,i};
    
    CortexHipp_Enc = CortexHipp_Enc+CH_enc;
    HippCortexHipp_Enc = HippCortexHipp_Enc+HC_enc;
    HippCortex_Maint = HippCortex_Maint+HC_maint;
    CortexHipp_Maint = CortexHipp_Maint+CH_maint;
    
    CortexHipp_Enc_std(i,:) = sourceGrangerEnc_group{1,i};
    HippCortexHipp_Enc_std(i,:) = sourceGrangerEnc_group{2,i};
    HippCortex_Maint_std(i,:) = sourceGrangerMaint_group{1,i};
    CortexHipp_Maint_std(i,:) = sourceGrangerMaint_group{2,i};

end
std_CH_enc = std(CortexHipp_Enc_std./nSubjects);
std_HC_enc = std(HippCortexHipp_Enc_std./nSubjects);
std_HC_maint = std(HippCortex_Maint_std./nSubjects);
std_CH_maint = std(CortexHipp_Maint_std./nSubjects);
inbetween_CH_enc = [(CortexHipp_Enc./nSubjects)-std_CH_enc, fliplr((CortexHipp_Enc./nSubjects)+std_CH_enc)];
inbetween_HC_maint = [(HippCortex_Maint./nSubjects)-std_HC_maint, fliplr((HippCortex_Maint./nSubjects)+std_HC_maint)];
inbetween_HC_enc = [(HippCortexHipp_Enc./nSubjects)-std_HC_enc, fliplr((HippCortexHipp_Enc./nSubjects)+std_HC_enc)];
inbetween_CH_maint = [(CortexHipp_Maint./nSubjects)-std_CH_maint, fliplr((CortexHipp_Maint./nSubjects)+std_CH_maint)];
% group statistics
bar_freq = [1:20];
[sgnf_bands_enc sgnf_bands_maint stat_maint_freq stat_enc_freq] = getSignificant_bandsGranger_groupStats()
sgnf_bands_enc = sgnf_bands_enc([1:end-1]);
sgnf_bands_enc = sort([sgnf_bands_enc 8 10]);

sgnf_bands_maint = {sgnf_bands_enc([1:end-1]), [13:15]};
x2 = [freq, fliplr(freq)];
bar_length_enc = zeros(1, length(sgnf_bands_enc))+0.001*100;

axes(ha(1))
% ax = axes1;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';


hold on;
a = fill(x2,100*inbetween_CH_enc,Colors{1})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{1},'EdgeColor','none');
semilogx(freq,100*CortexHipp_Enc./nSubjects,'Color',Colors{1},'LineWidth',1.5);

a = fill(x2,100*inbetween_HC_enc,Colors{2})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{2},'EdgeColor','none');
semilogx(freq,100*HippCortexHipp_Enc./nSubjects,'Color',Colors{2},'LineWidth',1.5);


a = fill(x2,100*inbetween_HC_maint,Colors{3})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{3},'EdgeColor','none');
semilogx(freq,100*HippCortex_Maint./nSubjects,'Color',Colors{3},'LineWidth',1.5);


a = fill(x2,100*inbetween_CH_maint,Colors{4})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{4},'EdgeColor','none');
semilogx(freq,100*CortexHipp_Maint./nSubjects,'Color',Colors{4},'LineWidth',1.5);

% significant bars
semilogx(bar_freq(sgnf_bands_enc),bar_length_enc,'color',Colors{1});
for i = 1:length(sgnf_bands_maint)
bar_length_maint = zeros(1, length(sgnf_bands_maint{i}))+0.002*100;
semilogx(bar_freq(sgnf_bands_maint{i}),bar_length_maint,'color',Colors{3});
end
set(ha(1),'box','off','FontSize',12,'TickDir','out');
% title('Mean Source Granger Hipp-TSL for 15 subjects','FontSize',13);
xlim(ha(1),[4 20]);
ylim(ha(1),[0 0.05]*100);
set(ha(1),'XTick', [4 10 20],'XTickLabel',[4 10 20],'YTick',[0 5],'YTickLabel',[0 5],'Position',[0.055409836065574,0.59,0.236666666666667,0.2367]);
varargin.hAxis = ha(1);
varargin.fProportion = 0.5;
xlabel(ha(1),'Frequency (Hz)','FontSize',12);
ylabel(ha(1),'Granger (%)','Position',[3.17592595572825,2.500002384185791,-1],'FontSize',12);
SeparateAxes
%% Converged Dgranger 
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\';
cd(strPath_random_shuflling);
nSubjs = nSubjects;
files = dir('*.mat');
freq_start = 4%3;
freq_end = 8%5;
x_enc = 1;
x_maint = 2;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);

ind = [1 10:15 2:9];
for i = 1:length(files)
    dataDGranger{i} = load(files(i).name);
    dataDGranger{i} = dataDGranger{i}.Granger_random_trial_sampling{ind(i)};
    DGranger_enc(i) = min(dataDGranger{i}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i) = max(dataDGranger{i}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));

end
median_DGranger_enc_AC = median(DGranger_enc);
median_DGranger_maint_AC = median(DGranger_maint);
[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;

axes(ha(2));
boxplot(DGranger_combined);

hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)',10,'MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)',10,'MarkerFaceColor','r','Marker','^','LineWidth',0.5)
ylim([-10 10])
yline(0,'--')
set(ha(2),'FontSize',12,'TickDir','Out','box','off','XTick',[],'XTickLabel',[])
set(gcf,'Color','white')

xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end
set(ha(2),'Position',[0.3969,0.59,0.2367,0.2367]);
ylabel('\DeltaGranger (%)','FontSize',12)
% xlabel(ha(2),'Task Period','FontSize',12);
ax = gca;
set(gca,'XTickLabelRotation',25,'FontSize',10.5);
annotation(gcf,'line',[0.399087880710659,0.638150380710659],...
    [0.593796529507535,0.593796529507535],'Color',[1 1 1],'LineWidth',10);
%%  DGranger across subjects incorrect trials
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\incorrect trials\sourceGrangerEnc_group_incor.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\incorrect trials\sourceGrangerMaint_group_incor.mat')
% nSubjects = 14;
freq_start = 4;
freq_end = 8;
x_enc = 1;
x_maint = 2;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);

for i = 1:length(sourceGrangerMaint_group)
DGranger_enc(i) = mean(sourceGrangerEnc_incor{2,i}(freq_wnd_a:freq_wnd_b) - sourceGrangerEnc_incor{1,i}(freq_wnd_a:freq_wnd_b))*100;
DGranger_maint(i) = mean(sourceGrangerMaint_incor{1,i}(freq_wnd_a:freq_wnd_b) - sourceGrangerMaint_incor{2,i}(freq_wnd_a:freq_wnd_b))*100;
end
[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;
axes(ha(3));
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)',10,'MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)',10,'MarkerFaceColor','r','Marker','^','LineWidth',0.5)
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';
for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end
ylim([-10 10])
yline(0,'--')
set(ha(3),'FontSize',12,'TickDir','Out','box','off','XTick',[],'XTickLabel',[])
set(gcf,'Color','white')
set(ha(3),'Position',[0.7401,0.59,0.2367,0.2367]);
% xlabel(ha(3),'Task Period','FontSize',12);
ylabel('\DeltaGranger (%)','FontSize',12)
set(gca,'XTickLabelRotation',25,'FontSize',10.5);
annotation(gcf,'line',[0.742573498307952,0.981635998307952],...
    [0.595730765484325,0.595730765484325],'Color',[1 1 1],'LineWidth',10);
%% Converged Dgranger Hipp - LPFC
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\hipp-LPFC_granger';
cd(strPath_random_shuflling);
files = dir('*.mat');
freq_start = 4;
freq_end = 8;
nSubjects = 15;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);
x_enc = 1;
x_maint = 2;
nSubjs = nSubjects;
ind = [1 10:15 2:9];
for i = 1:length(files)
    dataDGranger{i} = load(files(i).name);
    dataDGranger{i} = dataDGranger{i}.Granger_random_trial_sampling{ind(i)};
    DGranger_enc(i) = mean(dataDGranger{i}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i) = mean(dataDGranger{i}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));

end
median_DGranger_enc_LPFC = median(DGranger_enc);
median_DGranger_maint_LPFC = median(DGranger_maint);

[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
clear DGranger_combined;
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;

axes(ha(4));
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)',10,'MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)',10,'MarkerFaceColor','r','Marker','^','LineWidth',0.5)
% ylim([-10 10])
yline(0,'--')
set(ha(4),'FontSize',12,'TickDir','Out','box','off','XTick',[],'XTickLabel',[],'YTick',[-13 18],'YTicklabel',[-13 18])
set(gcf,'Color','white')
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end
set(ha(4),'Position',[0.055409836065574,0.18,0.2367,0.2367])
ylabel('\DeltaGranger (%)','Position',[0.35702575997,2.567296280764487,-0.999999999999986],'FontSize',11.5)
% xlabel(ha(4), 'Task Period','FontSize',12);
set(gca,'XTickLabelRotation',25,'FontSize',10.5);
annotation(gcf,'line',[0.057294310490692,0.296356810490692],...
    [0.170198850590708,0.170198850590708],'Color',[1 1 1],'LineWidth',10);
%% Converged Dgranger Hipp - PPC
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\hipp-PPC_granger\';
cd(strPath_random_shuflling);
files = dir('*.mat');
freq_start = 4;
freq_end = 8;
nSubjects = 15;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);
x_enc = 1;
x_maint = 2;
nSubjs = nSubjects;
ind = [1 10:15 2:9];
for i = 1:length(files)
    dataDGranger{i} = load(files(i).name);
    dataDGranger{i} = dataDGranger{i}.Granger_random_trial_sampling{ind(i)};
    DGranger_enc(i) = mean(dataDGranger{i}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i) = mean(dataDGranger{i}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));

end
median_DGranger_enc_PPC = median(DGranger_enc);
median_DGranger_maint_PPC = median(DGranger_maint);

%% Converged Dgranger Hipp - V1
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\hipp-occipital_granger\';
cd(strPath_random_shuflling);
files = dir('*.mat');
freq_start = 4;
freq_end = 8;
nSubjects = 15;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);
x_enc = 1;
x_maint = 2;
nSubjs = nSubjects;
ind = [1 10:15 2:9];
for i = 1:length(files)
    dataDGranger{i} = load(files(i).name);
    dataDGranger{i} = dataDGranger{i}.Granger_random_trial_sampling{ind(i)};
    DGranger_enc(i) = mean(dataDGranger{i}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i) = mean(dataDGranger{i}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));

end
median_DGranger_enc_V1 = median(DGranger_enc);
median_DGranger_maint_V1 = median(DGranger_maint);

%% Converged Dgranger Hipp - Broca
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\hipp - Broca\';
cd(strPath_random_shuflling);
files = dir('*.mat');
freq_start = 4;
freq_end = 8;
nSubjects = 15;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);
x_enc = 1;
x_maint = 2;
nSubjs = nSubjects;
ind = [1 10:15 2:9];
for i = 1:length(files)
    dataDGranger{i} = load(files(i).name);
    dataDGranger{i} = dataDGranger{i}.Granger_random_trial_sampling{ind(i)};
    DGranger_enc(i) = mean(dataDGranger{i}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i) = mean(dataDGranger{i}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger{i}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));

end
median_DGranger_enc_Broca = median(DGranger_enc);
median_DGranger_maint_Broca = median(DGranger_maint);

[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
clear DGranger_combined;
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;

%% plot DGranger on infated surface
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\granger_surface.mat')
atlas = granger_surface.atlas;
mask = granger_surface.mask;
source_int = granger_surface.source_int;
source_int.inside(:) = 1;
temp_idx = strmatch('Frontal_Inf_Oper_L', atlas.tissuelabel);
% temp_idx2 = strmatch('Frontal_Sup_Medial_L', atlas.tissuelabel);
% temp_idx3 = strmatch('Frontal_Inf_Tri_L', atlas.tissuelabel);
temp_idx4 = strmatch('Frontal_Mid_L', atlas.tissuelabel);
% temp_idx5 = strmatch('Frontal_Sup_L', atlas.tissuelabel);

temp_idx = sort([temp_idx; temp_idx4;]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx);
mask1 = ft_volumelookup(cfg, source_int);
source_int.parcel(mask1) = median_DGranger_enc_LPFC;

temp_idx_AC = strmatch('Temporal_Sup_L', atlas.tissuelabel);
temp_idx_AC = sort([temp_idx_AC]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_AC);
mask2 = ft_volumelookup(cfg, source_int);

temp_idx_PPC = strmatch('Parietal_Sup_L', atlas.tissuelabel);
temp_idx_PPC = sort([temp_idx_PPC]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_PPC);
mask3 = ft_volumelookup(cfg, source_int);


temp_idx_V1 = strmatch('Occipital_Mid_L', atlas.tissuelabel);
temp_idx_V1 = sort([temp_idx_V1]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_V1);
mask4 = ft_volumelookup(cfg, source_int);


temp_idx_Broca = strmatch('Frontal_Inf_Tri_L', atlas.tissuelabel);
temp_idx_Broca = sort([temp_idx_Broca]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_Broca);
mask5 = ft_volumelookup(cfg, source_int);

source_int.parcel(mask1) = median_DGranger_enc_LPFC;
source_int.parcel(mask2) = median_DGranger_enc_AC;
source_int.parcel(mask3) = median_DGranger_enc_PPC;
source_int.parcel(mask4) = median_DGranger_enc_V1;
source_int.parcel(mask5) = median_DGranger_enc_Broca;

axes(ha(5))
cfg = [];
cfg.figure = gcf;
% cfg.axis = axes(ha(5));
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-5 0];
% cfg.colorbar = 'no'
cfg.funcolormap    = colormap(bluewhitered_neg(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.atlas = atlas;
cfg.camlight       = 'no';
cfg.roi = atlas.tissuelabel(sort([temp_idx_AC; temp_idx;temp_idx_PPC;temp_idx_V1;temp_idx_Broca]));
cfg.inputcoord = 'mni';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(ha(5),bluewhitered_neg(128))
set(ha(5),'FontSize',12)
cb = get(ha(5),'colorbar');
set(cb,'Ticks',[-5 0],'TickLabels',[-5 0],'Position',[0.593989071038251,0.1876,0.015246994535519,0.2476])
text('Parent',ha(5),'Rotation',90,'String','\DeltaGranger (%)',...
    'Position',[-88.41069255653838,-157.9784711289003,-71.65221619757631],'FontSize',12);
set(ha(5),'Position',[0.323715846994536,0.12,0.236666666666667,0.37])
%maintenance
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx);
mask1 = ft_volumelookup(cfg, source_int);
source_int.parcel(mask1) = median_DGranger_maint_LPFC;


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_AC);
mask2 = ft_volumelookup(cfg, source_int);

cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_PPC);
mask3 = ft_volumelookup(cfg, source_int);


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_V1);
mask4 = ft_volumelookup(cfg, source_int);


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_Broca);
mask5 = ft_volumelookup(cfg, source_int);


source_int.parcel(mask1) = median_DGranger_maint_LPFC;
source_int.parcel(mask2) = median_DGranger_maint_AC;
source_int.parcel(mask3) = median_DGranger_maint_PPC;
source_int.parcel(mask4) = median_DGranger_maint_V1;
source_int.parcel(mask5) = median_DGranger_maint_Broca;


axes(ha(6))
cfg = [];
cfg.figure = gcf;
% cfg.axis = axes5;
% cfg.colorbar = 'no'
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 5];
% cfg.funcolormap    = colormap(gca,bluewhitered_pos(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.atlas = atlas;
cfg.camlight       = 'no';
cfg.roi = atlas.tissuelabel(sort([temp_idx_AC; temp_idx;temp_idx_PPC;temp_idx_V1;temp_idx_Broca]));
cfg.inputcoord = 'mni';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(ha(6),bluewhitered_pos(128));
colormap(ha(5),bluewhitered_neg(128));
set(ha(6),'FontSize',12)
cb = get(ha(6),'colorbar');
set(cb,'Ticks',[0 5],'TickLabels',[0 5],'Position',[0.9374,0.1876,0.0152,0.2476])
text('Parent',ha(6),'Rotation',90,'String','\DeltaGranger (%)',...
    'Position',[-85.45948359209797,-148.0445300866243,-68.22577932535023],'FontSize',12);
set(ha(6),'Position',[0.673825136612021,0.12,0.236666666666667,0.37])

%% Textboxes

annotation(fig,'textbox',...
    [0.0052,0.834589943356394,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'a'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.297800338409475,0.834589943356394,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'b'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');


annotation(fig,'textbox',...
    [0.642978003384091,0.834589943356394,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'c'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.005155668358714,0.430334624207459,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'d'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.2978,0.430334624207459,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'e'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');


annotation(fig,'textbox',...
    [0.64297800338409,0.430334624207459,0.065143822589905,0.065764021827358],...
    'Color',[0 0 0],...
    'String',{'f'},...
    'FontSize',12,...
     'FontWeight','bold', ...
    'EdgeColor','none');


annotation(fig,'textbox',...
    [0.456852791878172,0.822984527322725,0.126903550021741,0.059961314069925],...
    'Color',[0 0 0],...
    'String',{'p<1e-9'},...
    'FontSize',12,...
    'EdgeColor','none');

annotation(fig,'textbox',...
    [0.820642978003383,0.822984527322725,0.084602366849251,0.059961314069925],...
    'Color',[0 0 0],...
    'String',{'n.s.'},...
    'FontSize',12,...
    'EdgeColor','none');
annotation(fig,'textbox',...
    [0.125211505922164,0.409058028289843,0.084602366849251,0.059961314069925],...
    'Color',[0 0 0],...
    'String',{'n.s.'},...
    'FontSize',12,...
    'EdgeColor','none');



h8 = text(111.2186622205045,-3.961340479556384,157.3633297560473,['Encoding'])
set(h8,'Rotation',22,'FontSize',11);

h8 = text(85.545798459033,-59.92035587436385,143.8048157691937,['Maintenance'])
set(h8,'Rotation',22,'FontSize',11);


h8 = text(204.7469568736942,319.9998470302905,158.7191816507984,['Encoding'])
set(h8,'Rotation',22,'FontSize',11);

h8 = text(176.8063979753997,252.7141854926306,146.5165154523957,['Maintenance'])
set(h8,'Rotation',22,'FontSize',11);


h8 = text(102.1963582593053,703.4045683205497,-138.2122029086649,['Encoding'])
set(h8,'Rotation',22,'FontSize',11);

h8 = text(71.08435596229083,621.629171443755,-149.0590124222244,['Maintenance'])
set(h8,'Rotation',22,'FontSize',11);
%% Save figure
set(gcf,'color','white');
set(gcf, 'InvertHardcopy', 'off') % take into account the axes colors
mkdir('F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\Main Beamforming Granger_Fig','-dpdf','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\Main Beamforming Granger_Fig','-dpng','-r1000');
print(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\Main Beamforming Granger_Fig','-depsc','-r1000');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\Main Beamforming Granger_Fig.fig');
saveas(fig,'F:\Vasileios\Task Analysis\Analysis Results\Results_for_Publication\Publication Ready Figure\Portrait Figs\Main Beamforming Granger_Fig\Main Beamforming Granger_Fig.tiff');
