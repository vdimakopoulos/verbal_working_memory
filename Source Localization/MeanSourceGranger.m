%% load data
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerEnc_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerMaint_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\freq.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\Colors.mat')
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

figure;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


hold on;
a = fill(x2,100*inbetween_CH_enc,Colors{1})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{1},'EdgeColor','none');
semilogx(freq,100*CortexHipp_Enc./nSubjects,'Color',Colors{1},'LineWidth',3);

a = fill(x2,100*inbetween_HC_enc,Colors{2})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{2},'EdgeColor','none');
semilogx(freq,100*HippCortexHipp_Enc./nSubjects,'Color',Colors{2},'LineWidth',2);


a = fill(x2,100*inbetween_HC_maint,Colors{3})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{3},'EdgeColor','none');
semilogx(freq,100*HippCortex_Maint./nSubjects,'Color',Colors{3},'LineWidth',3);


a = fill(x2,100*inbetween_CH_maint,Colors{4})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{4},'EdgeColor','none');
semilogx(freq,100*CortexHipp_Maint./nSubjects,'Color',Colors{4},'LineWidth',2);

% significant bars
semilogx(bar_freq(sgnf_bands_enc),bar_length_enc,'color',Colors{1});
for i = 1:length(sgnf_bands_maint)
bar_length_maint = zeros(1, length(sgnf_bands_maint{i}))+0.002*100;
semilogx(bar_freq(sgnf_bands_maint{i}),bar_length_maint,'color',Colors{3});
end
ylabel('Granger');
xlabel('Frequency (Hz)');
set(gca,'box','off','FontSize',16);
% title('Mean Source Granger Hipp-TSL for 15 subjects','FontSize',13);
xlim([4 20]);
ylim([0 0.05]*100);
set(gca,'XTick', [4 10 20],'XTickLabel',[4 10 20]);
set(gcf,'Color','white');
%% mean Granger incorrect trials
%% load data
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\incorrect trials\sourceGrangerEnc_group_incor.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\incorrect trials\sourceGrangerMaint_group_incor.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\freq.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\Colors.mat')
sourceGrangerEnc_group = sourceGrangerEnc_incor;
sourceGrangerMaint_group = sourceGrangerMaint_incor;

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
bar_length_enc = zeros(1, length(sgnf_bands_enc))+0.001;
bar_length_maint = zeros(1, length(sgnf_bands_maint))+0.002;

x2 = [freq, fliplr(freq)];
figure;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


hold on;
a = fill(x2,inbetween_CH_enc,Colors{1})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{1},'EdgeColor','none');
semilogx(freq,CortexHipp_Enc./nSubjects,'Color',Colors{1},'LineWidth',3);

a = fill(x2,inbetween_HC_enc,Colors{2})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{2},'EdgeColor','none');
semilogx(freq,HippCortexHipp_Enc./nSubjects,'Color',Colors{2},'LineWidth',2);


a = fill(x2,inbetween_HC_maint,Colors{3})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{3},'EdgeColor','none');
semilogx(freq,HippCortex_Maint./nSubjects,'Color',Colors{3},'LineWidth',3);


a = fill(x2,inbetween_CH_maint,Colors{4})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{4},'EdgeColor','none');
semilogx(freq,CortexHipp_Maint./nSubjects,'Color',Colors{4},'LineWidth',2);

% significant bars
% semilogx(bar_freq(sgnf_bands_enc),bar_length_enc,'color',Colors{1});
% semilogx(bar_freq(sgnf_bands_maint),bar_length_maint,'color',Colors{3});

ylabel('Granger');
xlabel('Frequency (Hz)');
set(gca,'box','off','FontSize',16);
title('Mean Source Granger Hipp-TSL for 15 subjects','FontSize',13);
xlim([4 20]);
% ylim([0 0.05]);
set(gca,'XTick', [4 10 20],'XTickLabel',[4 10 20]);

%% median Granger

for j = 1:length(freq)
    tempMedian_CH_enc{j} = sourceGrangerEnc_group{1,1}(j);
    tempMedian_HC_enc{j} = sourceGrangerEnc_group{2,1}(j);
    tempMedian_HC_maint{j} = sourceGrangerMaint_group{1,1}(j);
    tempMedian_CH_maint{j} = sourceGrangerMaint_group{2,1}(j);
    for iSub =2:length(sourceGrangerEnc_group)
        tempMedian_CH_enc{j} = [tempMedian_CH_enc{j}  sourceGrangerEnc_group{1,i}(j)];
        tempMedian_HC_enc{j} = [tempMedian_HC_enc{j}  sourceGrangerEnc_group{2,i}(j)];
        tempMedian_HC_maint{j} = [tempMedian_HC_maint{j} sourceGrangerMaint_group{1,i}(j)];
        tempMedian_CH_maint{j} = [tempMedian_CH_maint{j} sourceGrangerMaint_group{2,1}(j)];
    end
end

for j = 1:length(freq)
    median_CH_enc(j) = median(tempMedian_CH_enc{j});
    median_HC_enc(j) = median(tempMedian_HC_enc{j});
    median_HC_maint(j) = median(tempMedian_HC_maint{j});
    median_CH_maint(j) = median(tempMedian_CH_maint{j});
    
end


figure;
hold on;
semilogx(freq,median_CH_enc,'Color',Colors{1},'LineWidth',3);
semilogx(freq,median_HC_enc,'Color',Colors{2},'LineWidth',3);
semilogx(freq,median_HC_maint,'Color',Colors{3},'LineWidth',3);
semilogx(freq,median_CH_maint,'Color',Colors{4},'LineWidth',3);
ylabel('Granger');
xlabel('Frequency (Hz)');
set(gca,'box','off','FontSize',16);
title('Median Source Granger Hipp-TSL for 15 subjects','FontSize',13);
xlim([4 20]);
% ylim([0 0.05]);
set(gca,'XTick', [4 10 20],'XTickLabel',[4 10 20]);


%% mean TFR Granger
strPath_enc = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\TFR_granger\encoding_new';
cd(strPath_enc)
files = dir('*.mat')
 max_length_enc= 0;
for i = 1:length(files)
    temp_TFR_enc{i} = load(files(i).name);
    temp_TFR_enc{i} = temp_TFR_enc{i}.Granger_SS_enc{4};
    
    %% find max length
    if size(temp_TFR_enc{i},2)>max_length_enc
        max_length_enc = size(temp_TFR_enc{i},2);
    end
  
end
%interpolation
for i = 1:length(temp_TFR_enc)
  if size(temp_TFR_enc{i},2)<max_length_enc
        temp_TFR_enc{i} = [temp_TFR_enc{i} zeros(size(temp_TFR_enc{i},1),max_length_enc-size(temp_TFR_enc{i},2))];
  end
end

strPath_maint = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\TFR_granger\maintenance_new';
cd(strPath_maint)
files = dir('*.mat')
 max_length_maint= 0;

for i = 1:length(files)
    temp_TFR_maint{i} = load(files(i).name);
    temp_TFR_maint{i} = temp_TFR_maint{i}.Granger_SS_maint{4};
    
        %% find max length
    if size(temp_TFR_maint{i},2)>max_length_maint
        max_length_maint = size(temp_TFR_maint{i},2);
    end

end

%interpolation
for i = 1:length(temp_TFR_maint)
    if size(temp_TFR_maint{i},2)<max_length_maint
        temp_TFR_maint{i} = [temp_TFR_maint{i} zeros(size(temp_TFR_maint{i},1),max_length_maint-size(temp_TFR_maint{i},2))];
    end
end
sum_TFR_enc = temp_TFR_enc{1};
sum_TFR_maint = temp_TFR_maint{1};
for i = 2:length(temp_TFR_enc)
    sum_TFR_enc = sum_TFR_enc+temp_TFR_enc{i};
    sum_TFR_maint = sum_TFR_maint+temp_TFR_maint{i};
    
end

mean_TFR_enc =sum_TFR_enc./length(temp_TFR_enc);
mean_TFR_maint = sum_TFR_maint/length(temp_TFR_maint);

timeAxis_enc = [-5:2/max_length_enc:-3-2/max_length_enc];
timeAxis_maint = [-2:2/max_length_maint:0-2/max_length_maint];
freq_Axis = [1:20];

figure;

clim = [-0.05 0.05]*100%[-0.1 0.15];%[-0.05 0.05];
contourf(timeAxis_enc,freq_Axis,mean_TFR_enc,100,'LineColor','none');

set(gca,'clim',clim,'yscale','linear');
set(gca,'ytick',[4,5 10 20] ,'YTickLabel',[4,5,10,20]) %background color of grid

colormap(bluewhitered(128));
set(gca,'XTick',[-5 -3], 'XTickLabel',[-5 -3])
% ylim([4 100]);
ylim([4 20])
colorbar;

figure;
clim = [-0.05 0.05]*100%[-0.1 0.15];%[-0.05 0.05];
contourf(timeAxis_maint,freq_Axis,mean_TFR_maint,100,'LineColor','none');

set(gca,'clim',clim,'yscale','linear');
set(gca,'ytick',[4,5 10 20] ,'YTickLabel',[4,5,10,20]) %background color of grid

colormap(bluewhitered(128));
set(gca,'XTick',[-5 -3], 'XTickLabel',[-5 -3])
% ylim([4 100]);
ylim([4 20])
colorbar;


%% mean statistical source power
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceStatistics_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceStatistics_group_maint.mat')
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
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 10];
cfg.funcolormap    = 'jet';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, mean_stat_enc);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 10];
cfg.funcolormap    = 'jet';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, mean_stat_maint);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

%% DGranger across subjects
nSubjects = 15;
nSubjs =nSubjects;
x_enc = 1;
x_maint = 2;
freq_start = 4;
freq_end =8;

freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);

for i = 1:length(sourceGrangerMaint_group)
DGranger_enc(i) = min(sourceGrangerEnc_group{2,i}(freq_wnd_a:freq_wnd_b) - sourceGrangerEnc_group{1,i}(freq_wnd_a:freq_wnd_b))*100;
DGranger_maint(i) = max(sourceGrangerMaint_group{1,i}(freq_wnd_a:freq_wnd_b) - sourceGrangerMaint_group{2,i}(freq_wnd_a:freq_wnd_b))*100;
end
[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;
figure;
boxplot(DGranger_combined);

hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
% ylim([-20 20])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
hold on;
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end

%% DGranger across subjects incorrect trials
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
figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end
ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
%% Converged Dgranger 
strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\';
cd(strPath_random_shuflling);
nSubjs = 15;
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

figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')

xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end

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

figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
% ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end

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

[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
clear DGranger_combined;
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;

figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end
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

[EncDiff_sorted,Ienc] =  sort(DGranger_enc,'descend');
[MaintDiff_sorted,Imaint] =  sort(DGranger_maint,'descend');
clear DGranger_combined;
DGranger_combined(:,1) = EncDiff_sorted;
DGranger_combined(:,2) = MaintDiff_sorted;

figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
% ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end

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

figure;
boxplot(DGranger_combined);
hold on;
a = scatter(ones(1,length(DGranger_enc))',   DGranger_combined(:,1)','MarkerFaceColor','b','Marker','v','LineWidth',0.5)
b = scatter(ones(1,length(DGranger_enc))+1',   DGranger_combined(:,2)','MarkerFaceColor','r','Marker','^','LineWidth',0.5)
% ylim([-10 10])
ylabel('\DeltaGranger (%)')
yline(0,'--')
set(gca,'FontSize',16,'TickDir','Out','box','off','XTick',[1 2],'XTickLabel',{'Encoding' 'Maintenance'})
set(gcf,'Color','white')
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';

for i=1:nSubjects
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','Color',[0.65,0.65,0.65], 'LineStyle','-','LineWidth',1)
    hold on
end

%% Mean Granger Hipp -LPFC
CortexHipp_Enc = dataDGranger{1}.Median_gdata.Enc{2};
HippCortexHipp_Enc = dataDGranger{1}.Median_gdata.Enc{1};
HippCortex_Maint = dataDGranger{1}.Median_gdata.Maint{2};
CortexHipp_Maint = dataDGranger{1}.Median_gdata.Maint{1};
nSubjects =  length(dataDGranger);

CortexHipp_Enc_std(1,:) = dataDGranger{1}.Median_gdata.Enc{2};
HippCortexHipp_Enc_std(1,:) = dataDGranger{1}.Median_gdata.Enc{1};
HippCortex_Maint_std(1,:) = dataDGranger{1}.Median_gdata.Maint{2};
CortexHipp_Maint_std(1,:) = dataDGranger{1}.Median_gdata.Maint{1};
for i = 2:nSubjects
% if i ~=[13]
    CH_enc = dataDGranger{i}.Median_gdata.Enc{2};
    HC_enc = dataDGranger{i}.Median_gdata.Enc{1};
    HC_maint = dataDGranger{i}.Median_gdata.Maint{2};
    CH_maint = dataDGranger{i}.Median_gdata.Maint{1}
    
    CortexHipp_Enc = CortexHipp_Enc+CH_enc;
    HippCortexHipp_Enc = HippCortexHipp_Enc+HC_enc;
    HippCortex_Maint = HippCortex_Maint+HC_maint;
    CortexHipp_Maint = CortexHipp_Maint+CH_maint;
    
    CortexHipp_Enc_std(i,:) = dataDGranger{i}.Median_gdata.Enc{2};
    HippCortexHipp_Enc_std(i,:) = dataDGranger{i}.Median_gdata.Enc{1};
    HippCortex_Maint_std(i,:) = dataDGranger{i}.Median_gdata.Maint{2};
    CortexHipp_Maint_std(i,:) = dataDGranger{i}.Median_gdata.Maint{1};
% end
end
nSubjects = nSubjects-1;
std_CH_enc = std(CortexHipp_Enc_std./nSubjects);
std_HC_enc = std(HippCortexHipp_Enc_std./nSubjects);
std_HC_maint = std(HippCortex_Maint_std./nSubjects);
std_CH_maint = std(CortexHipp_Maint_std./nSubjects);
inbetween_CH_enc = [(CortexHipp_Enc./nSubjects)-std_CH_enc, fliplr((CortexHipp_Enc./nSubjects)+std_CH_enc)];
inbetween_HC_maint = [(HippCortex_Maint./nSubjects)-std_HC_maint, fliplr((HippCortex_Maint./nSubjects)+std_HC_maint)];
inbetween_HC_enc = [(HippCortexHipp_Enc./nSubjects)-std_HC_enc, fliplr((HippCortexHipp_Enc./nSubjects)+std_HC_enc)];
inbetween_CH_maint = [(CortexHipp_Maint./nSubjects)-std_CH_maint, fliplr((CortexHipp_Maint./nSubjects)+std_CH_maint)];
% group statistics
% bar_freq = [1:20];
% [sgnf_bands_enc sgnf_bands_maint stat_maint_freq stat_enc_freq] = getSignificant_bandsGranger_groupStats()
% bar_length_enc = zeros(1, length(sgnf_bands_enc))+0.001;
% bar_length_maint = zeros(1, length(sgnf_bands_maint))+0.002;

x2 = [freq, fliplr(freq)];
figure;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


hold on;
a = fill(x2,inbetween_CH_enc,Colors{1})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{1},'EdgeColor','none');
semilogx(freq,CortexHipp_Enc./nSubjects,'Color',Colors{1},'LineWidth',3);

a = fill(x2,inbetween_HC_enc,Colors{2})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{2},'EdgeColor','none');
semilogx(freq,HippCortexHipp_Enc./nSubjects,'Color',Colors{2},'LineWidth',2);


a = fill(x2,inbetween_HC_maint,Colors{3})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{3},'EdgeColor','none');
semilogx(freq,HippCortex_Maint./nSubjects,'Color',Colors{3},'LineWidth',3);


a = fill(x2,inbetween_CH_maint,Colors{4})
set(a,'FaceAlpha',0.3,'FaceColor',Colors{4},'EdgeColor','none');
semilogx(freq,CortexHipp_Maint./nSubjects,'Color',Colors{4},'LineWidth',2);

% significant bars
% semilogx(bar_freq(sgnf_bands_enc),bar_length_enc,'color',Colors{1});
% semilogx(bar_freq(sgnf_bands_maint),bar_length_maint,'color',Colors{3});

ylabel('Granger');
xlabel('Frequency (Hz)');
set(gca,'box','off','FontSize',16);
title('Mean Source Granger Hipp-LPFC for 15 subjects','FontSize',13);
xlim([4 20]);
% ylim([0 0.05]);
set(gca,'XTick', [4 10 20],'XTickLabel',[4 10 20]);


%% Convergence Plot Granger w.r.t. #trials
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\convergence_plot.mat')
x = convergence_plot.x;
y_enc = convergence_plot.y_enc;
y_maint = convergence_plot.y_maint;
error_enc = convergence_plot.error_enc*10;
error_maint = convergence_plot.error_maint*10;

figure;
% subplot(2,1,1)
hold on;
l1 =errorbar(x,y_enc,error_enc,'-','color','b');
a = fill([x(2:end) fliplr(x(2:end))],[y_enc(2:end)-error_enc(2:end) fliplr(y_enc(2:end)+error_enc(2:end))], 'b');
set(a,'FaceAlpha',0.3,'FaceColor','b','linestyle','-','EdgeColor','none','LineWidth',3)
xlim([1 40])
yline(0, '--','color','k','LineWidth',2)

xlabel('Number of trials')
ylabel('\DeltaGranger (%)')
set(gca,'Fontsize',14,'XTick',[0:5:40],'XTickLabel',[0:5:40],'TickDir','out')

hold on;
l2 = errorbar(x,y_maint,error_maint,'-','color','r')
a = fill([x(2:end) fliplr(x(2:end))],[y_maint(2:end)-error_maint(2:end) fliplr(y_maint(2:end)+error_maint(2:end))], 'r');
set(a,'FaceAlpha',0.3,'FaceColor','r','linestyle','-','EdgeColor','none','LineWidth',3)
xlim([1 40])
legend([l1 l2],{'Encoding','Maintenance'},'EdgeColor','none')
set(gca,'Fontsize',14,'XTick',[0:5:40],'XTickLabel',[0:5:40],'TickDir','out','YTick',[-50 -25 0 25 50],'YTickLabel',[-50 -25 0 25 50] )
xlabel('Number of trials')
ylabel('\DeltaGranger (%)')
ylim([-50 50])
%% Project DGranger onto brain surface
%% plot roi activity
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


cfg = [];
% cfg.figure = gcf;
cfg.axis = gca%ha(1);
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-5 0];
cfg.funcolormap    = 'bluewhitered_neg(128)';
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

set(gca,'FontSize',16)
cb = get(gca,'colorbar');
set(cb,'Ticks',[-5 0],'TickLabels',[-5 0])
text('Parent',gca,'Rotation',90,'String','\DeltaGranger (%)',...
    'Position',[-79.96295300909969,-193.821745203596,-34.41982876541988],'FontSize',16);

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


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 5];
cfg.funcolormap    = 'bluewhitered_pos(128)';
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

set(gca,'FontSize',16)
cb = get(gca,'colorbar');
set(cb,'Ticks',[0 5],'TickLabels',[0 5])
text('Parent',gca,'Rotation',90,'String','\DeltaGranger (%)',...
    'Position',[-79.96295300909969,-193.821745203596,-34.41982876541988],'FontSize',16);
