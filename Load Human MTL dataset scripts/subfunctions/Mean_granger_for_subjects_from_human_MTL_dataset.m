load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Final_SubjGranger_spectra.mat');
% load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Granger_spectra_1_1_100_ss_68.mat');
SubjGranger_spectra = Granger_spectra;
% SubjGranger_spectra = Granger_spectra_1_1_100_ss_68;

dark_Green = [0.02 0.47 0.04];
orange       = [1 0.45 0.45];
Colors     = {dark_Green,'g','r',orange};
freq_ax    = [4:100];

CortexHipp_Enc = SubjGranger_spectra{1}.Enc.grangerspctrm(1,:);
HippCortexHipp_Enc = SubjGranger_spectra{1}.Enc.grangerspctrm(2,:);
HippCortex_Maint = SubjGranger_spectra{1}.Maint.grangerspctrm(1,:);
CortexHipp_Maint = SubjGranger_spectra{1}.Maint.grangerspctrm(2,:);
nSubjects =  length(SubjGranger_spectra);
for i = 2:nSubjects
    %     figure;
    %    semilogx(freq_ax, SubjGranger_spectra{i}.Enc.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3)
    %    hold on;
    %    semilogx(freq_ax, SubjGranger_spectra{ia}.Enc.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3)
    %    semilogx(freq_ax, SubjGranger_spectra{i}.Maint.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3)
    %    semilogx(freq_ax, SubjGranger_spectra{i}.Maint.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3)
    %    xlim([4 30])
    %
%     if i>=10
%         fr =[1:97]
%     else
%         fr= [4:100];
%     end
    CH_enc = SubjGranger_spectra{i}.Enc.grangerspctrm(1,:);
    HC_enc = SubjGranger_spectra{i}.Enc.grangerspctrm(2,:);
    HC_maint =SubjGranger_spectra{i}.Maint.grangerspctrm(1,:);
    CH_maint = SubjGranger_spectra{i}.Maint.grangerspctrm(2,:);
    
    %    Max_spectra_value(i) = max([max(CH_enc) max(HC_maint)])
    
    
    CortexHipp_Enc = CortexHipp_Enc+CH_enc;
    HippCortexHipp_Enc = HippCortexHipp_Enc+HC_enc;
    HippCortex_Maint = HippCortex_Maint+HC_maint;
    CortexHipp_Maint = CortexHipp_Maint+CH_maint;
end
% group statistics
[sgnf_bands_enc sgnf_bands_maint] = getSignificant_bandsGranger_groupStats()

figure;
% axes(ha(16))
semilogx(freq_ax,CortexHipp_Enc./nSubjects,'Color',Colors{1},'LineWidth',3)
hold on;
semilogx(freq_ax,HippCortexHipp_Enc./nSubjects,'Color',Colors{2},'LineWidth',3)
semilogx(freq_ax,HippCortex_Maint./nSubjects,'Color',Colors{3},'LineWidth',3)
semilogx(freq_ax,CortexHipp_Maint./nSubjects,'Color',Colors{4},'LineWidth',3)
ylabel('Granger');
xlabel('Frequency (Hz)')
set(gca,'box','off','FontSize',16)
title('Mean Granger between hipp and scalp for 15 subjects','FontSize',13)
xlim([4 30]);
ylim([0 0.05])
strSavePath = ['F:\Vasileios\Task Analysis\Analysis Results\Human MTL dataset results Granger\'];
strFilename = [strSavePath,'Mean Granger between hipp and scalp for 13 subjects'];
saveas(gcf,strFilename)
saveas(gcf,[strFilename,'.png'])

%% Median Granger spectra
for j = 1:length(freq_ax)
    tempMedian_CH_enc{j} = SubjGranger_spectra{1}.Enc.grangerspctrm(1,j);
    tempMedian_HC_enc{j} = SubjGranger_spectra{1}.Enc.grangerspctrm(2,j);
    tempMedian_HC_maint{j} = SubjGranger_spectra{1}.Maint.grangerspctrm(1,j);
    tempMedian_CH_maint{j} = SubjGranger_spectra{1}.Maint.grangerspctrm(2,j);
    for iSub =2:length(SubjGranger_spectra)
        tempMedian_CH_enc{j} = [tempMedian_CH_enc{j} SubjGranger_spectra{iSub}.Enc.grangerspctrm(1,j)];
        tempMedian_HC_enc{j} = [tempMedian_HC_enc{j} SubjGranger_spectra{iSub}.Enc.grangerspctrm(2,j)];
        tempMedian_HC_maint{j} = [tempMedian_HC_maint{j} SubjGranger_spectra{iSub}.Maint.grangerspctrm(1,j)];
        tempMedian_CH_maint{j} = [tempMedian_CH_maint{j} SubjGranger_spectra{iSub}.Maint.grangerspctrm(2,j)];
    end
end

for j = 1:length(freq_ax)
    median_CH_enc(j) = median(tempMedian_CH_enc{j});
    median_HC_enc(j) = median(tempMedian_HC_enc{j});
    median_HC_maint(j) = median(tempMedian_HC_maint{j});
    median_CH_maint(j) = median(tempMedian_CH_maint{j});
    
end

figure;
% axes(ha(16))
semilogx(freq_ax,median_CH_enc,'Color',Colors{1},'LineWidth',3)
hold on;
semilogx(freq_ax,median_HC_enc,'Color',Colors{2},'LineWidth',3)
semilogx(freq_ax,median_HC_maint,'Color',Colors{3},'LineWidth',3)
semilogx(freq_ax,median_CH_maint,'Color',Colors{4},'LineWidth',3)
ylabel('Granger');
xlabel('Frequency (Hz)')
set(gca,'box','off','FontSize',16)
xlim([4 30])

%% Mean Granger TFR
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\SubjGranger_TFR2')
MeanGrangerTFR = SubjGranger_TFR{1};
nSubjects = length(SubjGranger_TFR);
for i = 2:length(SubjGranger_TFR)
    MeanGrangerTFR = MeanGrangerTFR+SubjGranger_TFR{i};
end
% Plot the mean Granger TFR
grangerTimeAxis = [-6:0.25:2];
grangerFreqAxis = [1:100];
figure;
clim = [-0.02 0.02];%[-0.05 0.05];
% contourf(grangerTimeAxis,grangerFreqAxis,SubjGranger_TFR{i},100,'LineColor','none');

contourf(grangerTimeAxis,grangerFreqAxis,MeanGrangerTFR./nSubjects,100,'LineColor','none');

set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
colormap(bluewhitered(128));
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
ylim([4 100]);
ylim([4 30])
colorbar;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
set(gca,'FontSize',16,'box','off');
set(gcf,'Color','white')


%% Median TFR
for k =1:size(SubjGranger_TFR{1},1)
    for j = 1:size(SubjGranger_TFR{1},2)
        tempMedianTFR{k,j} = SubjGranger_TFR{1}(k,j);
       for iSub =2:length(SubjGranger_spectra)
            tempMedianTFR{k,j} = [tempMedianTFR{k,j} SubjGranger_TFR{iSub}(k,j)];
        end
    end
    
end

for k =1:size(SubjGranger_TFR{1},1)
    for j = 1:size(SubjGranger_TFR{1},2)
        MedianTFR_Granger(k,j) = median(tempMedianTFR{k,j});
               
    end
end


grangerTimeAxis = [-6:0.25:2];
grangerFreqAxis = [1:100];
figure;
clim = [-0.03 0.03];%[-0.05 0.05];
% contourf(grangerTimeAxis,grangerFreqAxis,SubjGranger_TFR{i},100,'LineColor','none');

contourf(grangerTimeAxis,grangerFreqAxis,MedianTFR_Granger,100,'LineColor','none');

set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
colormap(bluewhitered(128));
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
ylim([4 100]);
ylim([1 30])
colorbar;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
set(gca,'FontSize',16,'box','off');
set(gcf,'Color','white')
