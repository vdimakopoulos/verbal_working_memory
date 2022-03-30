strPath = 'F:/Vasileios/Task Analysis/Data/Analysis Data/Figure Data/210122 Granger_spectra/Granger_spectra_Fig.mat';
strPath_incorrect = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\210122 Granger_spectra\Delta_Granger_Incorrect_Trials_EEG.mat'
flag_incorrect = 1;
[data_maint, data_enc] = ft_prepare_GrangerData_permutation_test(strPath_incorrect,flag_incorrect);

cfg = [];
cfg.event_comparisons = {[1 2]}
stat_maint = statistics_fr(cfg,data_maint);
stat_enc = statistics_fr(cfg,data_enc);

sgnf_bands_maint =stat_maint{1}.hipp_cortex_maint.freq(stat_maint{1}.hipp_cortex_maint.mask)
sgnf_bands_enc = stat_enc{1}.hipp_cortex_maint.freq(stat_enc{1}.hipp_cortex_maint.mask);
clc;
if length(sgnf_bands_maint) == 0
    disp('No significance frequency bands for maintenance')
else
    for i = 1:length(sgnf_bands_maint)
        disp(sprintf('Significant frequency bins maintenance: %d Hz',sgnf_bands_maint(i)))
    end
end
if length(sgnf_bands_enc) == 0
    disp('No significance frequency bands for encoding')
else
    for j = 1:length(sgnf_bands_enc)
        disp(sprintf('Significant frequency bins encoding: %d Hz',sgnf_bands_enc(j)))
    end
end

freq  = stat_maint{1}.hipp_cortex_maint.freq;
t_maint = stat_maint{1}.hipp_cortex_maint.stat;
foi_ind = [1:27]; %corresponds to a foi of [4 30]

if flag_incorrect
    
    t_crit_val = max(t_maint(stat_maint{1}.hipp_cortex_maint.mask))
else
    t_crit_val = min(t_maint(stat_maint{1}.hipp_cortex_maint.mask));
end
% t_crit_val = 2.33;
figure('Position',[680,379,377,512]);
set(gcf,'Color','white')
subplot(211);
semilogx(freq,t_maint,'r*','MarkerSize',6);
hold on;
semilogx(freq(stat_maint{1}.hipp_cortex_maint.mask),t_maint(stat_maint{1}.hipp_cortex_maint.mask),'ko','MarkerSize',  10)
if all(stat_maint{1}.hipp_cortex_maint.mask==0)
    disp('No crtical t-value to plot')
else
    semilogx(freq,ones(length(freq),1)*t_crit_val,'k--')
end
xlim([4 30]);
xlabel('Frequency (Hz)');
ylabel('T-value')
set(gca,'FontSize',16,'box','off','tickdir','Out')
set(gca,'XTick',[4 10 20 30],'XTickLabel',[4 10 20 30])


freq  = stat_enc{1}.hipp_cortex_maint.freq;
t_enc = stat_enc{1}.hipp_cortex_maint.stat;
foi_ind = [1:27];
if flag_incorrect
    t_crit_val = max(t_enc(stat_enc{1}.hipp_cortex_maint.mask))
    t_enc(t_enc>2) =t_enc(t_enc>2)-2;
    t_enc(t_enc<-2) =t_enc(t_enc<-2)+2
else
    t_crit_val = min(t_maint(stat_enc{1}.hipp_cortex_maint.mask));
end
% t_crit_val = 2.24;
subplot(212);
semilogx(freq,t_enc,'b*','MarkerSize',6);
hold on;
semilogx(freq(stat_enc{1}.hipp_cortex_maint.mask),t_enc(stat_enc{1}.hipp_cortex_maint.mask),'ko','MarkerSize',  10)
if all(stat_enc{1}.hipp_cortex_maint.mask==0)
    disp('No crtical t-value to plot')
     

else
    
    semilogx(freq,ones(length(freq),1)*t_crit_val,'k--')

end
xlim([4 30])
xlabel('Frequency (Hz)');
ylabel('T-value')
set(gca,'FontSize',16,'box','off','tickdir','Out')
set(gca,'XTick',[4 10 20 30],'XTickLabel',[4 10 20 30])

