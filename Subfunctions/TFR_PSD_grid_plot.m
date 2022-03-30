figure('units','normalized','outerposition',[0 0 1 1])
set(gcf,'color','white')
ha = tight_subplot(8,8,[.05 .05],[.05 .05],[.05 .05])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';

FreqBand = [60,70];
freqAxis = [4:100];
[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));


%time points
time_enc = [13:14]%10;
time_maint = 19;
if strcmp(TaskPeriod,'encod')
    time = time_enc 
elseif strcmp(TaskPeriod,'maint')
    time = time_maint;
    
end

for i =9:72
Chan_to_plot = i%Bip_chans(1)
clim = [0 1];
subplot_flag = 0;
% axes(ha(plot_order(i-8)))

cfg = [];
cfg.baseline = [TFR2.time(1) TFR2.time(5)]
cfg.baselinetype = 'absolute';
TFR_baselined = ft_freqbaseline(cfg,TFR2)
PS = squeeze(TFR_baselined.powspctrm(Chan_to_plot,:,:));
if i == 26
    baseline = 1.5;
else
    baseline = 3;
end
TFR_psd_base = (squeeze(TFR_baselined.powspctrm(Chan_to_plot,:,:)))/baseline;%-0.5)/1.9%/3-1%
% contourf(TFR_psd.time,TFR_psd.freq,TFR_psd_base,100,'LineColor','none')

%Axes properties
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap 'jet' % to achieve correct colorbar for negative red peaks use flipud(jet) and -TFR_psd_base
colorbar
ylabel(strrep(dataBipolar.label{i},'_',''))

PSD_temp = TFR_psd_base(indFreq1:indFreq2,time);%TFR_baselined.powspctrm(Chan_to_plot,indFreq1:indFreq2,time)                                           -0.7;%abs(fr_maint.powspctrm(iCh,indFreq1:indFreq2)-fr_fix.powspctrm(iCh,indFreq1:indFreq2))./fr_fix.powspctrm(iCh,indFreq1:indFreq2);
PSD_In_Band_SS(Chan_to_plot-8) = mean(mean(PSD_temp));
end



PSD_In_Band_ToPlot = reshape(PSD_In_Band_SS,8,8);
figure;
imagesc(PSD_In_Band_ToPlot,[0.5 1])

ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'


colorbar
colormap(winter)
set(gca,'Fontsize',16)