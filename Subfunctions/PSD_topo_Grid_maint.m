
FreqBand = [11,14];
freqAxis = [4:100];
[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));


%time points
time_enc = 10;
time_maint = 19;
if strcmp(TaskPeriod,'encod')
    time = time_enc 
elseif strcmp(TaskPeriod,'maint')
    time = time_maint;
    
end
time_fix = 3;

%% Plot
PSD_In_Band_SS = [];
iSS_to_Plot = 4; 

for iCh = 9:72
    
    PSD_temp = TFR_baselined.powspctrm(iCh,indFreq1:indFreq2,time)                                           -0.7;%abs(fr_maint.powspctrm(iCh,indFreq1:indFreq2)-fr_fix.powspctrm(iCh,indFreq1:indFreq2))./fr_fix.powspctrm(iCh,indFreq1:indFreq2);
    PSD_In_Band_SS(iCh-8) = mean(PSD_temp);
end

%% Plot the heatmap of the PLV for a certain bipolar channel and for a certain set size

PSD_In_Band_ToPlot = reshape(PSD_In_Band_SS,8,8);
figure;
imagesc(PSD_In_Band_ToPlot,[0 1])

ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'


colorbar
colormap(bluewhitered)
set(gca,'Fontsize',16)