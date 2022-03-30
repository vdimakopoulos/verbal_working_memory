load('F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\210122 Granger_spectra\Granger_spectra_Fig.mat')
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Granger_spectra_incorrectTrials.mat')
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\PercentsCorr_IncorrGranger.mat')


SubjGranger_spectra_Correct = Granger_Spectra_Fig.spectra_recalculated;
SubjGranger_spectra_Incorrect = Granger_spectra_incorrectTrials;


dark_Green = [0.02 0.47 0.04];
orange       = [1 0.45 0.45];
Colors     = {dark_Green,'g','r',orange};
close all;
figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(4,4,[.07 .07],[.1 .04],[.07 .03])

for i =1:length(SubjGranger_spectra_Correct)

    HippCortex_MaintCorr = SubjGranger_spectra_Correct{i}.Maint.grangerspctrm(1,:);
    HippCortex_MaintIncorr = Granger_spectra_incorrectTrials{i}.Maint.grangerspctrm(1,:);
    if size(HippCortex_MaintCorr,2) == 100
        HippCortex_MaintCorr = HippCortex_MaintCorr(4:100);
    elseif size(HippCortex_MaintIncorr,2) == 100
        HippCortex_MaintIncorr = HippCortex_MaintIncorr(4:100);
    end
        
    freq_ax = [4:100];%SubjGranger_spectra{i}.Enc.freq;%
    percentMaint = Percents(i);
    percentMaintInc = Percents(i);
    Prcntile_maintCorr{i,95} = repelem(prctile(HippCortex_MaintCorr-HippCortex_MaintIncorr,percentMaint),size(HippCortex_MaintCorr,2));
    Prcntile_maintIncorr{i,95} = repelem(prctile(HippCortex_MaintCorr-HippCortex_MaintIncorr,percentMaintInc),size(HippCortex_MaintCorr,2));
    
    axes(ha(i));
    semilogx(freq_ax, HippCortex_MaintCorr,'Color','r','LineWidth',4);
    hold on;
    semilogx(freq_ax, HippCortex_MaintIncorr,'Color','k','LineWidth',4);
    
     indMaxBandMaintCorr_wrong{i} = Difference_Bar(HippCortex_MaintCorr-HippCortex_MaintIncorr,Prcntile_maintCorr{i,95},...
        freq_ax,[0 0.2], max(max(HippCortex_MaintCorr),max(HippCortex_MaintIncorr))+0.01,'r');
    indMaxBandMaintWrong_Corr{i} = Difference_Bar(HippCortex_MaintIncorr-HippCortex_MaintCorr,Prcntile_maintIncorr{i,95},...
        freq_ax,[0 0.2], max(max(HippCortex_MaintCorr),max(HippCortex_MaintIncorr))+0.02,'k');
    
    
    xlim([4 30])
    upperLim = max(max(HippCortex_MaintCorr),max(HippCortex_MaintIncorr));
    set(gca, 'FontSize',16,'box','off');
    ylabel('Granger')
    xlabel('Frequency (Hz)')

end