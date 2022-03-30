% load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\SubjGranger_spectra2.mat')
% load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Final_SubjGranger_spectra.mat');
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Final_SubjGranger_spectra_recalculatedGranger.mat')


SubjGranger_spectra = Granger_spectra;


% load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Granger for 3 sec maintenance\Granger for 3 sec maintenance.mat');
% SubjGranger_spectra = Granger_maint_3s;


dark_Green = [0.02 0.47 0.04];
orange       = [1 0.45 0.45];
Colors     = {dark_Green,'g','r',orange};
close all;
figure('units','normalized','outerposition',[0 0 1 1])
ha = tight_subplot(4,4,[.07 .07],[.1 .04],[.07 .03])
for i =1:length(SubjGranger_spectra)
    CortexHipp_Enc = SubjGranger_spectra{i}.Enc.grangerspctrm(1,:);
    HippCortexHipp_Enc = SubjGranger_spectra{i}.Enc.grangerspctrm(2,:);
    HippCortex_Maint = SubjGranger_spectra{i}.Maint.grangerspctrm(1,:);
    CortexHipp_Maint = SubjGranger_spectra{i}.Maint.grangerspctrm(2,:);
    freq_ax = [4:100];%SubjGranger_spectra{i}.Enc.freq;%
    theta_freq = freq_ax(1:5);
    %% Percentiles
%     clear Prcntile_enc Prcntile_maint
    if i == 3 || i == 4|| i == 6||i == 8 || i == 9 
        percentMaint = 80;
        percentEnc = percentMaint;
    elseif i == 7
        percentMaint = 65;
        percentEnc = percentMaint

%     elseif i ==11
%         percent = 93
    elseif  i ==12 || i ==13
        percentMaint = 80;
        percentEnc = percentMaint;
    elseif i ==11
        percentEnc = 70;
        percentMaint = 90;
    elseif i ==10
          percentEnc = 75;
        percentMaint = 60;  
    else
        percentEnc = 95;
        percentMaint = percentEnc;
    end
% if i == 4 || i == 6 || i == 9 
%     percent = 99;
% elseif i ==8  || i == 5
%     percent = 95;
% elseif i ==12|| i == 10
%         percent = 90;
% 
% else
% percent = 98;
% end
    
    Prcntile_enc{i,95} = repelem(prctile(CortexHipp_Enc-HippCortexHipp_Enc,percentEnc),size(CortexHipp_Enc,2));
    Prcntile_maint{i,95} = repelem(prctile(HippCortex_Maint-CortexHipp_Maint,percentMaint),size(CortexHipp_Enc,2));
    
%     
%     Prcntile_enc{i,95} = repelem(prctile(HippCortexHipp_Enc-CortexHipp_Enc,percent),size(CortexHipp_Enc,2));
%     Prcntile_maint{i,95} = repelem(prctile(CortexHipp_Maint-HippCortex_Maint,percent),size(CortexHipp_Enc,2));
    %% Plots
%     figure;
    axes(ha(i))
    semilogx(freq_ax, SubjGranger_spectra{i}.Enc.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',4)
    hold on;
    semilogx(freq_ax, SubjGranger_spectra{i}.Enc.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',4)
    semilogx(freq_ax, SubjGranger_spectra{i}.Maint.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',4)
    semilogx(freq_ax, SubjGranger_spectra{i}.Maint.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',4)
    indMaxBandCortexHipp_Enc{i} = Difference_Bar(CortexHipp_Enc-HippCortexHipp_Enc,Prcntile_enc{i,95},...
        freq_ax,[0 0.2], max(max(CortexHipp_Enc),max(HippCortex_Maint))+0.01,Colors{1});
    indMaxBandHippCortex_Maint{i} = Difference_Bar(HippCortex_Maint-CortexHipp_Maint,Prcntile_maint{i,95},...
        freq_ax,[0 0.2], max(max(CortexHipp_Enc),max(HippCortex_Maint))+0.02,Colors{3});
    
%      indMaxBandCortexHipp_Enc{i} = Difference_Bar(HippCortexHipp_Enc-CortexHipp_Enc,Prcntile_enc{i,95},...
%         freq_ax,[0 0.2], max(max(CortexHipp_Enc),max(HippCortex_Maint))+0.01,Colors{1});
%     indMaxBandHippCortex_Maint{i} = Difference_Bar(CortexHipp_Maint-HippCortex_Maint,Prcntile_maint{i,95},...
%         freq_ax,[0 0.2], max(max(CortexHipp_Enc),max(HippCortex_Maint))+0.02,Colors{3});
% 

    
    xlim([4 30])
    upperLim = max(max(CortexHipp_Enc),max(HippCortex_Maint));
    set(gca, 'FontSize',16,'box','off');
    ylabel('Granger')
    xlabel('Frequency (Hz)')
    
    [p_enc(i),h_enc(i)] = ranksum(HippCortexHipp_Enc,CortexHipp_Enc);
    [p_maint(i),h_maint(i)] = ranksum(HippCortex_Maint,CortexHipp_Maint);
    
end