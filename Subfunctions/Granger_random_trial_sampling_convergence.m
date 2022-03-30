%% Here we select randomly the same number of correct trials as the incorrect
%% We iterate 50 times where this random selection is performed
%% We calculate the granger for each task condition and for each iteration
%% We compute the median of Granger across the iterations


%% Calculate Granger for different number of trials - convergence
iSS_current = 5; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]
k=1;
nIterations = 200;
hipp_chanPairs = nChannelPairs%hipp_chanPairs;
iPair = 82
% for iTrial = 2:2:length(dataBipolar_SS{iSS_current+1}.trial)
for iRand = 1:nIterations
    trialpool = 1:length(dataBipolar_SS{iSS_current}.trial);
    nTrials_toSelect = 21%size(TrialInformationTable.Macro,1)-size(TrialInformationTable_iEEG_clean,1);
    trialsSelection = randi([trialpool(1) trialpool(end)],1,nTrials_toSelect);
    nPid = 42; %patient ID
    nSet_Size = size(dataBipolar_SS,2);
    Maint_Data = {};
    for iSS = iSS_current% 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    nSet_Size = size(dataBipolar_SS,2);
    Enc_Data = {};
    for iSS = iSS_current%1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Enc_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    nSet_Size = size(dataBipolar_SS,2);
    Fix_Data = {};
    for iSS = iSS_current%1:nSet_Size
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Fix_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    freq                   = [];
    freq.freqcfg           = [];
    freq.freqcfg.method    = 'mtmfft';
    freq.freqcfg.foi       = 4:1:100;
    freq.freqcfg.output    = 'fourier';
    freq.freqcfg.tapsmofrq = 2;
    Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS_current});
    Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS_current});
    Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS_current});
    
    
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    grngChan1 = strChannelNameList{hipp_chanPairs(iPair,1)};
    grngChan2 = strChannelNameList{hipp_chanPairs(iPair,2)};
    grangercfg.channelcmb = {grngChan1,grngChan2};
    grangercfg.granger.conditional = 'no';
    grangercfg.granger.sfmethod = 'bivariate';
    
    gdata = [];
    gdata.Maint    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata.Enc      = ft_connectivityanalysis(grangercfg, Encoding_freq);
    gdata.Fix      = ft_connectivityanalysis(grangercfg, Fixation_freq);
    
    % close all;
    cfg           = [];
    cfg.parameter = 'grangerspctrm';
    cfg.zlim      = [0 0.2];
    % figure;ft_connectivityplot(cfg, gdata.Maint);
    % figure;ft_connectivityplot(cfg, gdata.Enc);
    Granger_trial_cortex_hipp_enc(k,:) = gdata.Enc.grangerspctrm(1,:);
    Granger_trial_hipp_cortex_enc(k,:) =  gdata.Enc.grangerspctrm(2,:);
    Granger_trial_hipp_cortex_maint(k,:) = gdata.Maint.grangerspctrm(1,:);
    Granger_trial_cortex_hipp_maint(k,:) = gdata.Maint.grangerspctrm(2,:)
    k = k+1;
    
end
%%
% Visualize for maintenance
figure;
for iRand = 1:nIterations
    light_blue = [0.30,0.75,0.93];
    light_red       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',light_red};
    freq_ax    = gdata.Fix.freq;
    semilogx(gdata.Maint.freq,Granger_trial_hipp_cortex_maint(iRand,:)*100,'Color','r','LineWidth',3);
    hold on;
    xlim([4 30])
    set(gca,'FontSize',16);
    xlabel('Frequency (Hz)');
    ylabel('Granger')
    set(gca,'box', 'off')
end

for j = 1:size(Granger_trial_hipp_cortex_maint,2)
    medianGranger_randTrials_maint_hc(j) = median(Granger_trial_hipp_cortex_maint(:,j));
end
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_hc*100,'k','LineWidth',3);
xlim([4 30])

for j = 1:size(Granger_trial_cortex_hipp_maint,2)
    medianGranger_randTrials_maint_ch(j) = median(Granger_trial_cortex_hipp_maint(:,j));
end

% Visualize for encoding
figure;
for iRand = 1:nIterations
    light_blue = [0.30,0.75,0.93];
    light_red       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',light_red};
    freq_ax    = gdata.Fix.freq;
    semilogx(gdata.Enc.freq,Granger_trial_cortex_hipp_enc(iRand,:)*100,'Color','b','LineWidth',3);
    hold on;
    xlim([4 30])
    set(gca,'FontSize',16);
    xlabel('Frequency (Hz)');
    ylabel('Granger')
    set(gca,'box', 'off')
    k = k+1;
end

for j = 1:size(Granger_trial_cortex_hipp_enc,2)
    medianGranger_randTrials_enc_ch(j) = median(Granger_trial_cortex_hipp_enc(:,j));
end
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_ch*100,'k','LineWidth',3);
xlim([4 30])


for j = 1:size(Granger_trial_hipp_cortex_enc,2)
    medianGranger_randTrials_enc_hc(j) = median(Granger_trial_hipp_cortex_enc(:,j));
end

figure;
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_ch*100,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_hc*100,'Color',Colors{2},'LineWidth',3);
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_hc*100,'Color',Colors{3},'LineWidth',3);
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_ch*100,'Color',Colors{4},'LineWidth',3);

xlim([4 30])



%%
id = 15;
Granger_random_trial_sampling_scalp{id} =Granger_random_trial_sampling_scalp{id-1};
Granger_random_trial_sampling_scalp{id}.pID = 44
Granger_random_trial_sampling_scalp{id}.Median_gdata.Maint{2} = medianGranger_randTrials_maint_hc*100;
Granger_random_trial_sampling_scalp{id}.Median_gdata.Maint{1} = medianGranger_randTrials_maint_ch*100;
Granger_random_trial_sampling_scalp{id}.Median_gdata.Enc{1} = medianGranger_randTrials_enc_hc*100;
Granger_random_trial_sampling_scalp{id}.Median_gdata.Enc{2} = medianGranger_randTrials_enc_ch*100;
