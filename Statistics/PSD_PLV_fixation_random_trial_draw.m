%% Calculate PSD for different number of trials - zscore
iSS_current = 5; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]
k=1;
nIterations = 20;

% for iTrial = 2:2:length(dataBipolar_SS{iSS_current+1}.trial)
for iRand = 1:nIterations
    trialpool = 1:length(dataBipolar_SS{iSS_current}.trial);
    nTrials_toSelect = round(0.05*length(trialpool));
    trialsSelection = randi([trialpool(1) trialpool(end)],1,nTrials_toSelect)
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
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = [1:100];
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_fix{iRand} = ft_freqanalysis(cfg,Fix_Data{iSS})
    fr_enc{iRand} = ft_freqanalysis(cfg,Enc_Data{iSS})
    fr_maint{iRand} = ft_freqanalysis(cfg,Maint_Data{iSS})
    
    
end

for iSS = iSS_current% 1:nSet_Size
    cfg = [];
    cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
    cfg.trials = trialsSelection;
    Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

%% Scrambled labels
cfg = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';
cfg.output      = 'pow';
cfg.foi         = [1:100];
cfg.tapsmofrq   = 2;
cfg.channel = 'all';
fr_maint_real = ft_freqanalysis(cfg,Maint_Data{iSS});
fr_fix_real  = ft_freqanalysis(cfg,Fix_Data{iSS});
nIterations = 200;

% for iTrial = 2:2:length(dataBipolar_SS{iSS_current+1}.trial)
for iRand = 1:nIterations
nLabels = length(dataBipolar.label);
randlabels = [randi([1 nLabels],1,round(0.3*nLabels));randi([1 nLabels],1,round(0.3*nLabels))];

for j = 1:size(randlabels,2)
    fr_maint{iRand} = fr_maint_real;
    fr_fix{iRand} = fr_fix_real;
    if i>2
    fr_maint{iRand}.powspctrm(randlabels(1,j),:) = fr_maint_real.powspctrm(randlabels(2,j),:);
    fr_fix{iRand}.powspctrm(randlabels(1,j),:) = fr_fix_real.powspctrm(randlabels(2,j),:);
    end
end

end
%%
figure;
for i = 1:nIterations
    semilogx(fr_fix{i}.freq, fr_fix{i}.powspctrm(18,:),'k');%C2
    hold on;
end
hold off;
figure;
for i = 1:nIterations
    semilogx(fr_enc{i}.freq, fr_enc{i}.powspctrm(18,:),'g'); %C2
    hold on;
end
%%
FreqBand = [11,14];
freqAxis = fr_maint{iRand}.freq;
[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));


PSD_In_Band_SS = [];
iSS_to_Plot = 4;
for iRand = 1:nIterations
    for iCh = 9:72
%         PSD_temp = abs(fr_enc{iRand}.powspctrm(iCh,indFreq1:indFreq2)-fr_fix{iRand}.powspctrm(iCh,indFreq1:indFreq2))./fr_fix{iRand}.powspctrm(iCh,indFreq1:indFreq2);
        PSD_temp2 = abs(fr_maint{iRand}.powspctrm(iCh,indFreq1:indFreq2)-fr_fix{iRand}.powspctrm(iCh,indFreq1:indFreq2))./fr_fix{iRand}.powspctrm(iCh,indFreq1:indFreq2);
        
        PSD_In_Band_SS{iRand}(iCh-8) = median(zscore(PSD_temp));
        PSD_vector(iCh-8,iRand) = mean(PSD_temp2);
    end
    PSD_In_Band_ToPlot = reshape(PSD_In_Band_SS{iRand},8,8);
    %     figure;
    %     imagesc(PSD_In_Band_ToPlot,[0 1])
    
    ax1 = gca;
    ax1.XTick = 1:8;
    ax1.YTick = 1:8;
    ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
    ax1.YTickLabel = 1:8;
    ax1.YDir = 'normal'
    
    
    colorbar
    colormap(bluewhitered)
    set(gca,'Fontsize',16)
end


%% PLV
iSS_current = 4; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]
k=1;
nIterations = 10;

for iRand = 1:nIterations
    trialpool = 1:length(dataBipolar_SS{iSS_current}.trial);
    nTrials_toSelect = round(0.3*length(trialpool));
    trialsSelection = randi([trialpool(1) trialpool(end)],1,nTrials_toSelect);
    nSet_Size = size(dataBipolar_SS,2);
    Maint_Data = {};
    for iSS = iSS_current% 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    for iPair = 257:320%1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(iPair,1);
        nChannel_2 = nChannelPairs(iPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,5);
        for iSS = iSS_current
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,Maint_Data{iSS});%ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
        end
        
        %% Phase coherence
        for iSS = iSS_current
            
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.foi         = [1:29 30:5:100];
            
            
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1:2;
            Freq            = ft_freqanalysis(cfg,data_SinglePair_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'none';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            freqAxis = [1:29 30:5:100];
            FreqBand = [16 29];
            [~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
            [~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
            
            PLV_spectrum = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))'
            PLV_vector(iPair,iRand) = mean(abs(PLV_spectrum(indFreq1:indFreq2)));
        end
    end
end


