function Connectivity_metrics = compute_Connectivity_metrics(taskPeriodData,strChannelNameList, nChannelPairs,foi,nSet_Size,metric)

if nargin<=5
    metric = {'psi','ppc','wpli_debiased'};
end
for iPair = 1:size(nChannelPairs,1)
    %% Channel numbers
    nChannel_1 = nChannelPairs(iPair,1);
    nChannel_2 = nChannelPairs(iPair,2);
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,5);
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS} = ft_preprocessing(cfg,taskPeriodData{iSS})
    end
    for iSS = 1:nSet_Size
        
        cfg = [];
        cfg.method      = 'mtmfft';
        cfg.output      = 'fourier';
        cfg.foi         = foi;
        cfg.tapsmofrq   = 2;
        cfg.channel = 1:2;
        Freq =ft_freqanalysis(cfg,   data_SinglePair_SS{iSS});
        cfg.output = 'powandcsd';
        cfg.channelcmb = {strChannelNameList(nChannel_1) strChannelNameList(nChannel_1)};
        Freq_csd   = ft_freqanalysis(cfg,    data_SinglePair_SS{iSS});
        %% PSI
        cfg           = [];
        cfg.method    =  metric{1};
        cfg.bandwidth = 5;
        psi{iSS}           = ft_connectivityanalysis(cfg, Freq);
        psi2{iSS}         = ft_connectivityanalysis(cfg, Freq_csd);
        
        if size(metric,2)>=2
            %% PPC
            cfg           = [];
            cfg.method    =  metric{2};
            ppc{iSS}           = ft_connectivityanalysis(cfg, Freq);
            if size(metric,2)>=3
                %% PLI
                cfg           = [];
                cfg.method    =  metric{3};
                pli{iSS}           = ft_connectivityanalysis(cfg, Freq);
            end
        end
    end
end
switch size(metric,2)
    case 1
        Connectivity_metrics.psi = psi;
        Connectivity_metrics.psi_csd  = psi2; 
    case 2
        Connectivity_metrics.psi = psi;
        Connectivity_metrics.psi_csd  = psi2;
        Connectivity_metrics.ppc = ppc;
        
    case 3
        Connectivity_metrics.psi = psi;
        Connectivity_metrics.psi_csd  = psi2;
        Connectivity_metrics.ppc = ppc;
        Connectivity_metrics.pli = pli;
end
