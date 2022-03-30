function [Granger_struct, Granger_TFR_ss] = ...
    get_TFR_Granger(freq,nSubj,nChannelPairs,strChannelNameList,dataBipolar_SS)

Frequency = freq;
Logspace_freq = logspace(log10(Frequency(1)),log10(Frequency(end)),length(Frequency));
%% specify the pairs for each subject
switch nSubj
    case 1
        channel_cmb = {'AHL2-AHL3','T3'};
    case 2
        channel_cmb = {'PHL1-PHL3','C3'};
    case 3
        channel_cmb = {'AHR1-AHR5','C3'};
    case 4
        channel_cmb = {'ECL1-ECL3','T3'}%{'PHR1-PHR3','T4'};
    case 5
        channel_cmb = {'PHR1-PHR3','T4'};
    case 6
        channel_cmb = {'AHR1-AHR4','C4'};
    case 7
        channel_cmb = {'ECL1-ECL4','C3'}%{'AHL1-AHL4','C4'};
    case 8
        channel_cmb = {'AHR2-AHR5','C3'};
    case 9
        channel_cmb = {'AHL2-AHL4','C3'};
        
end

%% find the number of the pair
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))',channel_cmb{1})));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))',channel_cmb{2})));
nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);

for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,5);
    for iSS = 1:5
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_SS{iSS});%ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
    end
    
    %% Create dataset
    data_SinglePair_Rand_SS = data_SinglePair_SS;
    %% Granger tfr
    for iSS = 1:5
        
        cfg           = [];
        cfg.method    = 'mtmconvol';
        cfg.output    = 'fourier';
        cfg.taper     = 'hanning';
        cfg.foi       = 1:100;
        cfg.t_ftimwin = 0.15.*ones(size(4./cfg.foi',1), 1);
        cfg.tapsmofrq = 0.2*cfg.foi;
        cfg.toi       = -6:0.25:2;
        CrossFreq     = ft_freqanalysis(cfg, data_SinglePair_Rand_SS{iSS});
        
        grangercfg = [];
        grangercfg.method  = 'granger';
        grangercfg.channelcmb = channel_cmb;
     
        Granger_Rand_SS{iSS}= ft_connectivityanalysis(grangercfg,CrossFreq);
        CortexHipp_Enc = squeeze(Granger_Rand_SS{iSS}.grangerspctrm(2,1,:,:));
        HippCortex_Maint = (squeeze(Granger_Rand_SS{iSS}.grangerspctrm(1,2,:,:)));
        Granger_SS{iSS} = (HippCortex_Maint-CortexHipp_Enc);%CortexHipp_Enc - HippCortex_Maint;%
        
    end
    
end

Granger_TFR_ss = Granger_SS;
Granger_struct = Granger_Rand_SS;