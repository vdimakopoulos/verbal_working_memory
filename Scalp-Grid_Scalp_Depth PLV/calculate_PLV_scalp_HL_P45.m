function [PLV_Pairs_Scalp,nChannelPairs,freqAxis,dataBipolar_Scalp] = calculate_PLV_scalp_HL_P45(dataBipolar_SS,dataBipolarScalp_SS,...
    montage,chans_to_be_bipolar,PLV_pair_flag,Task_Period_Flag)
%% Merge scalp with macro data
cfg = [];
for i = 1:size(dataBipolar_SS,2)
   dataBipolar_SS{i} = ft_appenddata(cfg,dataBipolarScalp_SS{i},dataBipolar_SS{i});
end
%% Extract 2 seconds of encoding
if Task_Period_Flag == 1
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
%% Extract 2 seconds of maintenance
elseif Task_Period_Flag == 2
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end

else
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end

    
end

%% Channel pairs to run
%To be configured for every combination of electrode pairs you have

% clear nChannelPairs;

nStripChansGAR = length(find(contains(dataBipolar_SS{4}.label,'GAR')));
nStripChansGSR = length(find(contains(dataBipolar_SS{4}.label,'GAR')));
nStripChans = nStripChansGAR+nStripChansGSR;
ind_start =find(contains(dataBipolar_SS{4}.label,'GAR1'));

nScalpChans = size(dataBipolarScalp_SS{1}.label,1);

Depth_electrodes = length(montage.labelold);

if PLV_pair_flag == 1
    nChannelPairs = [];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+Depth_electrodes+nScalpChans+zeros(nScalpChans,1),(1:nScalpChans)']];
    end

else
        nChannelPairs = [];
        Pz_loc =find(contains(dataBipolar_SS{1}.label, 'Pz')==1);
        P3_loc =find(contains(dataBipolar_SS{1}.label, 'P3')==1);
        Fz_loc =find(contains(dataBipolar_SS{1}.label, 'Fz')==1);
        scalp_chans_grid_pair = [Fz_loc P3_loc Pz_loc]
        

    for i=1:length(scalp_chans_grid_pair)
        nChannelPairs = [nChannelPairs;[scalp_chans_grid_pair(i)+zeros(nStripChans,1),([ind_start:ind_start+nStripChans-1])']];
    end
end

%% Loop over pairs
tStart = tic;
for iPair = 1:size(nChannelPairs,1)
    %% Channel numbers
    nChannel_1 = nChannelPairs(iPair,1);
    nChannel_2 = nChannelPairs(iPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,5);
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
    end
    
    %% Phase coherence
    %     tic
    clear PLV_SS Coh_SS
    for iSS = 1:nSet_Size
        %% Frequency analysis
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
        cfg.foi         = [1:1:29 30:5:100];
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
        
        %% Coherence
        cfg             = [];
        cfg.method      = 'coh';
        cfg.complex     = 'complex';
        cfg.feedback    = 'none';
        Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
        
    end
    %     toc
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}.label';
    freqAxis = PLV_SS{iSS}.freq;
    for iSS = 1:nSet_Size
        PLV_Spectrum_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
        Coh_Spectrum_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
    end
    TimeInterval = [-2,-1/dataBipolar_SS{iSS}.fsample];
    
    %%
    PLV_Pairs(iPair,:) = PLV_Spectrum_SS;
    Coh_Pairs(iPair,:) = Coh_Spectrum_SS;
    
    tStop = toc(tStart);
    fprintf('\n\n\n\n\n\n\n\n\n\n Pair %d %f seconds elapsed \n',iPair,tStop)
end


PLV_Pairs_Scalp = PLV_Pairs;
dataBipolar_Scalp = dataBipolar_SS;