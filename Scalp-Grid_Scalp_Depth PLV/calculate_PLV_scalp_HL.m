function [PLV_Pairs_Scalp,nChannelPairs,freqAxis] = calculate_PLV_scalp_HL(macro_data,scalp_S1,scalp_S2,Scalp_TrialInfoTable,...
    montage,chans_to_be_bipolar,PLV_pair_flag,Task_Period_Flag)

PLV_scalp_ECOG_all_channels = 0;
IncorrectTrialFlag = 1;
%PLV_pair_flag: indicates the type of coupling between the recording modalities 
    %PLV_pair_flag = 1 ==> Scalp - AHL,PHL coupling
    %PLV_pair_flag = 2 ==> Scalp - ECoG Grid coupling

%% Merge sessions
cfg = [];
dataBipolar_scalp = ft_appenddata(cfg,scalp_S1,scalp_S2);

%% Reref based on the median of auricalar channels
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
% cfg               = []
% cfg.channel       = {'all', '-_Submm', '-_Submp'}
% cfg.reref         =  'yes' 
% cfg.refchannel    =   {'_A1' '_A2'}
% cfg.refmethod     =  'median'
% dataBipolar       = ft_preprocessing(cfg,dataBipolar_scalp);

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
macro_data = ft_resampledata(cfg,macro_data);
dataBipolar = ft_resampledata(cfg,dataBipolar_scalp);
% dataBipolar = ft_resampledata(cfg,dataBipolar);

%% Only for incorrect trials
if IncorrectTrialFlag
    IncorrectTrials = Scalp_TrialInfoTable(find(~Scalp_TrialInfoTable.Correct),:);
    cfg = [];
    cfg.trials = find(~Scalp_TrialInfoTable.Correct);
    dataBipolar = ft_selectdata(cfg,dataBipolar);
end
%% ICA
cfg = [];
cfg.method       = 'runica'
cfg.channel      = {'all' '-_Submm','-_Submp'}%,'-_A1','-_A2'};
if IncorrectTrialFlag
    cfg.trials = 'all';
else
    cfg.trials       = setdiff([1:77],[8 14 34]); % Reject trials: 8,14,34
end
cfg.numcomponent = 'all';
cfg.demean       = 'yes';
cfg.feedback     =  'text'
IC_components = ft_componentanalysis(cfg,dataBipolar)

%% Reject IC components and reconstruct
cfg = [];
cfg.component = [1 2 3 4 5 8 10 21 23];
cfg.demean = 'yes';
dataBipolar = ft_rejectcomponent(cfg,IC_components);

cfg = [];
cfg.channel = {'all'} % '-C3','-C4','-Cz','-O1','-P3','-T1','-T2'}
dataBipolarScalp_Reconstructed = ft_preprocessing(cfg,dataBipolar);
if ~IncorrectTrialFlag
    cfg = []
    cfg.trials = setdiff([1:77],[8 14 34]);
    macro_data = ft_preprocessing(cfg,macro_data);
end
%% Merge Scalp and macro data
    cfg = [];
    cfg.keepsampleinfo = 'no' 
    dataBipolar_M = ft_appenddata(cfg,dataBipolarScalp_Reconstructed,macro_data);



%% Select only correct trials
if IncorrectTrialFlag
    dataBipolar_SS = dataBipolar_M; %% For the incorrect trials
    TrialInformationTable = IncorrectTrials;
else
    [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar_M,Scalp_TrialInfoTable);
end
%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable);



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
        cfg.latency = [-6,-5-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end

end
%% Channel pairs to run
%To be configured for every combination of electrode pairs you have

% clear nChannelPairs;

nGridChans = 64;
% nScalpChans = size(scalp_S1.label,1)-2;
nScalpChans = size(dataBipolar.label,1);

AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));

if PLV_pair_flag == 1
    nChannelPairs = [];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+nGridChans+AHL_chans+PHL_chans+nScalpChans+zeros(nScalpChans,1),(1:nScalpChans)']];
    end
else 
    nChannelPairs = [];

    if ~PLV_scalp_ECOG_all_channels
        Pz_loc =find(contains(dataBipolar_SS{1}.label, '_Pz')==1);
        P3_loc =find(contains(dataBipolar_SS{1}.label, '_P3')==1);
        Fz_loc =find(contains(dataBipolar_SS{1}.label, '_Fz')==1);
        scalp_chans_grid_pair = [Fz_loc P3_loc Pz_loc]%[Fz_loc Pz_loc]
        %      scalp_chans_grid_pair = [Pz_loc P3_loc]
    else
%         Pz_loc =find(contains(dataBipolar_SS{1}.label, '_Pz')==1);
%         P3_loc =find(contains(dataBipolar_SS{1}.label, '_P3')==1);
%         Fz_loc =find(contains(dataBipolar_SS{1}.label, '_Fz')==1);
%         scalp_chans_grid_pair = [Fz_loc P3_loc Pz_loc]
        scalp_chans_grid_pair = [1:size(dataBipolarScalp_Reconstructed.label)];
        
    end

    for i=1:length(scalp_chans_grid_pair)
        nChannelPairs = [nChannelPairs;[scalp_chans_grid_pair(i)+zeros(nGridChans,1),nScalpChans+AHL_chans+(1:nGridChans)']];
    end
    
end

%% Loop over pairs
% PLV_Pairs = {};
% Coh_Pairs = {};
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
        cfg.foi         = [1:1:30 30:5:100];%1:1:30;
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
end