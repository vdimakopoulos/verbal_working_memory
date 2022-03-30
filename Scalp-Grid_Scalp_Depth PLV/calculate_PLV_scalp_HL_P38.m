function [PLV_Pairs_Scalp,nChannelPairs,freqAxis,dataBipolar_Scalp] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_Session,Scalp_TrialInfoTable,...
    montage,chans_to_be_bipolar,PLV_pair_flag,Task_Period_Flag)

%% Merge sessions
cfg = [];
dataBipolarScalp = data_Scalp_Session(1).data;
for nSes = 2:length(data_Scalp_Session)
    cfg = [];
    dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_Session(nSes).data);
end

Trial_Information_Table = Scalp_TrialInfoTable;

%% Resample
cfg = [];
cfg.resamplefs = 500;
macro_data = ft_resampledata(cfg,macro_data);
dataBipolar = ft_resampledata(cfg,dataBipolarScalp);

cfg = [];
cfg.channel = {'all' '-Subm1','-Subm2'};
dataBipolar = ft_selectdata(cfg,dataBipolar)


%% macro data
% Select only correct trials
[macro_SS,Trial_Information_Table] = Get_Only_Correct_Trials_FieldTrip(macro_data,Trial_Information_Table);

cfg = [];
cfg.trials = setdiff([1:length(macro_SS.trialinfo)],[177 178])
macro_SS = ft_selectdata(cfg,macro_SS);
%Update trial information table
Trial_Information_Table = Trial_Information_Table([1:176,179:size(Trial_Information_Table,1)],:);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[macro_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(macro_SS,Trial_Information_Table);

%% Reject trial
cfg =[];
cfg.trials = setdiff([1:length(macro_SS{4}.trialinfo)],[1 2:2:8 9 11:1:22 25 27 28 31 32 34:1:36 39:1:42 47 52 54:1:57 59 62 65 66 69:1:71 73 75 79 85 88 102 107 108 115 116 121 125 131 138 143 144])';
macro_SS{4} =  ft_selectdata(cfg,macro_SS{4});


%% EEG Preproc

%% Select only correct trials
[dataBipolar_Scalp_SS,Scalp_TrialInfoTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,Scalp_TrialInfoTable);

cfg = [];
cfg.trials = setdiff([1:length(dataBipolar_Scalp_SS.trialinfo)],[177 178])
dataBipolar_Scalp_SS = ft_selectdata(cfg,dataBipolar_Scalp_SS);
%Update trial information table
Scalp_TrialInfoTable = Scalp_TrialInfoTable([1:176,179:size(Scalp_TrialInfoTable,1)],:)



%% ICA
cfg = [];
cfg.method       = 'runica'
cfg.channel      = {'all' '-_Submm','-_Submp'};
cfg.numcomponent = 'all';
cfg.demean       = 'yes';
cfg.feedback     =  'text'
IC_components = ft_componentanalysis(cfg,dataBipolar_Scalp_SS)


%% Reject IC componentss
cfg = [];
cfg.component = [1:7 10 12 13 17 20];
cfg.demean = 'yes';
dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);


%% re-ref EEG
cfg               = []
cfg.channel       = {'all'}
cfg.reref         =  'yes'
cfg.refchannel    =   [1, 2]
cfg.refmethod     =  'avg'
dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp_Reconstructed);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5

Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_Scalp_SS,TrialInformationTableScalp_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolarScalp,Scalp_TrialInfoTable);

%% Reject trials
cfg =[];
cfg.trials = setdiff([1:length(dataBipolar_Scalp_SS{4}.trialinfo)],[1 2:2:8 9 11:1:22 25 27 28 31 32 34:1:36 39:1:42 47 52 54:1:57 59 62 65 66 69:1:71 73 75 79 85 88 102 107 108 115 116 121 125 131 138 143 144])';
dataBipolar_Scalp_SS{4} =  ft_selectdata(cfg,dataBipolar_Scalp_SS{4});


%% Merge Scalp and macro data
for i = 1:5
    cfg = [];
    cfg.keepsampleinfo = 'no'
    dataBipolar_SS{i} = ft_appenddata(cfg,dataBipolar_Scalp_SS{i},macro_SS{i});
    
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
        cfg.latency = [-6,-5-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end

    
end

%% Channel pairs to run
%To be configured for every combination of electrode pairs you have

% clear nChannelPairs;

nStripChans = length(find(contains(macro_data.label,'T')));
% nScalpChans = size(scalp_S1.label,1)-2;
nScalpChans = size(dataBipolarScalp.label,1);

PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
AHR_chans = length(find(contains(montage.labelold, 'AHR')==1));
PHR_chans = length(find(contains(montage.labelold, 'PHR')==1));
Depth_electrodes = +PHL_chans+AHR_chans+PHR_chans;

if PLV_pair_flag == 1
    nChannelPairs = [];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+nStripChans+Depth_electrodes+nScalpChans+zeros(nScalpChans,1),(1:nScalpChans)']];
    end

else
    nChannelPairs = [];

  
        Pz_loc =find(contains(dataBipolar_SS{1}.label, 'Pz')==1);
        P3_loc =find(contains(dataBipolar_SS{1}.label, 'P3')==1);
        Fz_loc =find(contains(dataBipolar_SS{1}.label, 'Fz')==1);
        scalp_chans_grid_pair = [1:21]%[Fz_loc P3_loc Pz_loc]
   

    for i=1:length(scalp_chans_grid_pair)
        nChannelPairs = [nChannelPairs;[scalp_chans_grid_pair(i)+zeros(nStripChans,1),nScalpChans+Depth_electrodes+(1:nStripChans)']];
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
        cfg.foi         = [1:1:30 30:5:100];
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

end

    

