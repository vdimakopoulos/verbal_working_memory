%% Paths
strPaths.Main = 'F:\Vasileios\';
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
strPaths.GeneralFunctions = 'F:\Vasileios\Task Analysis\Code\';
strPaths.Data = [strPaths.Main, 'Task Analysis\Data\'];
strPaths.ExtractedData = [strPaths.Main, 'Task Analysis\Extracted_Data Sternberg\'];
strPaths.Statistics = [strPaths.Main, 'Task Analysis\Code\Statistics\'];
strPaths.ChanLoc = [strPaths.Main, 'Task Analysis\Code\Channel Localization\'];
strPaths.EEGLAB_Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\EEGLAB Subfunctions\'];
strPaths.Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\']
% Results
strPaths.Results = [strPaths.Main,'Task Analysis\Analysis Results\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20191126\';
% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = 'F:\Vasileios\Toolboxes\eeglab14_1_1b\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(strPaths.Main)
addpath(strPaths.Project)
addpath(genpath(strPaths.GeneralFunctions))
addpath(strPaths.Data)
addpath(strPaths.ExtractedData)
addpath(strPaths.Statistics)
addpath(strPaths.ChanLoc)
addpath(strPaths.Subfunctions)
addpath(strPaths.EEGLAB_Subfunctions)
addpath(strPaths.Results)
addpath(strPaths.Toolboxes.FieldTrip)
% Remove EEGLAB from path

rmpath(genpath(strPaths.Toolboxes.EEGLAB))

Figures_Path = 'F:\Vasileios\Task Analysis\Analysis Results\Statistics\PLV_significance\'

strPlotColors = {'b','g','r','c','k'};
ft_defaults

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Load data for 2 sessions
Session1                   =      'F:\Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_01_Part_01';
Session2                   =      'F:\Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_02_Part_01';
data_s1                    =      load(Session1);
data_s2                    =      load(Session2);
TrialInformationTable      =      [data_s1.TrialInformationTable;data_s2.TrialInformationTable] ;
data_s1                    =      data_s1.dataMacro;
data_s2                    =      data_s2.dataMacro;

%% Load scalp data for 2 sessions
Sclp_Session1                    =      'F:\Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_01_Part_01';
Sclp_Session2                    =      'F:\Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_02_Part_01';
data_s1_scalp                    =      load(Sclp_Session1);
data_s2_scalp                    =      load(Sclp_Session2);
TrialInformationTable_Scalp      =      [data_s1_scalp.TrialInformationTable;data_s2_scalp.TrialInformationTable] ;
data_scalp_s1                    =      data_s1_scalp.dataScalp;
data_scalp_s2                    =      data_s2_scalp.dataScalp;

%% Flags for coupling types and other flags
Hippocampus_Grid_Coupling           =   1;
Hippocampus_Coupling                =   0;
Scalp_Hippocampus_coupling          =   0;
Scalp_ECoG_coupling                 =   0;
Task_Period = 2; %Values 1,2,3 for Encoding,Maintenance,Retrieval period correspondigly
PLV_significance = 0;
PLV_significance_maint_vs_encod = 1;
PLV_significance_maint_vs_fix = 1;
PSD_significance = 1;
SetSizeComparison = 1;
Gamma_Freq = 1;

if Hippocampus_Grid_Coupling
    coupling_type = 'Hippocampus_grid_coupling';
    Path_to_load = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\PLV_Pairs_grid_HC_all_setSizes_maintenance';
    Path_for_PLV_sgnf = [Figures_Path,'Hippocampus Grid Coupling\'];

elseif Hippocampus_Coupling
    coupling_type = 'Hippocampus_coupling';
    Path_to_load = 'F:\Vasileios\Task Analysis\Data\Analysis Data\AHL-PHL coupling\PLV_Pairs_AHL-PHL_maintenance';
    Path_for_PLV_sgnf = [Figures_Path,'Hippocampus Coupling\'];
elseif Scalp_Hippocampus_coupling
    coupling_type = 'Scalp_Hippocampus_coupling';
    Path_to_load = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\PLV_Pairs__scalp_AHL_PHL_all_setSizes_maintenance';
    Path_for_PLV_sgnf = [Figures_Path,'Scalp-Hippocampus Coupling\'];

elseif Scalp_ECoG_coupling
    coupling_type = 'Scalp_ECoG_coupling';
    Path_to_load = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\PLV_Pairs_scalp_grid_all_setSizes_maintenance';
    Path_for_PLV_sgnf = [Figures_Path,'Scalp-Grid Coupling\'];

end
if PLV_significance_maint_vs_encod
    %File Paths for the PLV Pairs in maintenance and encoding
    if Hippocampus_Grid_Coupling
        if Gamma_Freq
            PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_gamma.mat';
            PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_gamma.mat';
        else
            PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance.mat';
            PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding.mat';
        end
    elseif Hippocampus_Coupling
        PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\AHL-PHL coupling\PLV_Pairs_AHL-PHL_maintenance';
        PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\AHL-PHL coupling\PLV_Pairs_AHL-PHL_encoding';
    elseif Scalp_Hippocampus_coupling
        PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_maintenance.mat';
        PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_encoding.mat';
    elseif Scalp_ECoG_coupling
         PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_maintenance.mat';
        PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_encoding.mat';
    end
end
mkdir(Path_for_PLV_sgnf)

%% Merge sessions
cfg = [];
dataBipolar = ft_appenddata(cfg,data_s1,data_s2);
dataBipolarScalp = ft_appenddata(cfg,data_scalp_s1,data_scalp_s2);


[dataBipolar_test,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable);

%% Reref based on the average of the contacts in White matter
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
cfg               = []
cfg.reref         =  'yes' 
cfg.refchannel    =  [3,4] %white matter contacts referencing
cfg.refmethod     =  'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

%Reref scalp data based on the median of auricular channels
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
% cfg               =  []
% cfg.channel       =  {'all', '-_Submm', '-_Submp'}
% cfg.reref         =  'yes' 
% cfg.refchannel    =  {'_A1' '_A2'}
% cfg.refmethod     =  'median'
% dataBipolarSclp   =  ft_preprocessing(cfg,dataBipolarScalp);

%% Apply montage
clear montage;
montage.labelold = dataBipolar.label;
num_bipolar_chans = 6;
num_reference_chans = size(montage.labelold,1);
%prepare the montage matrix
montage_matrix = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar = [1 2 1 3 2 3 73 74 73 75 74 75 ]; %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end
%Append the bipolar channel labels to the reference channels labels
for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew = {dataBipolar.label{:},montage.labelnew{:}};
montage.tra = montage_matrix;
% dataBipolar = ft_apply_montage(dataBipolar,montage);

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod = 'bipolar'
cfg.refchannel = [82 85]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));


%% ICA
cfg = [];
cfg.method       = 'runica'
cfg.channel      = {'all' '-_Submm','-_Submp'}%,'-_A1','-_A2'};
cfg.trials       = setdiff([1:77],[8 14 34]); % Reject trials: 8,14,34
cfg.numcomponent = 'all';
cfg.demean       = 'yes';
cfg.feedback     =  'text'
IC_components = ft_componentanalysis(cfg,dataBipolarScalp)

%% Reject IC components and reconstruct
cfg = [];
cfg.component = [1 2 3 4 5 8 10 21 23];
cfg.demean = 'yes';
dataBipolarScalp = ft_rejectcomponent(cfg,IC_components);

cfg = [];
cfg.channel = {'all'} % '-C3','-C4','-Cz','-O1','-P3','-T1','-T2'}
dataBipolarSclp = ft_preprocessing(cfg,dataBipolarScalp);

cfg = []
cfg.trials = setdiff([1:77],[8 14 34]);
dataBipolar = ft_preprocessing(cfg,dataBipolar);

%% Merge Scalp and macro data
if Scalp_ECoG_coupling || Scalp_Hippocampus_coupling

    cfg = [];
    cfg.keepsampleinfo = 'no' 
    dataBipolar = ft_appenddata(cfg,dataBipolarSclp,dataBipolar);
end

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);

% %% Low Pass filtering - LPfreq = 125 Hz
% cfg = [];
% cfg.lpfilter = 'yes';
% cfg.lpfreq = 125;
% dataBipolar = ft_preprocessing(cfg,dataBipolar);

%% Select only correct trials
strChannelNameList = dataBipolar.label;
nChannelList = 1:length(strChannelNameList);

[dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable);

%% Extract the corresponding task period latency   
    if Task_Period == 1
        %% Extract 2 seconds of encoding latency
        nSet_Size = size(dataBipolar_SS,2);
        dataBipolar_Ret_SS = {};
        for iSS = 1:nSet_Size
            cfg = [];
            cfg.latency = [-5,-3-1/dataBipolar.fsample];
            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        end
    elseif Task_Period == 2
        %% Extract 2 seconds of maintenance
        nSet_Size = size(dataBipolar_SS,2);
        dataBipolar_Ret_SS = {};
        for iSS = 1:nSet_Size
            cfg = [];
            cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        end

    else
        %% Extract 2 seconds of retrieval latency
        nSet_Size = size(dataBipolar_SS,2);
        dataBipolar_Ret_SS = {};
        for iSS = 1:nSet_Size
            cfg = [];
            cfg.latency = [0,2-1/dataBipolar.fsample];
            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        end

    end
%% Channel pairs to run
%To be configured for every combination of electrode pairs you have

clear nChannelPairs;

if Hippocampus_Coupling
    nChannelPairs = [81 83;81 84; 82 83; 82 84];
end

if Hippocampus_Grid_Coupling
    nGridChans = 64;
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
%     nChannelPairs = [];
    nChannelPairs = [[ones(nGridChans,1),(1:nGridChans)'+AHL_chans];[nGridChans+AHL_chans+ones(nGridChans,1),(1:nGridChans)'+AHL_chans]];
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChans,1),(1:nGridChans)'+ AHL_chans]];
    end
    
end

if Scalp_Hippocampus_coupling
    nGridChans = 64;
    nScalpChans = size(dataBipolarSclp.label,1);
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
    nChannelPairs = [];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+nGridChans+AHL_chans+PHL_chans+nScalpChans+zeros(nScalpChans,1),(1:nScalpChans)']];
    end
end

if Scalp_ECoG_coupling
    nGridChans = 64;
    nScalpChans = size(dataBipolarSclp.label,1);
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
    nChannelPairs = [];
    Pz_loc =find(contains(dataBipolarSclp.label, '_Pz')==1);
    Fz_loc =find(contains(dataBipolarSclp.label, '_Fz')==1);
    P3_loc =find(contains(dataBipolarSclp.label, '_P3')==1);
    P4_loc =find(contains(dataBipolarSclp.label, '_P4')==1);

%     scalp_chans_grid_pair = [Fz_loc Pz_loc];
    scalp_chans_grid_pair = [P3_loc P4_loc Pz_loc];

    
    for i=1:length(scalp_chans_grid_pair)
        nChannelPairs = [nChannelPairs;[scalp_chans_grid_pair(i)+zeros(nGridChans,1),nScalpChans+AHL_chans+(1:nGridChans)']];
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Number of trials for each set size
nNumberOfTrialsSetSizes = zeros(1,size(dataBipolar_SS,2));
for iSS = 1:size(dataBipolar_SS,2)
    nNumberOfTrialsSetSizes(iSS) = length(dataBipolar_SS{iSS}.trial);
end


%% PLV significance for a specific task period 

if PLV_significance
    %% Random Permutations for the trials on every set size
    nNumberOfPermutations = 200;
    nRandomTrialNumberList_Pair_SS = cell(1,size(dataBipolar_SS,2));
    for iSS = 1:size(dataBipolar_SS,2)
        nTrialNumbersAll = 1:nNumberOfTrialsSetSizes(iSS);
        for nPair = 1:size(nChannelPairs,1)
            nRandomTrialNumberList = zeros(nNumberOfPermutations,nNumberOfTrialsSetSizes(iSS));
            for nRand = 1:nNumberOfPermutations
                cond = 1;
                while(cond)
                    indRand = randperm(nNumberOfTrialsSetSizes(iSS));
                    cond = ~isempty(find((indRand-(1:length(indRand)))==0,1));
                end
                nRandomTrialNumberList(nRand,:) = nTrialNumbersAll(indRand);
            end
            nRandomTrialNumberList = [nTrialNumbersAll;nRandomTrialNumberList];
            nRandomTrialNumberList_Pair_SS{iSS}{nPair} = nRandomTrialNumberList;
        end
    end
    
    %% Save variables
    strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\'
    mkdir(strPaths.Results)
    strVariableFolder = [strPaths.Results,'Statistics\Rand_Trials_PLV\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms.mat'],'nRandomTrialNumberList_Pair_SS','nNumberOfPermutations','-v7.3')
    
    %% Select pairs for analysis
    clear PLV_Pairs;
    load(Path_to_load);
    thrPLV = 0.2;  %???
    iSS = 4;
    flagRunAnalysis_Pair = zeros(size(nChannelPairs,1),1);
    for iPair = 1:size(nChannelPairs,1)
        flagRunAnalysis_Pair(iPair) = ~isempty(find(abs(PLV_Pairs{iPair,iSS})>thrPLV,1));
    end
    nPairs_ToRun = find(flagRunAnalysis_Pair);
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_'))); %% to be checked
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    
    
    %%%%
    %%
    iSS_ToCompute = 4;
    tStartRand = cputime;
    clear PLV_Rand_SS;
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,size(dataBipolar_SS,2));
        for iSS = iSS_ToCompute
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
            % For plotting purposes, otherwise there would be an error in line:522
            if nChannel_1 > nChannel_2
                tmp1 = data_SinglePair_SS{iSS}.label{1};
                data_SinglePair_SS{iSS}.label{1} = data_SinglePair_SS{iSS}.label{2};
                data_SinglePair_SS{iSS}.label{2} = tmp1;
                for tt = 1:size(data_SinglePair_SS{iSS}.trial,2)
                    tmp2 = data_SinglePair_SS{1,iSS}.trial{tt}(1,:);
                    data_SinglePair_SS{iSS}.trial{tt}(1,:) = data_SinglePair_SS{iSS}.trial{tt}(2,:);
                    data_SinglePair_SS{iSS}.trial{tt}(2,:) = tmp2;
                    
                end
            end
            
        end
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = cell(1,size(dataBipolar_SS,2));
        nChannel_ToChange = 2;
        for iSS = iSS_ToCompute
            for nRand = 1:nNumberOfPermutations+1
                data_SinglePair_Rand_SS{iSS}{nRand} = data_SinglePair_SS{iSS};
                for nTrial = 1:length(data_SinglePair_SS{iSS}.trial)
                    nTrialInRand = nRandomTrialNumberList_Pair_SS{iSS}{nPair}(nRand,nTrial);
                    data_SinglePair_Rand_SS{iSS}{nRand}.trial{nTrial}(nChannel_ToChange,:) = ...
                        data_SinglePair_SS{iSS}.trial{nTrialInRand}(nChannel_ToChange,:);
                end
            end
        end
        
        %% Phase coherence for randomized trials
        clear Freq PLV_SS
        for nRand = 1:nNumberOfPermutations+1
            for iSS = iSS_ToCompute
                %% Frequency analysis
                cfg             = [];
                cfg.method      = 'mtmfft';
                cfg.taper       = 'dpss';
                cfg.output      = 'fourier';
                cfg.foi         = 0.5:0.5:30;
                cfg.tapsmofrq   = 2;
                cfg.channel     = 1:2;
                %             cfg.pad         = 'nextpow2';
                Freq    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS}{nRand});
                
                cfgFreq = cfg;
                
                %% PLV analysis
                cfg             = [];
                cfg.method      = 'plv';
                cfg.complex     = 'complex';
                PLV_Rand_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq);
                
            end
        end
        clear Freq
        fprintf('\nChannel pair %d/%d complete\n',iPair,length(nPairs_ToRun))
        
        %% Store results for the channel pair
        iSS = iSS_ToCompute;
        strPairLabels = PLV_Rand_SS{iSS}{1}.label';
        freqAxis = PLV_Rand_SS{iSS}{1}.freq;
        for nRand = 1:length(PLV_Rand_SS{iSS})
            PLV_Rand_SS{iSS}{nRand} = squeeze(PLV_Rand_SS{iSS}{nRand}.plvspctrm(1,2,:))';
        end
        TimeInterval = [-2,-1/dataBipolar_Ret_SS{iSS}.fsample];
        
        %% Convert cell to double
        for iSS = 1:size(dataBipolar_Ret_SS,2)
            try
                PLV_Rand_SS{iSS} = cell2mat(PLV_Rand_SS{iSS}');
            catch
                PLV_Rand_SS{iSS} = [];
            end
        end
        
        %% Save values
        iSS = iSS_ToCompute;
        strPaths.var = [strPaths.Results,'Statistics\Rand_PLV_pairs\']
        mkdir(strPaths.var)
        save([strPaths.var ,'Patient_',num2str(42,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_ss',num2str(iSS),'.mat'],...
            'PLV_Rand_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
        
        %% Real values
        nRand = 1;
        Real_PLV = PLV_Rand_SS{iSS}(nRand,:);
        
        %% Null distribution
        PLV_Rand_Dist_SingleSS = PLV_Rand_SS{iSS}(2:end,:); %??
        
        %% Percentiles
        prcList = [1,5,50,95,99];
        RandPrc = {};
        RandPrc_Im = {};
        for iPrc = 1:length(prcList)
            val  = prctile(abs(PLV_Rand_Dist_SingleSS),prcList(iPrc));
            RandPrc{prcList(iPrc)} = val;
            val  = prctile(imag(PLV_Rand_Dist_SingleSS),prcList(iPrc));
            RandPrc_Im{prcList(iPrc)} = val;
        end
        
        %% Min max
        RandStats.Min = min(abs(PLV_Rand_Dist_SingleSS));
        RandStats.Max = max(abs(PLV_Rand_Dist_SingleSS));
        RandStats.Mean = mean(abs(PLV_Rand_Dist_SingleSS));
        RandStats.Std = std(abs(PLV_Rand_Dist_SingleSS));
        RandStats.NoPermutations = size(PLV_Rand_Dist_SingleSS,1);
        
        RandStats_Im.Min = min(imag(PLV_Rand_Dist_SingleSS));
        RandStats_Im.Max = max(imag(PLV_Rand_Dist_SingleSS));
        RandStats_Im.Mean = mean(imag(PLV_Rand_Dist_SingleSS));
        RandStats_Im.Std = std(imag(PLV_Rand_Dist_SingleSS));
        RandStats_Im.NoPermutations = size(PLV_Rand_Dist_SingleSS,1);
        
        %% Save results
        strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\']
        mkdir(strPaths.stats)
        save([strPaths.stats,'Rand_Stats_Patient_',num2str(42,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'],...
            'Real_PLV','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
        clear Real_PLV RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis cfgFreq TimeInterval PLV_Rand_SS
        
    end
    tStopRand = cputime-tStartRand;
    
    
    %% Load results / PLV for all pairs / Randomized trials / 2 seconds of retention
    % strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Stats 180405\'];
    PLV_Rand_Ret_Variables_Pair_SS = cell(1,size(dataBipolar_Ret_SS,2));
    for iSS = 1:size(dataBipolar_Ret_SS,2)
        PLV_Rand_Ret_Variables_Pair_SS{iSS} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            strVariablePath = [strPaths.stats,'Rand_stats_Patient_',num2str(42,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'];
            try
                Vars = load(strVariablePath);
                PLV_Rand_Ret_Variables_Pair_SS{iSS}{nPair} = Vars;
                clear Vars
            catch
            end
        end
    end
    
    %% Plot results for each depth electrode for randomized values
    iSS = 4;
    % strImageFolder = [strPaths.Sternberg,'Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Scalp Layout Plots 180320\'];
    figLayout = Get_Figure_Scalp_Layout_Information();
    % Initialize plot
    Max_PLV_Pairs_Table = zeros(size(nChannelPairs,1),6);
    Max_PLV_Pairs_Table = array2table(Max_PLV_Pairs_Table);
    Max_PLV_Pairs_Table.Properties.VariableNames = {'Scalp','ECoG','Max_PLV','freq_Max_PLV','AreaUnder_PLV','freqBand_AreaUnder_PLV'};
    Max_PLV_Pairs_Table.Scalp = cell(size(nChannelPairs,1),1);
    Max_PLV_Pairs_Table.Depth = cell(size(nChannelPairs,1),1);
    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV = cell(size(nChannelPairs,1),1);
    %% Find the significance comparing to the precalculated PLV_Pairs
    %Check the script Phase_Coherence_Patient_42_DS
    if Hippocampus_Grid_Coupling
        PLV_Pairs = PLV_Pairs; % Change this to include the Pairs for a specific task period or for a specific coupling type
    end
    if Scalp_Hippocampus_coupling
        PLV_Pairs = PLV_Pairs_Scalp_maint_HC;
    end
    nChannelList_Grid = unique([nChannelPairs(nPairs_ToRun,2)]); %%tbc;
    clear SubAxes;
    SubAxes = cell(1,6);
    
    for iChannel_grid = 1:length(nChannelList_Grid) % 13
        % fig = figure;
        % Make_Plot_Secondary_Display_Fullscreen_Gca;
        % Grid channel
        nChannel_grid = nChannelList_Grid(iChannel_grid);
        strGridChannelName = strChannelNameList{nChannel_grid};
        % Scalp channels
        indGridChannelPairs = find(nChannelPairs(:,2)==nChannel_grid);
        nChannelList_Scalp = nChannelPairs(indGridChannelPairs,1);
        strScalpChannelNameList = strChannelNameList(nChannelList_Scalp);
        %%  Plot results on layout
        for iScalp = 1:length(strScalpChannelNameList)
            if iChannel_grid == 1
                %             Make_Plot_Secondary_Display_Fullscreen_Gca;
                if Hippocampus_Grid_Coupling
                    figure
                    ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
                    SubAxes{iScalp} = ha;
                end
                if Scalp_Hippocampus_coupling
                    %                 clear nSubplot_psd;
                    figure;
                    set(gcf,'color','white')
                    ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
                    SubAxes{iScalp} = ha;
                end
            end
            % Scalp channel number and label
            nChannel_Scalp = nChannelList_Scalp(iScalp);
            strScalpChannelName = strChannelNameList{nChannel_Scalp};
            %% Pair
            nPair = find(nChannelPairs(:,1)== nChannel_Scalp&nChannelPairs(:,2)==nChannel_grid);
            %% Subplot in layout
            %         indSubplot = find(strcmpi(figLayout(:,1),strrep(sprintf('%s',strScalpChannelName),'_','')));
            %         nSubplot = figLayout{indSubplot,2};
            %         subplot(5,5,nSubplot)
            Vars = PLV_Rand_Ret_Variables_Pair_SS{iSS}{nPair};
            fig = figure(iScalp);
            
            if(~isempty(Vars))
                % Frequency axis
                freqAxis = Vars.freqAxis;
                % Real PLV
                %             Real_PLV = Vars.Real_PLV;
                Real_PLV = PLV_Pairs{nPair,iSS};
                
                
                %% Percentiles
                prcList = [1,5,50,95,99];
                for iPrc = 1:length(prcList)
                    Vars.RandPrc{prcList(iPrc)};
                    Vars.RandPrc{prcList(iPrc)};
                end
                %
                %% Plot
                if Hippocampus_Grid_Coupling
                    axes(SubAxes{iScalp}(iChannel_grid));
                end
                if Scalp_Hippocampus_coupling
                    Scalp_channel = dataBipolar.label{iChannel_grid};
                    Scalp_channel = strrep(sprintf('%s',Scalp_channel),'_',' ');
                    indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
                    nSubplot_psd(iChannel_grid,iScalp) = figLayout{indSubplot,2};
                    axes(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)));
                end
                plot(freqAxis,abs(Real_PLV),strPlotColors{iSS}); ylim([0 1])
                if Hippocampus_Grid_Coupling
                    set(SubAxes{iScalp}(iChannel_grid),'XTickLabel',''); set(SubAxes{iScalp}(iChannel_grid),'YTickLabel','')
                end
                if Scalp_Hippocampus_coupling
                    set(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)),'XTickLabel',''); set(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)),'YTickLabel','')
                end
                
                if mod(iChannel_grid,64) == 0
                    spt = [strrep(sprintf('%s',strScalpChannelName),'_',' '),' - Grid Coupling'];
                    suptitle(spt);
                end
                ylabel(sprintf('%s',strrep(sprintf('%s',strGridChannelName),'_',' ')));
                
                hold on
                %             shadedErrorBar(freqAxis,Vars.RandPrc{50},[Vars.RandPrc{95}-Vars.RandPrc{50};Vars.RandPrc{50}-Vars.RandPrc{5}],{'Color',0.5*[1,1,1]});
                %             plot(freqAxis,Vars.RandStats.Mean)
                %             plot(freqAxis,Vars.RandStats.Mean+Vars.RandStats.Mean,'m')
                %             plot(freqAxis,Vars.RandStats.Mean-Vars.RandStats.Mean,'m')
                %             plot(freqAxis,Vars.RandPrc{99},'g')
                %             plot(freqAxis,abs(Real_PLV),strPlotColors{iSS})
                
                if(~(strcmpi(Vars.strPairLabels{1},strScalpChannelName)&&strcmpi(Vars.strPairLabels{2},strGridChannelName)))
                    error('Channel labels do not match')
                end
                
                %% Axis limits
                %             ylim([0,1])
                %             xlim([1,30])
                %             set(gca,'XTick',[1,5:5:30])
                %             title(sprintf('%s - %s',strrep(sprintf('%s',strScalpChannelName),'_',' '),strrep(sprintf('%s',strGridChannelName),'_',' ')))
                
                %% Channel names
                Max_PLV_Pairs_Table.Scalp{nPair} = strScalpChannelName;
                Max_PLV_Pairs_Table.Depth{nPair} = strGridChannelName;
                
                %% Max PLV
                [valMax,indMax] = max(abs(Real_PLV));
                Max_PLV_Pairs_Table.Max_PLV(nPair) = valMax;
                Max_PLV_Pairs_Table.freq_Max_PLV(nPair) = Vars.freqAxis(indMax);
                
                %% Area under curve
                thrPrc = 95;
                thrNeighborPts = 3;
                thrFreq = 0;
                Real_PLV = abs(Real_PLV);
                Thr_PLV = Vars.RandPrc{thrPrc};
                indAboveThr = Real_PLV>=Thr_PLV;
                cc = bwconncomp(indAboveThr);
                if(cc.NumObjects~=0)
                    % Threshold number of points
                    numPixels = cellfun(@numel,cc.PixelIdxList);
                    cc.PixelIdxList = cc.PixelIdxList(numPixels>=thrNeighborPts);
                    % Threshold lower bound for frequency
                    cc.PixelIdxList = cc.PixelIdxList(freqAxis(cellfun(@min,cc.PixelIdxList))>=thrFreq);
                    % Largest cluster
                    %                 numPixels = cellfun(@numel,cc.PixelIdxList);
                    %                 [~,indMaxCluster] = max(numPixels);
                    % Cluster with the highest area under
                    plv_all = cellfun(@(x) Real_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    thr_plv_all = cellfun(@(x) Thr_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    [~,indMaxCluster] = max(cellfun(@sum,cellfun(@(x,y) x-y,plv_all,thr_plv_all,'UniformOutput',0)));
                    if(~isempty(indMaxCluster))
                        indMaxBand = cc.PixelIdxList{indMaxCluster};
                        freqBand = freqAxis(indMaxBand);
                        AreaUnder_PLV_Prc = sum(Real_PLV(indMaxBand)-Thr_PLV(indMaxBand));
                        % stem(freqAxis(indMaxBand),Real_PLV(indMaxBand),'Color','magenta')
                        % Plot the bar for PLV significance
                        %                   axes(SubAxes{iScalp}(iChannel_grid));
                        
                        
                        if Hippocampus_Grid_Coupling
                            axes(SubAxes{iScalp}(iChannel_grid));
                        end
                        if Scalp_Hippocampus_coupling
                            Scalp_channel = dataBipolar.label{iChannel_grid};
                            Scalp_channel = strrep(sprintf('%s',Scalp_channel),'_',' ');
                            indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
                            nSubplot_psd(iChannel_grid,iScalp) = figLayout{indSubplot,2};
                            axes(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)));
                        end
                        Plv_sgnf_struct.plv_sign_band = freqAxis(indMaxBand);
                        Plv_sgnf_struct.label = strScalpChannelName;
                        PLV_significance{iScalp,iChannel_grid} = Plv_sgnf_struct;
                        
                        plot(freqAxis(indMaxBand),0.05*ones(length(freqAxis(indMaxBand)),1)','Color','magenta','LineWidth',2)
                        
                        % If a frequency band is given
                        FreqBand = [6.5,18.5];
                        indFreqBand = freqAxis>=FreqBand(1)&freqAxis<=FreqBand(2);
                        AreaUnder_PLV_Prc_SelectedBand = sum(Real_PLV(indFreqBand)-Thr_PLV(indFreqBand));
                    else
                        freqBand = [];
                        AreaUnder_PLV_Prc = 0;
                        AreaUnder_PLV_Prc_SelectedBand = 0;
                    end
                    Max_PLV_Pairs_Table.AreaUnder_PLV(nPair) = AreaUnder_PLV_Prc;
                    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV{nPair} = freqBand;
                    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV_Min{nPair} = min(freqBand);
                    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV_Max{nPair} = max(freqBand);
                    Max_PLV_Pairs_Table.AreaUnder_PLV_Prc_SelectedBand(nPair) = AreaUnder_PLV_Prc_SelectedBand;
                    latency = [dataBipolar_Ret_SS{iSS}.time{1}(1) dataBipolar_Ret_SS{iSS}.time{1}(end) ];
                    %                 pl = area(freqAxis(indMaxBand),x(indMaxBand));
                    %                 title(sprintf('%s - %s - %.2f %s %s %s %s %s %s',strrep(sprintf('%s',strScalpChannelName),'_',' '),strrep(sprintf('%s',strGridChannelName),'_',' '),AreaUnder_PLV_Prc ....
                    %                     , ' set size [', num2str(Set_Sizes{iSS}),'] Time [',num2str(latency(1)),num2str(latency(2)),']'));
                    ylabel(sprintf('%s',strrep(sprintf('%s',strGridChannelName),'_',' ')));
                    flagSig_SS6(nPair) = length(cc.PixelIdxList);
                    flagSig_SS6(nPair) = length(cc.PixelIdxList);
                    if Hippocampus_Grid_Coupling
                        set(SubAxes{iScalp}(iChannel_grid),'XTickLabel',''); set(SubAxes{iScalp}(iChannel_grid),'YTickLabel','')
                    end
                    if Scalp_Hippocampus_coupling
                        set(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)),'XTickLabel',''); set(SubAxes{iScalp}(nSubplot_psd(iChannel_grid,iScalp)),'YTickLabel','')
                    end
                    str1 = sprintf('%s - %s',strrep(sprintf('%s',strScalpChannelName),'_',' '),strrep(sprintf('%s',strGridChannelName),'_',' '));
                    str_legend1 = ['PLV spectrum of the pair'];
                    str_legend2 = ['PLV significance ',num2str(thrPrc),'-th percentile'];
                    %                 legend(str_legend1,str_legend2)
                    ylim([0 1])
                    
                    Figure_name = [Path_for_PLV_sgnf,'PLV_significance_',num2str(thrPrc),'-th percentile',strScalpChannelName,'-',strGridChannelName];
                    %                 saveas(fig,Figure_name,'png')
                    
                end
                
                
            end
            if iChannel_grid == 1 && Scalp_Hippocampus_coupling
                strScalpChannelName = strrep(sprintf('%s',strScalpChannelName),'_',' ');
                suptitle(['PLV significance ',strScalpChannelName,'- Scalp',' Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
                for i = 1:size(ha,1)
                    
                    if ~ismember(i,nSubplot_psd)
                        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
                    end
                    
                end
            end
        end
        
    end
end

%% Maintenance vs Encoding and Maint vs Fix
if PLV_significance_maint_vs_encod
   %% Dont forget to adapt the names in case the cfg.foi = 0.5:10:120
    nSS_Ran = [4,4];
  
    nNumberOfPermutations = 200;
    nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
    nTrialNumbersAll_AcrossSS = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
    nRandomTrialNumberList_AcrossTime_Pairs = cell(1,size(nChannelPairs,1));
    for iPair = 1:size(nChannelPairs,1)
        nPair = iPair; % nPairs_ToRun(iPair);
        nRandomTrialNumberList = zeros(nNumberOfPermutations,sum(nNumberOfSetSizes_ToCompare));
        for nRand = 1:nNumberOfPermutations
            cond = 1;
            while(cond)
                indRand = randperm(sum(nNumberOfSetSizes_ToCompare));
                cond = ~isempty(find((indRand-(1:length(indRand)))==0,1));
            end
            nRandomTrialNumberList(nRand,:) = nTrialNumbersAll_AcrossSS(indRand);
        end
        nRandomTrialNumberList = [nTrialNumbersAll_AcrossSS;nRandomTrialNumberList];
        nRandomTrialNumberList_AcrossTime_Pairs{nPair} = nRandomTrialNumberList;
    end
    
    %% Save Permutation numbers
    strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\'
    mkdir(strPaths.Results)
    strVariableFolder = [strPaths.Results,'Statistics\Rand_Trials_PLV\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms_maint_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')
    
    %% Extract 2 seconds of encoding or 1s of fixation
        dataBipolar_Enc_SS = cell(1,size(Set_Sizes,1));
        for iSS = 1:size(Set_Sizes,1)
            cfg = [];
%             cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
            cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample]; % fixation vs maintenance
            dataBipolar_Enc_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
            if PLV_significance_maint_vs_fix
                cfg =[];
                cfg.resamplefs = 1000;
                dataBipolar_Enc_SS{iSS} = ft_resampledata(cfg,dataBipolar_Enc_SS{iSS});
            end
        end
        
        
    %% Pairs to run analysis for
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','_')));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_GL')));
    
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    strVariableFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_vs_Encoding\';
    mkdir(strVariableFolder)
    %%
    tStartRand = cputime;
    for iPair = 1:length(nPairs_ToRun)
        disp(iPair)

        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        

        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,2);
        Period = 1;
        iSS = nSS_Ran(Period);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{Period} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
        
        Period = 2;
        iSS = nSS_Ran(Period);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{Period} = ft_preprocessing(cfg,dataBipolar_Enc_SS{iSS});
    
        
        %% Randomly mix trials between encoding & maintenance

        for nRand = 1:nNumberOfPermutations+1
            
            %% Create datasets with randomized trials
            data_SinglePair_Rand_SS = data_SinglePair_SS;
            iTrialList = nRandomTrialNumberList_AcrossTime_Pairs{nPair}(nRand,:);
            for iTrial = 1:length(iTrialList)
                nTrial = iTrialList(iTrial);
                if(iTrial<=nNumberOfSetSizes_ToCompare(1)) % assign values for set size 1
                    iSS_Rand = 1;
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                else % comes from set size 2
                    iSS_Rand = 2;
                    iTrial = iTrial-nNumberOfSetSizes_ToCompare(1);
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                end
                data_SinglePair_Rand_SS{iSS_Rand}.trial{iTrial} = ...
                    data_SinglePair_SS{iSS_Real}.trial{abs(nTrial)};
            end
            
            %% Phase coherence for randomized trials
            clear Freq_SS
            fprintf('Calculating %s with  %s right now. Elapsed pairs are %d..... \n',strChannelNameList{nChannel_1},strChannelNameList{nChannel_2},length(nPairs_ToRun)-iPair);

            for iSS = 1:2
                %% Frequency analysis
                cfg             = [];
                cfg.method      = 'mtmfft';
                cfg.taper       = 'dpss';
                cfg.output      = 'fourier';
                cfg.tapsmofrq   = 2;
                cfg.pad         = 2;
                if PLV_significance_maint_vs_fix
%                     cfg.foi         = 0.5:0.5:120;
                cfg.foi         = 1:1:30;
                else
                    cfg.foi         = 0.5:0.5:30;
                end
                cfg.channel     = 1:2;
                Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
                
                %% PLV analysis
                cfg             = [];
                cfg.method      = 'plv';
                cfg.complex     = 'complex';
                PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            end
        end
          clear Freq_SS
          
          %% Store results for the channel pair
          strPairLabels = PLV_SS{iSS}{1}.label';
          freqAxis = PLV_SS{iSS}{1}.freq;
          for iSS = 1:2
              for nRand = 1:nNumberOfPermutations+1
                  PLV_AcrossPeriods{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
              end
          end
          clear PLV_SS
          
          %% Save results
          strPLVResultFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_vs_Encoding\';
        
          if Scalp_Hippocampus_coupling
             strVariableFolder = [strPLVResultFolder,'Scalp_Hippocampus\'];
             mkdir(strVariableFolder);
             save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_freq_1_1_30_Maint_vs_Fix_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
              'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
          elseif Hippocampus_Grid_Coupling
              if PLV_significance_maint_vs_fix
                  strVariableFolder = [strPLVResultFolder,'Grid_Hippocampus\Fix_vs_maint\'];
                  mkdir(strVariableFolder)
                  save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Fix_vs_maint_mt_tsf_2_Set_Size_freq_1_30_',num2str(2*iSS+2,'%.1d'),'.mat'],...
                      'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
              else
                  strVariableFolder = [strPLVResultFolder,'Grid_Hippocampus\'];
                  mkdir(strVariableFolder)
                  save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
                      'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
              end
         elseif Scalp_ECoG_coupling  
               if PLV_significance_maint_vs_fix
                  strVariableFolder = [strPLVResultFolder,'Scalp_Grid\Fix_vs_maint\'];
                  mkdir(strVariableFolder)
                  save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_Fix_vs_maint_mt_tsf_2_Set_Size_freq_1_30_',num2str(2*iSS+2,'%.1d'),'.mat'],...
                      'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
               end
          end
                  
          
          %% Real values
          nRand = 1;
          for iSS = 1:2
              Real_PLV_AcrossPeriods{iSS} = PLV_AcrossPeriods{iSS}{nRand};
          end
          %% Null distribution
          for iSS = 1:2
              RandDist_PLV_All_AcrossPeriods{iSS} = [];
              for nRand = 1:length(PLV_AcrossPeriods{iSS})-1
                  RandDist_PLV_All_AcrossPeriods{iSS}(nRand,:) = PLV_AcrossPeriods{iSS}{nRand+1};
              end
          end
          clear PLV_AcrossPeriods
          
          PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossPeriods{1})-abs(RandDist_PLV_All_AcrossPeriods{2});
          RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossPeriods{1})-imag(RandDist_PLV_All_AcrossPeriods{2});
         
          %% Percentiles
          prcList = [1,5,50,95,99];
          RandPrc = {};
          RandPrc_Im = {};
          for iPrc = 1:length(prcList)
              val  = prctile(PLV_Rand_Dist_SingleSS,prcList(iPrc));
              RandPrc{prcList(iPrc)} = val;
              val  = prctile(RandDist_PLV_All_Im,prcList(iPrc));
              RandPrc_Im{prcList(iPrc)} = val;
          end
          
          %% Min max
          RandStats.Min = min(PLV_Rand_Dist_SingleSS);
          RandStats.Max = max(PLV_Rand_Dist_SingleSS);
          RandStats.Mean = mean(PLV_Rand_Dist_SingleSS);
          RandStats.Std = std(PLV_Rand_Dist_SingleSS);
          RandStats.NoPermutations = size(PLV_Rand_Dist_SingleSS,1);
          
          RandStats_Im.Min = min(RandDist_PLV_All_Im);
          RandStats_Im.Max = max(RandDist_PLV_All_Im);
          RandStats_Im.Mean = mean(RandDist_PLV_All_Im);
          RandStats_Im.Std = std(RandDist_PLV_All_Im);
          RandStats_Im.NoPermutations = size(RandDist_PLV_All_Im,1);
          
          %% Save results
          if Scalp_Hippocampus_coupling
              if PLV_significance_maint_vs_fix
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Hippocampus\Maint_vs_Fix\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
              else
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Hippocampus\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
              end
          elseif Hippocampus_Grid_Coupling
              if PLV_significance_maint_vs_fix
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
                  
              else
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
              end
          elseif Scalp_ECoG_coupling
              if PLV_significance_maint_vs_fix
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Grid\Fix_vs_Maint\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Scalp',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
                  
              end
          end
          
    end
    clear PLV_SS;
    tStopRand = cputime-tStartRand;
   %% Load the results 
   if Scalp_Hippocampus_coupling
       if PLV_significance_maint_vs_fix
           strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Hippocampus\Maint_vs_Fix\']
       else
           strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\']
       end
%        strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\']
   elseif Hippocampus_Grid_Coupling
       if PLV_significance_maint_vs_fix
           strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
           
       else
           strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\']
       end
       elseif Scalp_ECoG_coupling
           if PLV_significance_maint_vs_fix
                 strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Grid\Fix_vs_Maint\'];
           end
           
   end
    PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,size(Set_Sizes,1));
    for iSS = 1:2%size(Set_Sizes,1)
        PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            if Scalp_Hippocampus_coupling
                if PLV_significance_maint_vs_fix
                    strVariablePath =  [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation','.mat']
                else
                    strVariablePath =  [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'];

                end
            elseif Hippocampus_Grid_Coupling
                if PLV_significance_maint_vs_fix
                    strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat'];
                else
                    strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'];

                end
%                 strVariablePath =  [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'];
            
            elseif Scalp_ECoG_coupling
                if PLV_significance_maint_vs_fix
                    strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Scalp',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat'];
                end
            
            end

%             [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'];
            try
                Vars = load(strVariablePath);
                PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPair} = Vars;
                clear Vars
            catch
            end
        end
    end
    %% Load PLV_Pairs
   if Scalp_Hippocampus_coupling
       
       if PLV_significance_maint_vs_fix
           PLV_Path_maint_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_maintenance_cleanEEG.mat';
           PLV_Path_fix_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_fixation_cleanEEG.mat';
            load(PLV_Path_maint_sclp);
            load(PLV_Path_fix_sclp);
       elseif PLV_significance_maint_vs_encod
            load(PLV_Path_maint);
            load(PLV_Path_enc);
       end
   elseif Scalp_ECoG_coupling
       if PLV_significance_maint_vs_fix
           PLV_Path_maint_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_maintenance_clean_EEG_f_1_30_all_scalp_chans.mat';
           PLV_Path_fix_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_fixation_clean_EEG_f_1_30_all_scalp_chans.mat';
           load(PLV_Path_maint_sclp);
           load(PLV_Path_fix_sclp);
       end
   end
    
    %% Display results for every pair Maint vs Enc Scalp Hippocampus
%     nPair = 58%nPairs_ToRun;
    saveFlag = 0;
    figLayout = Get_Figure_Scalp_Layout_Information();
    freqBand_topo = [8 12];
    [~,indFreq1] = min(abs(freqAxis_sclp-freqBand_topo(1)));
    [~,indFreq2] = min(abs(freqAxis_sclp-freqBand_topo(2)));

    for nPair = 1:size(nPairs_ToRun,2)
        iSS = 4;
        if Scalp_Hippocampus_coupling
            Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair};
            Maint_Encod_Diff = abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS})-abs(PLV_Pairs_Scalp_enc_HC{nPair,iSS});
            datavector(nPair) = median(Maint_Encod_Diff(indFreq1:indFreq2))
        end
        if  nChannelPairs(nPair,2) == 1 %For every different hippocampal channel paired with Scalp channels
            fig = figure
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            set(gcf,'color','white')
            ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
            
        end
        Diff_significance = find(Maint_Encod_Diff>=Vars.RandPrc{95}); %Comparing with 95th percentile
        % Area under 95th percentile
        thrPrc = 95;
        thrNeighborPts = 1;
        thrFreq = 0;
        Real_PLV = Maint_Encod_Diff;
        Thr_PLV = Vars.RandPrc{thrPrc};
        indAboveThr = Real_PLV>=Vars.RandPrc{95};
        %     fig = figure;
        Scalp_channel = strrep(dataBipolarSclp.label{nChannelPairs(nPair,2)},'_',' ');
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_maint(nPair) = figLayout{indSubplot,2};
        axes(ha(nSubplot_maint(nPair)));
        cc = bwconncomp(indAboveThr);
        if(cc.NumObjects~=0)
            % Threshold number of points
            numPixels = cellfun(@numel,cc.PixelIdxList);
            cc.PixelIdxList = cc.PixelIdxList(numPixels>=thrNeighborPts);
            % Threshold lower bound for frequency
            cc.PixelIdxList = cc.PixelIdxList(freqAxis_sclp(cellfun(@min,cc.PixelIdxList))>=thrFreq);
            % Largest cluster
            % numPixels = cellfun(@numel,cc.PixelIdxList);
            % [~,indMaxCluster] = max(numPixels);
            % Cluster with the highest area under
            plv_all = cellfun(@(x) Real_PLV(x),cc.PixelIdxList,'UniformOutput',0);
            thr_plv_all = cellfun(@(x) Thr_PLV(x),cc.PixelIdxList,'UniformOutput',0);
            [~,indMaxCluster] = max(cellfun(@sum,cellfun(@(x,y) x-y,plv_all,thr_plv_all,'UniformOutput',0)));
            if(~isempty(indMaxCluster))
                indMaxBand = cc.PixelIdxList{indMaxCluster};
                freqBand = freqAxis_sclp(indMaxBand);
                AreaUnder_PLV_Prc = sum(Real_PLV(indMaxBand)-Thr_PLV(indMaxBand));
                plot(freqAxis_sclp(indMaxBand),-0.05*ones(length(freqAxis_sclp(indMaxBand)),1)','Color','k','LineWidth',2)
                ylim([-0.1,1]);
                hold on;
                
            else
                freqBand = [];
                AreaUnder_PLV_Prc = 0;
                AreaUnder_PLV_Prc_SelectedBand = 0;
            end
            
        end
        
        nChannel1 = strrep( strChannelNameList(nChannelPairs(nPair,1)),'_',' ');
        nChannel2 = strrep( strChannelNameList(nChannelPairs(nPair,2)),'_',' ');
        
        %Visualization
        %     plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossPeriods{1})-abs(Vars.Real_PLV_AcrossPeriods{2}),'LineWidth',2)
        plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS}),'r','LineWidth',2)
        ylim([-0.1,1]);

        hold on
        plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{nPair,iSS}),'g','LineWidth',2)
        %     plot(freqAxis_sclp(indMaxBand),0.05*ones(length(freqAxis_sclp(indMaxBand)),1)','Color','k','LineWidth',2)
        %     plot(freqAxis_sclp(Diff_significance),0.05*ones(length(freqAxis_sclp(Diff_significance)),1)','Color','k','LineWidth',2)
        ylim([-0.1,1])
        
        %     plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossPeriods{1}))
        %     plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossPeriods{2}))
        
        %     plot(Vars.freqAxis,Vars.RandPrc{95});
        strTitle = strcat('Maint/Encod Comparison Pair: ',nChannel1,'-', ...
            'Scalp',' Set Size [',num2str(Set_Sizes{iSS}),' ]');
        ylabel(nChannel2);
              
        %         title(strTitle);
        %         ylabel('PLV');
        %         xlabel('Frequency(Hz)');
        %         legend('Difference Significant','Maintenance','Encoding')



        %Save figure
        Fig_dir = [Path_for_PLV_sgnf,'Patient_42_Maint_Encod_diff\'];
        mkdir(Fig_dir);
        Figure_Name = strcat(Fig_dir,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),strjoin(nChannel1),'_','Scalp','_Enc_Maint');
        
        
        
        if ~mod(nPair,size(dataBipolarSclp.label,1)) %% When a hippocampal is paired with every scalp channel 
            suptitle(strTitle)
            set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
            set(ha(1:size(ha,1)),'FontSize',18);
            for i = 1:size(ha,1)
                
                if ~ismember(i,nSubplot_maint)
                    set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
                end
                
            end
            if saveFlag
                saveas(fig,Figure_Name,'png');
                saveas(fig,Figure_Name,'fig');
            end
        end
    end
   
    % Topoplot of the maint - encod in a frequency band of interest
    nScalpChans = size(dataBipolarSclp.label,1);
    nHippChans = length(Bip_chans);
    datavector = reshape(datavector',nScalpChans,nHippChans);
    dataBipolarSclp.label =  strrep(dataBipolarSclp.label,'_',' ');
%     Scalp_topoplot(strPaths.Toolboxes,freqBand_topo,freqAxis_sclp,datavector(:,3),dataBipolarSclp,4,nChannelPairs,1)

    %% Scalp Hippocampus visualize for every pair maint vs fixation
     freqBand_topo = [5 9];
    [~,indFreq1] = min(abs(freqAxis_sclp-freqBand_topo(1)));
    [~,indFreq2] = min(abs(freqAxis_sclp-freqBand_topo(2)));
    
    saveFlag = 0;
    figLayout = Get_Figure_Scalp_Layout_Information();
     for nPair = 61%1:size(nPairs_ToRun,2)%61
        iSS = 4;
        if Scalp_Hippocampus_coupling
            Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair};
            Maint_Fix_Diff = abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS})-abs(PLV_Pairs_Scalp_fix_HC{nPair,iSS});
            datavector(nPair) = median(abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS}(indFreq1:indFreq2)))
        end
        if  nChannelPairs(nPair,2) == 1 %For every different hippocampal channel paired with Scalp channels
            fig = figure
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            set(gcf,'color','white')
%             ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
            
        end
        Scalp_channel = strrep(dataBipolarSclp.label{nChannelPairs(nPair,2)},'_',' ');
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_maint(nPair) = figLayout{indSubplot,2};
%         axes(ha(nSubplot_maint(nPair)));
        nChannel1 = strrep( strChannelNameList(nChannelPairs(nPair,1)),'_',' ');
        nChannel2 = strrep( strChannelNameList(nChannelPairs(nPair,2)),'_',' ');
        
        %Visualization
        %     plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossPeriods{1})-abs(Vars.Real_PLV_AcrossPeriods{2}),'LineWidth',2)
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS}),'k','LineWidth',3)
        ylim([0,1]);

        hold on
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_HC{nPair,iSS}),'r','LineWidth',3)
        hold on;
        indMaxBandMaintFix = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis_sclp,[0 1], 0.05,'k');
        ylim([0,1])
        strTitle = strcat('Maint/Fix Significance Pair: ',nChannel1,'-', ...
            'Scalp',' Set Size [',num2str(Set_Sizes{iSS}),' ]');
        ylabel(nChannel2);
        xlim([3 100])
       Fig_dir = [Path_for_PLV_sgnf,'Maint vs Fix\','Patient_42_Maint_Encod_diff\'];
        mkdir(Fig_dir);
        Figure_Name = strcat(Fig_dir,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),strjoin(nChannel1),'_','Scalp','_Fix_Maint');
        
        
        
        if ~mod(nPair,size(dataBipolarSclp.label,1)) %% When a hippocampal is paired with every scalp channel 
%             suptitle(strTitle)
%             set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
%             set(ha(1:size(ha,1)),'FontSize',18);
            for i = 1:size(ha,1)
                
                if ~ismember(i,nSubplot_maint)
                    set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
                end
                
            end
            if saveFlag
                saveas(fig,Figure_Name,'png');
                saveas(fig,Figure_Name,'fig');
            end
        end
        
     end
     
     
datavector = [0.57 0.56  0.50 0.47 0.27 0.33 0.35 0.4 0.51 0.35 0.57 0.37 0.49 0.51 0.51 0.47 0.33 0.27 0.27 0.54 0.49 0.58 0.56]'     
     dataBipolarSclp.label = strrep(dataBipolarSclp.label,'_',' ');
%      Scalp_topoplot(strPaths.Toolboxes,freqBand_topo,freqAxis_sclp,datavector(47:69),dataBipolarSclp,4,nChannelPairs,1)
     Scalp_topoplot(strPaths.Toolboxes,freqBand_topo,freqAxis_sclp,datavector,dataBipolarSclp,4,nChannelPairs,1)

     %% Scalp Grid Visualization of every pair
     iSS_to_plot =4;
     for i = 1%1:size(nPairs_ToRun,2)/nGridChans
         
         fig = figure;
%          ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
         plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
        for j =30%1:nGridChans
            Scalp_channel = strrep(dataBipolar.label{nChannelPairs((j+(i-1)*nGridChans),1)},'_',' ');
            s = nChannelPairs_Sclp(find(nChannelPairs_Sclp== nChannelPairs((j+(i-1)*nGridChans),1)),1);
            s_loc = s(1);
            Grid_Channel = strrep(dataBipolar.label{nChannelPairs((j+(i-1)*nGridChans),2)},'_',' ');
%             axes(ha(plot_order(j)));
            semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{(j+(s_loc-1)*nGridChans),iSS_to_plot}),'k','LineWidth',3);
            hold on;
%             semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{(j+(s_loc-1)*nGridChans),iSS_to_plot}),'g','LineWidth',3);

            semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix{(j+(s_loc-1)*nGridChans),iSS_to_plot}),'r','LineWidth',3);
            ylim([0 1])
            xlim([3 100]);
            hold on;
            Maint_Fix_Diff = abs(PLV_Pairs_Scalp_maint{(j+(s_loc-1)*nGridChans),iSS_to_plot}) - abs(PLV_Pairs_Scalp_fix{(j+(s_loc-1)*nGridChans),iSS_to_plot});
            Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS_to_plot}{(j+(i-1)*nGridChans)};
            indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis_sclp,[0 1],0.05,'k');
            ylabel(Grid_Channel);
        end
        titl = ['PLV significance Maintenance - Fixation Scalp Channel',Scalp_channel,' to Grid'];
        suptitle(titl);
        strAnalysisResults = 'F:\Vasileios\Task Analysis\Analysis Results\';
        strFigFolder = [strAnalysisResults,'Maint_vs_Encod_vs_Fix\Scalp-Grid\'];
        Fig_name = strcat(strFigFolder,'PLV significance Maintenance - Fixation Scalp Channel',Scalp_channel,' to Grid');
        mkdir(strFigFolder);
        saveas(fig,Fig_name,'png')
        saveas(fig,Fig_name,'fig')
     end
     
     
     
     
     
     
     
     
    %% Grid Hippocampus PLV significance visualization for all pairs
    %% Load PLV_Pairs
    PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_gamma_smooth.mat'
    PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_gamma_smooth.mat'
    PLV_Pairs_maint_grid =  load(PLV_Path_maint);
    PLV_Pairs_maint_grid = PLV_Pairs_maint_grid.PLV_Pairs;
    PLV_Pairs_enc_grid = load(PLV_Path_enc);
    PLV_Pairs_enc_grid = PLV_Pairs_enc_grid.PLV_Pairs;
    strPlotColors = {'r','g','b','k','c'};
    for j = 4;%3:8
    fig = figure;
    freqAxis = [0.5:10:120];%[0.5:0.5:120]
    ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
    
    for i = 1:nGridChans
        axes(ha(plot_order(i)));
        iSS = 4;
        plot(freqAxis,abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS}),strPlotColors{4},'LineWidth',1);       
        hold on;
        plot(freqAxis,abs(PLV_Pairs_enc_grid{i+((j-1)*nGridChans),iSS}),strPlotColors{2},'LineWidth',1);       
        Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{i+((j-1)*nGridChans)};
        Maint_Encod_Diff = abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS})-abs(PLV_Pairs_enc_grid{i+((j-1)*nGridChans),iSS});
        indMaxBand = Difference_Bar(Maint_Encod_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'m');
        Encod_Maint_Diff = abs(PLV_Pairs_enc_grid{i+((j-1)*nGridChans),iSS})-abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS});
        color2 = [0.75, 0.75, 0];

        indMaxBand = Difference_Bar(Encod_Maint_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'c');%color2)

        elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
        ylabel(elec_pair);
        ylim([0,1]);
    end
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    BipChanLabel = strrep(dataBipolar.label{Bip_chans(j-2)},'_',' ');
    titl = ['PLV significance Maintenance vs Encoding',BipChanLabel];
%     title(titl);
    suptitle(titl)
    end
    
    
    %% Visualization of PLV significance Maintenance vs Fixation
    
    PLV_Path_maint = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_f_1_30.mat'
    PLV_Path_fix = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_fixation_test_f_1_30.mat'
    PLV_Path_enc = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_test_f_1_30.mat'
    PLV_Pairs_maint_grid =  load(PLV_Path_maint);
    freqAxis  = PLV_Pairs_maint_grid.freqAxis
    PLV_Pairs_maint_grid = PLV_Pairs_maint_grid.PLV_Pairs;
    PLV_Pairs_fix_grid = load(PLV_Path_fix);
    PLV_Pairs_fix_grid = PLV_Pairs_fix_grid.PLV_Pairs;
    PLV_Pairs_enc_grid = load(PLV_Path_enc);
    PLV_Pairs_enc_grid = PLV_Pairs_enc_grid.PLV_Pairs;
    strPlotColors = {'r','g','b','k','c'};
    for j = 5%3:8
    fig = figure;
    ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
    
    for i = 1:nGridChans%18%58%
        axes(ha(plot_order(i)));
        iSS = 4;
        semilogx(freqAxis,abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS}),strPlotColors{4},'LineWidth',2);       
        hold on;
        semilogx(freqAxis,abs(PLV_Pairs_fix_grid{i+((j-1)*nGridChans),iSS}),strPlotColors{1},'LineWidth',2);    
        hold on;
%         semilogx(freqAxis,abs(PLV_Pairs_enc_grid{i+((j-1)*nGridChans),iSS}),strPlotColors{2},'LineWidth',2);

        Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{i+((j-1)*nGridChans)};
        Maint_Fix_Diff = abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS})-abs(PLV_Pairs_fix_grid{i+((j-1)*nGridChans),iSS});
        indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'k');
        if (any(indMaxBand==1))
            GridChans_sgnf_15_30(i) =1;
        else
            GridChans_sgnf_15_30(i) =0;
        end
        Grid_sgnf = [1 1 0 0 1 1 1 1;1 1 1 1 0 1 1 1;1 1 1 0 1 0 0 1;1 1 1 1 0 1 1 1;1 1 1 1 1 1 1 1 ;1 1 1 0 1 1 1 1;1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 0;];
        Fix_Maint_Diff = abs(PLV_Pairs_fix_grid{i+((j-1)*nGridChans),iSS})-abs(PLV_Pairs_maint_grid{i+((j-1)*nGridChans),iSS});
        indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'r');

        elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
        ylabel(elec_pair);
        ylim([0,1]);
        xlim([3 100])
    end
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    BipChanLabel = strrep(dataBipolar.label{Bip_chans(j-2)},'_',' ');
    titl = ['PLV significance Maintenance vs Fixation',BipChanLabel];
    title(titl);
    suptitle(titl)
    end
end
    
%% PSD significance in encoding vs maintenance
if PSD_significance
    
    %%
    nSS_Ran = [4,4];
    
    nNumberOfPermutations = 200;
    nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
    nTrialNumbersAll_AcrossSS = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
    nRandomTrialNumberList_AcrossTime_Pairs = cell(1,length(nChannelList));
    for iPair = 1:length(strChannelNameList)
        nPair = iPair; % nPairs_ToRun(iPair);
        nRandomTrialNumberList = zeros(nNumberOfPermutations,sum(nNumberOfSetSizes_ToCompare));
        for nRand = 1:nNumberOfPermutations
            cond = 1;
            while(cond)
                indRand = randperm(sum(nNumberOfSetSizes_ToCompare));
                cond = ~isempty(find((indRand-(1:length(indRand)))==0,1));
            end
            nRandomTrialNumberList(nRand,:) = nTrialNumbersAll_AcrossSS(indRand);
        end
        nRandomTrialNumberList = [nTrialNumbersAll_AcrossSS;nRandomTrialNumberList];
        nRandomTrialNumberList_AcrossTime_Pairs{nPair} = nRandomTrialNumberList;
    end
    
    %% Save Permutation numbers
    strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\'
    mkdir(strPaths.Results)
    strVariableFolder = [strPaths.Results,'Statistics\Rand_Trials_PSD\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms_PSD_maint_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')
    
    
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Enc_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    %% Extract 2 seconds of maintenance
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Maint_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
        
    %% Compute random PSD for Grid Channels
    strVariableFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PSD_pairs\Maintenance_vs_Encoding\';
    mkdir(strVariableFolder)
    %Grid Channels
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelList)','GL')));

    for iChannel = ind(1):ind(end)%1:length(nChannelList)
              
        nChannel = nChannelList(iChannel);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,2);
        iSS_AcrossSS = 1;
        iSS = nSS_Ran(iSS_AcrossSS);
        cfg = [];
        cfg.channel = nChannel;
        data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolar_Maint_SS{iSS});
        
        iSS_AcrossSS = 2;
        iSS = nSS_Ran(iSS_AcrossSS);
        cfg = [];
        cfg.channel = nChannel;
        data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolar_Enc_SS{iSS});
        
        
        %% Mix the trials
        Power_SS = {};
        for nRand = 1:nNumberOfPermutations+1
            
            %% Create datasets with randomized trials
            data_SinglePair_Rand_SS = data_SinglePair_SS;
            iTrialList = nRandomTrialNumberList_AcrossTime_Pairs{iChannel}(nRand,:);
            for iTrial = 1:length(iTrialList)
                nTrial = iTrialList(iTrial);
                if(iTrial<=nNumberOfSetSizes_ToCompare(1)) % assign values for set size 1
                    iSS_Rand = 1;
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                else % comes from set size 2
                    iSS_Rand = 2;
                    iTrial = iTrial-nNumberOfSetSizes_ToCompare(1);
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                end
                data_SinglePair_Rand_SS{iSS_Rand}.trial{iTrial} = ...
                    data_SinglePair_SS{iSS_Real}.trial{abs(nTrial)};
            end
            %% Power Spectrum Density for randomized trials
            for iSS = 1:2
                %% Frequency analysis
                cfg             = [];
                cfg.method      = 'mtmfft';
                cfg.taper       = 'dpss';
                cfg.output      = 'pow';
                cfg.foi         = 0.5:0.5:120;
                cfg.tapsmofrq   = 2;
                cfg.channel     = 1;
                cfg.keeptrials  = 'yes';
                Power_SS{iSS}{nRand}     = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
                
            end
            clear cfg
        end
        
        %% Store results for the channel pair
        strChannelLabel = Power_SS{1}{1}.label';
        freqAxis = Power_SS{1}{1}.freq;
        for iSS = 1:2
            for nRand = 1:nNumberOfPermutations+1
                Power_SS{iSS}{nRand} = mean(squeeze(Power_SS{iSS}{nRand}.powspctrm(:,1,:)));
            end
        end
        
        %% Save results
        iSS = 4;
        save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Maint_SS{iSS}.label{nChannel},'_mt_tsf_2_Maint_vs_Enc_PSD',num2str(2*iSS+2,'%.1d'),'.mat'],...
            'Power_SS','strChannelLabel','freqAxis','-v7.3')
        
        
        
        %% Save statistical summary
        %% Real values
        Real_Power_SS = {};
        nRand = 1;
        for iSS = 1:2
            Real_Power_SS{iSS} = Power_SS{iSS}{nRand};
        end
        
        %% Null distribution
        for iSS = 1:2
            RandDist_PSD_All{iSS} = [];
            for nRand = 1:length(Power_SS{iSS})-1
                RandDist_PSD_All{iSS}(nRand,:) = Power_SS{iSS}{nRand+1};
            end
        end
        
        PSD_Rand_Dist_SingleSS = abs(RandDist_PSD_All{1})-abs(RandDist_PSD_All{2});
        RandDist_PSD_All_Im = imag(RandDist_PSD_All{1})-imag(RandDist_PSD_All{2});
        
        %% Percentiles
        prcList = [1,5,50,95,99];
        RandPrc = {};
        RandPrc_Im = {};
        for iPrc = 1:length(prcList)
            val  = prctile(PSD_Rand_Dist_SingleSS,prcList(iPrc));
            RandPrc{prcList(iPrc)} = val;
            val  = prctile(RandDist_PSD_All_Im,prcList(iPrc));
            RandPrc_Im{prcList(iPrc)} = val;
        end
        
        %% Min max
        RandStats.Min = min(PSD_Rand_Dist_SingleSS);
        RandStats.Max = max(PSD_Rand_Dist_SingleSS);
        RandStats.Mean = mean(PSD_Rand_Dist_SingleSS);
        RandStats.Std = std(PSD_Rand_Dist_SingleSS);
        RandStats.NoPermutations = size(PSD_Rand_Dist_SingleSS,1);
        
        RandStats_Im.Min = min(RandDist_PSD_All_Im);
        RandStats_Im.Max = max(RandDist_PSD_All_Im);
        RandStats_Im.Mean = mean(RandDist_PSD_All_Im);
        RandStats_Im.Std = std(RandDist_PSD_All_Im);
        RandStats_Im.NoPermutations = size(RandDist_PSD_All_Im,1);
        
        %% Save stats
        strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\PSD_maint_vs_encoding\']
        mkdir(strPaths.stats)
        save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Maint_SS{iSS}.label{nChannel},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding_PSD','.mat'],...
        'Real_Power_SS','PSD_Rand_Dist_SingleSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strChannelLabel','freqAxis','-v7.3')  ;      
    
        clear Real_Power_SS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
        
    end
    
    %% Visualization for a single pair
    load('F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_Stats\PSD_maint_vs_encoding\Fixation_PSD.mat');
    
    %%% Save for publication !
    fig = figure;
%     ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
    for iChannel = 67%ind(1):ind(end)%67
        nChannel = iChannel;
        iSS = 4;
        strPaths.stats = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_Stats\PSD_maint_vs_encoding\';
        StatsVars = load([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Maint_SS{iSS}.label{nChannel},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding_PSD','.mat']);
        boolLinearScale = 0 ;
        
        if boolLinearScale
            Real_TP1 = StatsVars.Real_Power_SS{1};
            Real_TP2 = StatsVars.Real_Power_SS{2};
            ylimit = [-10 400];
            yShift = 15;
            
            
        else
            Real_TP1 = 10*log10(StatsVars.Real_Power_SS{1});
            Real_TP2 = 10*log10(StatsVars.Real_Power_SS{2});
            ylimit = [-15 40];
            yShift = -13.5;
        end
%         axes(ha(plot_order(iChannel-AHL_chans)))
        semilogx(StatsVars.freqAxis,Real_TP1,'k','LineWidth',3);
        hold on;
%         color_enc = [0.25, 0.25, 0.25]
        semilogx(StatsVars.freqAxis,Real_TP2,'g','LineWidth',3);
        hold on;
        semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(iChannel-AHL_chans,:)),'r','LineWidth',3);
        ylim([ylimit(1) ylimit(2)])
        xlim([3,100]);
        Diff1 = StatsVars.Real_Power_SS{1} - StatsVars.Real_Power_SS{2};
        Diff2 = StatsVars.Real_Power_SS{2} - StatsVars.Real_Power_SS{1};
        
        fAxis = StatsVars.freqAxis;
        color = 'm'
        color2 = [0.75, 0.75, 0];
        elecPair = strrep(dataBipolar_Maint_SS{4}.label{nChannel},'_',' ');
        ylabel(elecPair);
        
        indMaxBand = Difference_Bar(Diff1,StatsVars.RandPrc{95},fAxis,[ylimit(1) ylimit(2)],yShift,color);
        indMaxBand2 = Difference_Bar(Diff2,StatsVars.RandPrc{95},fAxis,[ylimit(1) ylimit(2)],yShift,color2);

    end
%     set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% ,strrep(strChannelNameList(nChannel),'_',' ')
    titl = ['PSD Maint-Encod significance Grid Channels']
%     suptitle(titl)
    
%     ylabel('Power 10*log10(\muV^2/Hz)')
%     xlabel('Frequency (Hz)')
%     
end
%% Set Size Comparison
if SetSizeComparison
     nSS_Ran = [4,1];
  
    nNumberOfPermutations = 200;
    nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
    nTrialNumbersAll_AcrossSS = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
    nRandomTrialNumberList_AcrossTime_Pairs = cell(1,size(nChannelPairs,1));
    for  iPair = 1:size(nChannelPairs,1)
        nPair = iPair; % nPairs_ToRun(iPair);
        nRandomTrialNumberList = zeros(nNumberOfPermutations,sum(nNumberOfSetSizes_ToCompare));
        for nRand = 1:nNumberOfPermutations
            cond = 1;
            while(cond)
                indRand = randperm(sum(nNumberOfSetSizes_ToCompare));
                cond = ~isempty(find((indRand-(1:length(indRand)))==0,1));
            end
            nRandomTrialNumberList(nRand,:) = nTrialNumbersAll_AcrossSS(indRand);
        end
        nRandomTrialNumberList = [nTrialNumbersAll_AcrossSS;nRandomTrialNumberList];
        nRandomTrialNumberList_AcrossTime_Pairs{nPair} = nRandomTrialNumberList;
    end
    
       %% Save Permutation numbers
    strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\'
    mkdir(strPaths.Results)
    strVariableFolder = [strPaths.Results,'Statistics\Rand_Trials_PLV\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms_maint_hi_vs_maint_lo.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')
 
    %% Extract 2 seconds of maintenance low
        dataBipolar_Maint_Low_SS = cell(1,size(Set_Sizes,1));
        for iSS = 1:size(Set_Sizes,1)
            cfg = [];
%             cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
            cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample]; % maintenance high vs maintenance low
            dataBipolar_Maint_Low_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
      
        end
    %% Pairs to run analysis for
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','-')));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_')));
    
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    strVariableFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_Hi_vs_Maintenance_Lo\';
    mkdir(strVariableFolder)
   %%     
     tStartRand = cputime;
    for iPair = 1:length(nPairs_ToRun)
        disp(iPair)

        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        

        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,2);
        Period = 1;
        iSS = nSS_Ran(Period);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{Period} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
        
        Period = 2;
        iSS = nSS_Ran(Period);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{Period} = ft_preprocessing(cfg,dataBipolar_Maint_Low_SS{iSS});
        
        
        %% Randomly mix trials between maintenance hi & maintenance lo

        for nRand = 1:nNumberOfPermutations+1
            
            %% Create datasets with randomized trials
            data_SinglePair_Rand_SS = data_SinglePair_SS;
            iTrialList = nRandomTrialNumberList_AcrossTime_Pairs{nPair}(nRand,:);
            for iTrial = 1:length(iTrialList)
                nTrial = iTrialList(iTrial);
                if(iTrial<=nNumberOfSetSizes_ToCompare(1)) % assign values for set size 1
                    iSS_Rand = 1;
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                else % comes from set size 2
                    iSS_Rand = 2;
                    iTrial = iTrial-nNumberOfSetSizes_ToCompare(1);
                    if(nTrial>0) % comes from set size 1
                        iSS_Real = 1;
                    else % comes from set size 2
                        iSS_Real = 2;
                    end
                end
                data_SinglePair_Rand_SS{iSS_Rand}.trial{iTrial} = ...
                    data_SinglePair_SS{iSS_Real}.trial{abs(nTrial)};
            end
            
            %% Phase coherence for randomized trials
            clear Freq_SS
            fprintf('Calculating %s with  %s right now. Elapsed pairs are %d..... \n',strChannelNameList{nChannel_1},strChannelNameList{nChannel_2},length(nPairs_ToRun)-iPair);

            for iSS = 1:2
                %% Frequency analysis
                cfg             = [];
                cfg.method      = 'mtmfft';
                cfg.taper       = 'dpss';
                cfg.output      = 'fourier';
                cfg.tapsmofrq   = 2;
                cfg.pad         = 2;
                cfg.foi         = 1:1:30;
                cfg.channel     = 1:2;
                Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
                
                %% PLV analysis
                cfg             = [];
                cfg.method      = 'plv';
                cfg.complex     = 'complex';
                PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            end
        end
          clear Freq_SS
           %% Store results for the channel pair
          strPairLabels = PLV_SS{iSS}{1}.label';
          freqAxis = PLV_SS{iSS}{1}.freq;
          for iSS = 1:2
              for nRand = 1:nNumberOfPermutations+1
                  PLV_AcrossPeriods{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
              end
          end
          clear PLV_SS
          
             %% Save results
          strPLVResultFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_Hi_vs_Maintenance_Lo\';
          mkdir(strPLVResultFolder);
          if Scalp_Hippocampus_coupling
             strVariableFolder = [strPLVResultFolder,'Scalp_Hippocampus\'];
             mkdir(strVariableFolder);
             save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_freq_1_1_30_Maint_Hi_vs_Maint_Lo_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
              'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
          elseif Hippocampus_Grid_Coupling
                  strVariableFolder = [strPLVResultFolder,'Grid_Hippocampus\'];
                  mkdir(strVariableFolder)
                  save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Maint_Hi_vs_Maint_Lo_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
                      'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
             
         elseif Scalp_ECoG_coupling  
                  strVariableFolder = [strPLVResultFolder,'Scalp_Grid\'];
                  mkdir(strVariableFolder)
                  save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_Maint_Hi_vs_Maint_Lo_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
                      'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
          end
               
          
          %% Real values
          nRand = 1;
          for iSS = 1:2
              Real_PLV_AcrossPeriods{iSS} = PLV_AcrossPeriods{iSS}{nRand};
          end
          %% Null distribution
          for iSS = 1:2
              RandDist_PLV_All_AcrossPeriods{iSS} = [];
              for nRand = 1:length(PLV_AcrossPeriods{iSS})-1
                  RandDist_PLV_All_AcrossPeriods{iSS}(nRand,:) = PLV_AcrossPeriods{iSS}{nRand+1};
              end
          end
          clear PLV_AcrossPeriods
          
          PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossPeriods{1})-abs(RandDist_PLV_All_AcrossPeriods{2});
          RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossPeriods{1})-imag(RandDist_PLV_All_AcrossPeriods{2});
         
          %% Percentiles
          prcList = [1,5,50,95,99];
          RandPrc = {};
          RandPrc_Im = {};
          for iPrc = 1:length(prcList)
              val  = prctile(PLV_Rand_Dist_SingleSS,prcList(iPrc));
              RandPrc{prcList(iPrc)} = val;
              val  = prctile(RandDist_PLV_All_Im,prcList(iPrc));
              RandPrc_Im{prcList(iPrc)} = val;
          end
          
          %% Min max
          RandStats.Min = min(PLV_Rand_Dist_SingleSS);
          RandStats.Max = max(PLV_Rand_Dist_SingleSS);
          RandStats.Mean = mean(PLV_Rand_Dist_SingleSS);
          RandStats.Std = std(PLV_Rand_Dist_SingleSS);
          RandStats.NoPermutations = size(PLV_Rand_Dist_SingleSS,1);
          
          RandStats_Im.Min = min(RandDist_PLV_All_Im);
          RandStats_Im.Max = max(RandDist_PLV_All_Im);
          RandStats_Im.Mean = mean(RandDist_PLV_All_Im);
          RandStats_Im.Std = std(RandDist_PLV_All_Im);
          RandStats_Im.NoPermutations = size(RandDist_PLV_All_Im,1);
          
          %% Save results
          if Scalp_Hippocampus_coupling
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Scalp_Hippocampus\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Maint_Lo','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
             
          elseif Hippocampus_Grid_Coupling
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Grid_Hippocampus\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_Hi_vs_Maint_Lo','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
                  
           
          elseif Scalp_ECoG_coupling
            
                  strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Scalp_Grid\']
                  mkdir(strPaths.stats)
                  save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Scalp',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_Hi_vs_Maint_Lo','.mat'],...
                      'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                  clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
                  
          end
          
    end
    tStopRand = cputime-tStartRand;
             
     %% Load the results 
   if Scalp_Hippocampus_coupling
       strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Scalp_Hippocampus\']
      
   elseif Hippocampus_Grid_Coupling
       strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Grid_Hippocampus\']
      
   elseif Scalp_ECoG_coupling
       strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_Hi_vs_Maintenance_Lo\Scalp_Grid\']
           
           
   end
    PLV_Rand_MaintHi_MaintLo_Variables_Pair_SS = cell(1,size(Set_Sizes,1));
    for iSS = 1:2%size(Set_Sizes,1)
        PLV_Rand_MaintHi_MaintLo_Variables_Pair_SS{nSS_Ran(iSS)} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            if Scalp_Hippocampus_coupling
                    strVariablePath =  [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Maint_Lo','.mat'];
            elseif Hippocampus_Grid_Coupling
                    strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_Hi_vs_Maint_Lo','.mat'];
               
            elseif Scalp_ECoG_coupling
                    strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Scalp',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_Hi_vs_Maint_Lo','.mat'];
               
            end

%             [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Encoding','.mat'];
            try
                Vars = load(strVariablePath);
                PLV_Rand_MaintHi_MaintLo_Variables_Pair_SS{nSS_Ran(iSS)}{nPair} = Vars;
                clear Vars
            catch
            end
        end
    end
    %% Load PLV_Pairs
   if Scalp_Hippocampus_coupling
       
           PLV_Path_maint_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_maintenance_cleanEEG.mat';
           PLV_Path_fix_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_fixation_cleanEEG.mat';
            load(PLV_Path_maint_sclp);
            load(PLV_Path_fix_sclp);
    
   elseif Scalp_ECoG_coupling
           PLV_Path_maint_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_maintenance_clean_EEG_f_1_30_all_scalp_chans.mat';
           PLV_Path_fix_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\Patient_42_Scalp-Grid_PLV_Pairs_fixation_clean_EEG_f_1_30_all_scalp_chans.mat';
           load(PLV_Path_maint_sclp);
           load(PLV_Path_fix_sclp);
       
   end
    %% Scalp Hippocampus Maint Hi vs Maint Lo
    freqBand_topo = [6 7];
    [~,indFreq1] = min(abs(freqAxis_sclp-freqBand_topo(1)));
    [~,indFreq2] = min(abs(freqAxis_sclp-freqBand_topo(2)));
    iSSLow = 1;
    saveFlag = 0;
%     figure;
    figLayout = Get_Figure_Scalp_Layout_Information();
     for nPair = 1:size(nPairs_ToRun,2)%61
        iSS = 4;
        if Scalp_Hippocampus_coupling
            Vars = PLV_Rand_MaintHi_MaintLo_Variables_Pair_SS{iSSLow}{nPair};
            MaintHi_vs_Lo_Diff = abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS})-abs(PLV_Pairs_Scalp_maint_HC{nPair,iSSLow});
            datavector(nPair) = median(abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS}(indFreq1:indFreq2)))
        end
        if  nChannelPairs(nPair,2) == 1 %For every different hippocampal channel paired with Scalp channels
            fig = figure
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            set(gcf,'color','white')
            ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
            
        end
        Scalp_channel = strrep(dataBipolarSclp.label{nChannelPairs(nPair,2)},'_',' ');
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_maint(nPair) = figLayout{indSubplot,2};
        axes(ha(nSubplot_maint(nPair)));
        nChannel1 = strrep( strChannelNameList(nChannelPairs(nPair,1)),'_',' ');
        nChannel2 = strrep( strChannelNameList(nChannelPairs(nPair,2)),'_',' ');
        
        %Visualization
        %     plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossPeriods{1})-abs(Vars.Real_PLV_AcrossPeriods{2}),'LineWidth',2)
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS}),'k','LineWidth',3)
        ylim([0,1]);

        hold on
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSSLow}),'b','LineWidth',3)
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_HC{nPair,iSS}),'r','LineWidth',3)
        hold on;
        indMaxBandMaintHivsLo = Difference_Bar(MaintHi_vs_Lo_Diff,Vars.RandPrc{95},freqAxis_sclp,[0 1], 0.05,'b');
        ylim([0,1])
        strTitle = strcat('Maint High/Maint Low Significance Pair: ',nChannel1,'-', ...
            'Scalp',' Set Size [',num2str(Set_Sizes{iSS}),' ]');
        ylabel(nChannel2);
        xlim([3 100])
       Fig_dir = [Path_for_PLV_sgnf,'MaintHi vs MaintLow\','Patient_42_MaintHi_MaintLo_diff\'];
        mkdir(Fig_dir);
%         Figure_Name = strcat(Fig_dir,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),strjoin(nChannel1),'_','Scalp','_MaintHi_MaintLo');
        
        
        
        if ~mod(nPair,size(dataBipolarSclp.label,1)) %% When a hippocampal is paired with every scalp channel 
            suptitle(strTitle)
%             set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
%             set(ha(1:size(ha,1)),'FontSize',18);
            for i = 1:size(ha,1)
                
                if ~ismember(i,nSubplot_maint)
                    set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
                end
                
            end
            if saveFlag
                saveas(fig,Figure_Name,'png');
                saveas(fig,Figure_Name,'fig');
            end
        end
        
     end
     
     
% datavector = [0.57 0.56  0.50 0.47 0.27 0.33 0.35 0.4 0.51 0.35 0.57 0.37 0.49 0.51 0.51 0.47 0.33 0.27 0.27 0.54 0.49 0.58 0.56]'     
     dataBipolarSclp.label = strrep(dataBipolarSclp.label,'_',' ');
%      Scalp_topoplot(strPaths.Toolboxes,freqBand_topo,freqAxis_sclp,datavector(47:69),dataBipolarSclp,4,nChannelPairs,1)
%      Scalp_topoplot(strPaths.Toolboxes,freqBand_topo,freqAxis_sclp,datavector(nPair),dataBipolarSclp,4,nChannelPairs,1)
   
        
    
        
%     %%
%     % Load PLV Pairs for Maint and Enc
%         load(PLV_Path_maint);
%         load(PLV_Path_enc);
%         
%         if Scalp_Hippocampus_coupling
%             nScalpChans = size(dataBipolarSclp.label,1);
%             nHippChans = size(nChannelPairs_Sclp_HC,1)/nScalpChans;
%             iSS_ToPLot = 4;
%             iSS_to_compare = 1;
%             figLayout = Get_Figure_Scalp_Layout_Information;
%             for j = 1:nHippChans
%                 
%                 fig = figure;
%                 set(fig,'units','normalized','outerposition',[0 0 1 1])
%                 set(gcf,'color','white')
%                 ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01]);
%                 
%                 for i = 1:nScalpChans
%                     strScalpChannel = strrep(dataBipolarSclp.label{i},'_',' ');
%                     strHippChannel = strrep(dataBipolar.label{nChannelPairs_Sclp_HC(nScalpChans*j,1)},'_',' ');
%                     indSubplot = find(strcmpi(figLayout(:,1),strScalpChannel));
%                     nSubplot(i) = figLayout{indSubplot,2};
%                     nPair = i+((j-1)*nScalpChans);
%                     axes(ha(nSubplot(i)));
%                     plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS_ToPLot}),'r','LineWidth',2);
%                     hold on;
%                     plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{nPair,iSS_ToPLot}),'g','LineWidth',2);
%                     hold on;
%                     plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPair,iSS_to_compare}),'b','LineWidth',2); %
%                     ylabel(strScalpChannel);
%                     ylim([0 1])
%                     
%                 end
%                 set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
%                 suptitle(['Scalp-Hippocampus ',strHippChannel,' Maintenance vs Encoding and Set Size Comparison, Set Size ','[', num2str(Set_Sizes{iSS_ToPLot}),'] vs [', num2str(Set_Sizes{iSS_to_compare}),']'])
% 
%                 for i = 1:size(ha,1)
%                     if ~ismember(i,nSubplot)
%                         set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
%                     end
%                     
%                 end
%                 Fig_dir_SS = [Path_for_PLV_sgnf,'Patient_42_SetSize_Comparison\'];
%                 mkdir(Fig_dir_SS);
%                 Figure_Name = strcat(Fig_dir_SS,'Patient_',num2str(42),'_Pair_',strHippChannel,'_SS_comparison ', ...
%                     '[ ',num2str(Set_Sizes{iSS_ToPLot}),' ] -','[ ',num2str(Set_Sizes{iSS_to_compare}),']');
%                 
%                 saveas(fig,Figure_Name,'png');
%                 saveas(fig,Figure_Name,'fig');
%             end
%         end
end