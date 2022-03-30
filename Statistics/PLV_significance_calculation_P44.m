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

%% Flags for coupling types and other flags
Hippocampus_Grid_Coupling           =   0;
Hippocampus_Coupling                =   0;
Scalp_Hippocampus_coupling          =   1;
Scalp_ECoG_coupling                 =   0;
Task_Period = 2; %Values 1,2,3 for Encoding,Maintenance,Retrieval period correspondigly
PLV_significance = 0;
PLV_significance_maint_vs_encod = 1;
PLV_significance_maint_vs_fix = 1;
PSD_significance = 0;
SetSizeComparison = 0;
Gamma_Freq = 0;

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

%% Load data for 4 sessions - Macro Data
strMacroDataDir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\44 MJ\';
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
 data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end


%% Load data for 4 sessions - Scalp Data
strScalpDataDir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\44 MJ\';
cd (strScalpDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
 data_Scalp_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable_Scalp = [];
for i=1:length(data_Scalp_ses)
    TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
end

%% Merge sessions
dataBipolar = data_ses(1).data;
for nSes = 2:length(data_ses)
    cfg = [];
    dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
    
end

dataBipolarScalp = data_Scalp_ses(1).data;
for nSes = 2:length(data_Scalp_ses)
       cfg = [];
    dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
    
end
cd(strPaths.Project)

%% Rereferencing
cfg               =   [];
cfg.reref         =  'yes' 
cfg.refchannel    =  [4] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 24;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 33 34 33 35 34 35 41 42 41 43 42 43 ...
    49 50 49 51 50 51 57 58 57 59 58 59] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end
%Append the bipolar channel labels to the reference channels labels
for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
montage.tra             = montage_matrix;
% dataBipolar = ft_apply_montage(dataBipolar,montage);

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod  = 'bipolar'
cfg.refchannel = [66:1:88]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);

 %% ICA  
 cfg = [];
 cfg.method       = 'runica'
 cfg.channel      = {'all' '-Submm','-Submp','-EKGm','-EKGp'};
 cfg.numcomponent = 'all';
 cfg.demean       = 'yes';
 cfg.feedback     =  'text'
 IC_components = ft_componentanalysis(cfg,dataBipolarScalp)
 
 %% Reject IC components
 cfg = [];
 cfg.component = [1 2 5 6 9 12 13 14 15 16 18 20 21];
 cfg.demean = 'yes';
 dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);
 
 cfg =[];
 cfg.channels = {'all' '-Submm','-Submp','-EKGm','-EKGp'};
 dataBipolarScalp_Reconstructed = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
 

 %% re-ref EEG
 cfg               = []
 cfg.channel       = {'all'}
 cfg.reref         =  'yes'
 cfg.refchannel    =   [1, 2]
 cfg.refmethod     =  'avg'
 dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp_Reconstructed);
 
 %% Merge Scalp and macro data
cfg = [];
cfg.keepsampleinfo = 'no'
dataBipolar_M = ft_appenddata(cfg,dataBipolarScalp_Reconstructed,dataBipolar);

%% Select only correct trials
[dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar_M,TrialInformationTable_Scalp);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable);

strChannelNameList = dataBipolar_SS{1}.label;


%% Extract 2 seconds of encoding
if Task_Period == 1
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
%% Extract 2 seconds of maintenance
elseif Task_Period == 2
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

%% Channel Pairs
nScalpChans = numel(dataBipolarScalp_Reconstructed.label);

AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
AHR_chans = length(find(contains(montage.labelold, 'AHR')==1));
PHR_chans = length(find(contains(montage.labelold, 'PHR')==1));
AL_chans = length(find(contains(montage.labelold, 'AL')==1));
AR_chans = length(find(contains(montage.labelold, 'AR')==1));
ECL_chans = length(find(contains(montage.labelold, 'ECL')==1));
ECR_chans = length(find(contains(montage.labelold, 'ECL')==1));
Depth_electrodes = AHL_chans+PHL_chans+AHR_chans+PHR_chans+AL_chans+AR_chans+ECL_chans+ECR_chans;

if Scalp_Hippocampus_coupling == 1
    
    nChannelPairs = [];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+Depth_electrodes+nScalpChans+zeros(nScalpChans,1),(1:nScalpChans)']];
    end
end

%% Number of trials for each set size
nNumberOfTrialsSetSizes = zeros(1,size(dataBipolar_SS,2));
for iSS = 1:size(dataBipolar_SS,2)
    nNumberOfTrialsSetSizes(iSS) = length(dataBipolar_SS{iSS}.trial);
end
%%
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
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms_maint_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')
    
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
        if isempty(find(contains(strChannelNameList{1},' ')))
            for i =1:nScalpChans
                strChannelNameList(i) = strcat({' '}, strChannelNameList{i});
            end
        end
        ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','PHR1-PHR2')));
        ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))',' ')));
        
        nPairs_ToRun = intersect(nPairs_ToRun,ind);
        nPairs_ToRun = intersect(nPairs_ToRun,ind2);
        
        strVariableFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_vs_Fixation\';
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
                cfg.foi         = [1:1:30 30:5:100];
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
        strPLVResultFolder = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_PLV_pairs\Maintenance_vs_Fixation\';
        
        strVariableFolder = [strPLVResultFolder,'Scalp_Hippocampus\'];
        mkdir(strVariableFolder);
        save([strVariableFolder,'Patient_',num2str(44),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_freq_1_1_30_Maint_vs_Fix_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
            'PLV_AcrossPeriods','strPairLabels','freqAxis','-v7.3')
        
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
        
        strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Hippocampus\Maint_vs_Fix\']
        mkdir(strPaths.stats)
        save([strPaths.stats,'Patient_',num2str(44),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation','.mat'],...
            'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
        clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    end
    clear PLV_SS;
    tStopRand = cputime-tStartRand;
    
    %% Load the results
    if Scalp_Hippocampus_coupling
        strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Scalp_Hippocampus\Maint_vs_Fix\']
    end
    PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,size(Set_Sizes,1));
    for iSS = 1:2%size(Set_Sizes,1)
        PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            if Scalp_Hippocampus_coupling
                if PLV_significance_maint_vs_fix
                    strVariablePath =  [strPaths.stats,'Patient_',num2str(44),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation','.mat']
                    
                end
            end
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
        
        PLV_Path_maint_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_44_Scalp-Hippocampus_PLV_Pairs_maintenance_f_1_1_30_5_100.mat';
        PLV_Path_fix_sclp = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_44_Scalp-Hippocampus_PLV_Pairs_fixation_f_1_1_30_5_100.mat';
        load(PLV_Path_maint_sclp);
        load(PLV_Path_fix_sclp);
        
    end
    
    
    %% Plot Scalp Depth PLV for one or more bipolar channels with significance bars
    Bipolar_Channels = unique(nChannelPairs_Sclp_depth(:,1));
    Bip_chans_to_plot = Bipolar_Channels(22); %%PHR1-2(22)
    iSS = 4; % Set Size [6 8]
    
    for i = 1:numel(Bip_chans_to_plot)
        nFig = figure
        indBipChan = find(nChannelPairs_Sclp_depth==Bip_chans_to_plot(i))
        ha =  tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
        set(gcf,'color','white')
        figLayout = Get_Figure_Scalp_Layout_Information();
        for j = 1:numel(indBipChan)
            Scalp_channel = strcat({' '},dataBipolar_appended{iSS}.label{j});
            indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
            nSubplot(j) = figLayout{indSubplot,2};
            axes(ha(nSubplot(j)));
            nPair = indBipChan(j);
            Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair};
            Maint_Fix_Diff = abs(PLV_Pairs_Scalp_maint_depth{nPair,iSS})-abs(PLV_Pairs_Scalp_fix_depth{nPair,iSS});
            semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_depth{indBipChan(j),iSS}),'k','LineWidth',3);
            hold on;
            semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_depth{indBipChan(j),iSS}),'r','LineWidth',3);
            ylim([0 1]);
            xlim([4 100]);
            indMaxBandMaintFix = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis_sclp,[0 1], 0.05,'r');
            ylabel(Scalp_channel);
            hold off;
        end
        suptitle(sprintf('Depth Electrode %s - Scalp PLV',dataBipolar_appended{iSS}.label{Bip_chans_to_plot(i)}))
        for k = 1:size(ha,1)
            
            if ~ismember(k,nSubplot)
                set(ha(k),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
            end
            
        end
    end
    
end