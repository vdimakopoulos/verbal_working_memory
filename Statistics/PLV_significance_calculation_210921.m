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

%% Reref to different sources
strChannelNameList = dataBipolar.label;
AhippChans = find(contains(strChannelNameList,'AH'));
PhippChans = find(contains(strChannelNameList,'PH'));
GridChans = find(contains(strChannelNameList,'GL'));
dataRerefAhipp = reref_toSeparate_Chan(dataBipolar, 4, 'all','avg');
dataRerefGrid = reref_toSeparate_Chan(dataBipolar, 5, 'all','avg');
dataRerefPhipp = reref_toSeparate_Chan(dataBipolar, 4, 'all','avg');

%selection of channels
cfg = [];
cfg.channel = AhippChans;
dataRerefAhipp = ft_selectdata(cfg,dataRerefAhipp);
cfg = [];
cfg.channel = GridChans;
dataRerefGrid = ft_selectdata(cfg,dataRerefGrid);
cfg = [];
cfg.channel = PhippChans;
dataRerefPhipp = ft_selectdata(cfg,dataRerefPhipp);

% append them again
cfg = [];
dataBipolar = ft_appenddata(cfg,dataRerefAhipp,dataRerefGrid);
dataBipolar = ft_appenddata(cfg,dataBipolar,dataRerefPhipp);

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);


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

chans_to_be_bipolar = [1:3 73:75]; %% Comment this if you have bipolar hipp chans
nChannelPairs = []
for i = 1:length(chans_to_be_bipolar)
    nChannelPairs = [nChannelPairs; [[repelem(chans_to_be_bipolar(i),length(GridChans))]' [GridChans]]]
end
strChannelNameList = dataBipolar.label;

%% Number of trials for each set size
nNumberOfTrialsSetSizes = zeros(1,size(dataBipolar_SS,2));
for iSS = 1:size(dataBipolar_SS,2)
    nNumberOfTrialsSetSizes(iSS) = length(dataBipolar_SS{iSS}.trial);
end


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
save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Hipp-Grid','_Reref_Randomized_Trial_Numbers_200_perms_maint_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')

%% Extract 1s of fixation
dataBipolar_Enc_SS = cell(1,size(Set_Sizes,1));
for iSS = 1:size(Set_Sizes,1)
    cfg = [];
    cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample]; % fixation vs maintenance
    dataBipolar_Enc_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    cfg =[];
    cfg.resamplefs = 1000;
    dataBipolar_Enc_SS{iSS} = ft_resampledata(cfg,dataBipolar_Enc_SS{iSS});
    
end

%% Pairs to run analysis for
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','_AHL2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_GL_C2')));

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
            %                     cfg.foi         = 0.5:0.5:120;
            cfg.foi         = [1:29 30:5:100];
            
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
    strVariableFolder = [strPLVResultFolder,'Grid_Hippocampus\Fix_vs_maint\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Reref_Fix_vs_maint_mt_tsf_2_Set_Size_freq_1_30_',num2str(2*iSS+2,'%.1d'),'.mat'],...
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
    strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
    mkdir(strPaths.stats)
    save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','_reref.mat'],...
        'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
    clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
end
clear PLV_SS;
tStopRand = cputime-tStartRand;
%% Load the results
strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,size(Set_Sizes,1));
for iSS = 1:2%size(Set_Sizes,1)
    PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)} = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','_reref.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPair} = Vars;
            clear Vars
        catch
        end
    end
end

%% plot the results
PLV_maint = abs(PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPairs_ToRun}.Real_PLV_AcrossPeriods{1});
PLV_fix = abs(PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPairs_ToRun}.Real_PLV_AcrossPeriods{2});

Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPairs_ToRun};
freq = Vars.freqAxis;
figure;
semilogx(freq,PLV_fix,'LineWidth',3,'Color','k')
hold on;
semilogx(freq,PLV_maint,'LineWidth',3,'Color','r')
ylim([-0.02 1])
xlim([4 100])
xlabel('Frequency (Hz)')
ylabel('PLV')
indMaintFixbar = Difference_Bar((PLV_maint-PLV_fix),Vars.RandPrc{95},freq,[0 1], -0.005,'r');
set(gca,'box','off','FontSize',16,'TickDir','out','XTick',[4 10 20 100],'XTickLabel',[4 10 20 100],...
    'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1])