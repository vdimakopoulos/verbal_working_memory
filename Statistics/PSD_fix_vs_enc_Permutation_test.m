%% Close all figures, clear variables and command window
close all
clear
clc

%% Paths
Drive_Letter = 'F:\';
strPaths.Main = [Drive_Letter,'Vasileios\'];
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
strPaths.GeneralFunctions = [Drive_Letter,'Vasileios\Task Analysis\Code\'];%'F:\Vasileios\Task Analysis\Code\';
strPaths.Data = [strPaths.Main, 'Task Analysis\Data\'];
strPaths.ExtractedData = [strPaths.Main, 'Task Analysis\Extracted_Data Sternberg\'];
strPaths.Statistics = [strPaths.Main, 'Task Analysis\Code\Statistics\'];
strPaths.ChanLoc = [strPaths.Main, 'Task Analysis\Code\Channel Localization\'];
strPaths.EEGLAB_Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\EEGLAB Subfunctions\'];
strPaths.Subfunctions = [strPaths.Main,'Task Analysis\Code\Subfunctions\']
% Results
strPaths.Results = [strPaths.Main,'Task Analysis\Analysis Results\'];

% FieldTrip toolbox
% strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20191126\';
strPaths.Toolboxes.FieldTrip            = [Drive_Letter,'Vasileios\Toolboxes\fieldtrip-20200315\'];

% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = [Drive_Letter,'Vasileios\Toolboxes\eeglab14_1_1b\'];

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


% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Flags for coupling types and other flags
Hippocampus_Grid_Coupling           =   1;
Hippocampus_Coupling                =   0;
Scalp_Hippocampus_coupling          =   0;
Scalp_ECoG_coupling                 =   0;
Task_Period = 2; %Values 1,2,3 for Encoding,Maintenance,Retrieval period correspondigly
PLV_significance = 0;
PLV_significance_maint_vs_encod = 0;
PLV_significance_maint_vs_fix = 0;
PSD_significance = 1;
SetSizeComparison = 0;
Gamma_Freq = 0;


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


%% Merge sessions
cfg = [];
dataBipolar = ft_appenddata(cfg,data_s1,data_s2);
dataBipolarScalp = ft_appenddata(cfg,data_scalp_s1,data_scalp_s2);


[dataBipolar_test,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable);

%% Reref based on the average of the signals
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =  [3 4] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

%% Apply montage
clear montage;
montage.labelold = dataBipolar.label;
num_bipolar_chans=6;
num_reference_chans=size(montage.labelold,1);
%prepare the montage matrix
montage_matrix = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar = [1 2 1 3 2 3 73 74 73 75 74 75 ] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
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
cfg.refmethod  = 'bipolar'
cfg.refchannel = [82 85]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));

cfg = [];
cfg.resamplefs = 500 ;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);

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
%% channel pairs
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;
chans_to_be_bipolar = [1 2 1 3 2 3 73 74 73 75 74 75 ];
if PLV_depth_electrodes == 1
    
    nChannelPairs = [81 83;81 84; 82 83; 82 84];
    
else
    nGridChans = 56;
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
    nChannelPairs = [];
    nChannelPairs = [[ones(nGridChans,1),(1:nGridChans)'+AHL_chans];[nGridChans+AHL_chans+ones(nGridChans,1),(1:nGridChans)'+AHL_chans]];
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChans,1),(1:nGridChans)'+ AHL_chans]];
        %         nChannelPairs = [nChannelPairs;[round(i/2)+72+zeros(nGridChans,1),(1:nGridChans)'+ 8]];
        
    end
    
    
end

%% Number of trials for each set size
nNumberOfTrialsSetSizes = zeros(1,size(dataBipolar_SS,2));
for iSS = 1:size(dataBipolar_SS,2)
    nNumberOfTrialsSetSizes(iSS) = length(dataBipolar_SS{iSS}.trial);
end

%% Permute trials
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
save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','_Randomized_Trial_Numbers_200_perms_PSD_fix_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')

%% Extract 2 seconds of encoding latency
nSet_Size = size(dataBipolar_SS,2);
dataBipolar_Ret_SS = {};
for iSS = 1:nSet_Size
    cfg = [];
    cfg.latency = [-3.5,-3-1/dataBipolar_SS{iSS}.fsample];
    dataBipolar_Enc_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end
%% Extract 1 seconds of fixation
for iSS = 1:size(Set_Sizes,1)
    cfg = [];
    cfg.latency = [-5.5,-5-1/dataBipolar_SS{iSS}.fsample];
    dataBipolar_Fix_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    cfg =[];
%     cfg.resamplefs = 1000;
%     dataBipolar_Fix_SS{iSS} = ft_resampledata(cfg,dataBipolar_Fix_SS{iSS});
    
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
    data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolar_Fix_SS{iSS});
    
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
            cfg.foi         = 1:100;
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1;
            cfg.keeptrials  = 'yes';
            Power_SS{iSS}{nRand}     = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
        end
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
    iSS = 5;
    save([strVariableFolder,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Fix_SS{iSS}.label{nChannel},'_mt_tsf_2_Fix_vs_Enc_PSD',num2str(2*iSS+2,'%.1d'),'.mat'],...
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
    strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\PSD_fix_vs_encoding\']
    mkdir(strPaths.stats)
    save([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Fix_SS{iSS}.label{nChannel},'_mt_tsf_2_Stats_Summary_Across_Periods_Fix_vs_Encoding_PSD','.mat'],...
        'Real_Power_SS','RandDist_PSD_All','PSD_Rand_Dist_SingleSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strChannelLabel','freqAxis','-v7.3')  ;
    
    clear Real_Power_SS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end


%% Visualization for one pair

%%% Save for publication !
fig = figure;
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for iChannel = ind(1):ind(end)
    nChannel = iChannel;
    iSS = 4;
    strPaths.stats = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\Statistics\Rand_Stats\PSD_fix_vs_encoding\';
    StatsVars = load([strPaths.stats,'Patient_',num2str(42),'_Pair_',num2str(nChannel,'%.4d'),'_',dataBipolar_Fix_SS{iSS}.label{nChannel},'_mt_tsf_2_Stats_Summary_Across_Periods_Fix_vs_Encoding_PSD','.mat']);
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
    axes(ha(plot_order(iChannel-8)))
    semilogx(StatsVars.freqAxis,Real_TP1,'k','LineWidth',3);
    hold on;
    %         color_enc = [0.25, 0.25, 0.25]
    semilogx(StatsVars.freqAxis,Real_TP2,'g','LineWidth',3);
    hold on;
    %         semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(iChannel-AHL_chans,:)),'r','LineWidth',3);
    ylim([ylimit(1) ylimit(2)])
    xlim([3,100]);
    Diff1 = StatsVars.Real_Power_SS{1} - StatsVars.Real_Power_SS{2};
    Diff2 = StatsVars.Real_Power_SS{2} - StatsVars.Real_Power_SS{1};
    
    fAxis = StatsVars.freqAxis;
    color = 'm'
    color2 = [0.75, 0.75, 0];
    elecPair = strrep(dataBipolar_Fix_SS{4}.label{nChannel},'_',' ');
    ylabel(elecPair);
    
    indMaxBand = Difference_Bar(Diff1,StatsVars.RandPrc{95},fAxis,[ylimit(1) ylimit(2)],yShift,color);
    indMaxBand2 = Difference_Bar(Diff2,StatsVars.RandPrc{95},fAxis,[ylimit(1) ylimit(2)],yShift,color2);
    
    surrogate = (Diff2) -  (StatsVars.RandPrc{95})%StatsVars.RandPrc{95};
    FreqBand = [60,100];
    freqAxis = fAxis
    [~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
    [~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
    PSD_temp = surrogate([indFreq1:indFreq2]);
    PSD_In_Band_SS(iChannel-8) = max(zscore(PSD_temp));

end
PSD_In_Band_ToPlot = reshape(PSD_In_Band_SS,8,8);
%%
figure;
imagesc(PSD_In_Band_ToPlot,[2 3])

ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'

colorbar
colormap(bluewhitered)
set(gca,'Fontsize',16)
