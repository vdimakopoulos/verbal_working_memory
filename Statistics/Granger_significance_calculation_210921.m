%%
clear;
close all;
clc;

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
strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\';
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
Task_Period                         = 1; %Values 1,2,3 for Encoding,Maintenance,Retrieval period correspondigly
Granger_significance   = 1;
Granger_significance_maint_vs_encod = 0;

%% Merge sessions
cfg = [];
dataBipolar = ft_appenddata(cfg,data_s1,data_s2);
dataBipolarScalp = ft_appenddata(cfg,data_scalp_s1,data_scalp_s2);

[dataBipolar_test,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable);

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
%%
if Granger_significance
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
    strVariableFolder = [strPaths.Results,'Statistics\Rand_Trials_Granger\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Grid-Hippocampus','_Randomized_Trial_Numbers_200_perms_reref.mat'],'nRandomTrialNumberList_Pair_SS','nNumberOfPermutations','-v7.3')
    %% Channel Pairs
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    ind =find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','PHR1')));
    %find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','_AHL2')));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','TLSL1')));
    %find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_GL_C2')));
    
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
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
                cfg.foi         = [4:1:100];
                cfg.tapsmofrq   = 2;
                cfg.channel     = 1:2;
                Freq    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS}{nRand});
                
                cfgFreq = cfg;
                
                %% Granger analysis
                grangercfg = [];
                grangercfg.method  = 'granger';
                grangercfg.channelcmb = {char(strChannelNameList(nChannel_1)), char(strChannelNameList(nChannel_2))};%{'_AHL2-_AHL3','_GL_C2'}; %{'_GL_C2','_AHL2-_AHL3'}
                grangercfg.granger.conditional = 'no';
                % grangercfg.complex = 'yes';
                grangercfg.granger.sfmethod = 'bivariate';
                gdata_Pairs_ss{iSS}{nRand}    = ft_connectivityanalysis(grangercfg, Freq);

                
            end
        end
        clear Freq
        fprintf('\nChannel pair %d/%d complete\n',iPair,length(nPairs_ToRun))
        
        %% Store results for the channel pair
        iSS = iSS_ToCompute;
        strPairLabels = gdata_Pairs_ss{iSS}{1}.labelcmb';
        freqAxis = gdata_Pairs_ss{iSS}{1}.freq;
        for nRand = 1:length(gdata_Pairs_ss{iSS})
            gdata_Pairs_ss{iSS}{nRand} = gdata_Pairs_ss{iSS}{nRand}.grangerspctrm;
        end
%         TimeInterval = [-2,-1/dataBipolar_Ret_SS{iSS}.fsample];
        
        %% Convert cell to double
        for iSS = 1:size(dataBipolar_Ret_SS,2)
            try
                gdata_Pairs_ss{iSS} = cell2mat(gdata_Pairs_ss{iSS}');
            catch
                gdata_Pairs_ss{iSS} = [];
            end
        end
        
        %% Save values
        iSS = iSS_ToCompute;
        strPaths.var = [strPaths.Results,'Statistics\Rand_Granger_pairs\']
        mkdir(strPaths.var)
        save([strPaths.var ,'Patient_',num2str(42,'%.2d'),'_Rand_RerefGranger_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_ss',num2str(iSS),'_reref.mat'],...
            'gdata_Pairs_ss','strPairLabels','freqAxis','cfgFreq','-v7.3')
        
        %% Real values
        nRand = 1;
        Real_Granger = gdata_Pairs_ss{iSS}(nRand,:);
        Real_Granger2 = gdata_Pairs_ss{iSS}(2,:);

        %% Null distribution
        Granger_Rand_Dist_SingleSS = gdata_Pairs_ss{iSS}(202:end,:)- gdata_Pairs_ss{iSS}(1:201,:); %??
%         for iSS = 1:2
%             Granger_Rand_Dist_SingleSS{iSS} = [];
%             for nRand = 1:length(gdata_Pairs_ss{iSS})-1
%                 Granger_Rand_Dist_SingleSS{iSS}(nRand,:) = gdata_Pairs_ss{iSS}{nRand+1};
%             end
%         end
        
%         Granger_Rand_Dist_SingleSS = abs(Granger_Rand_Dist_SingleSS{1})-abs(Granger_Rand_Dist_SingleSS{2});
        %% Percentiles
        prcList = [1,5,50,95,99];
        RandPrc = {};
        RandPrc_Im = {};
        for iPrc = 1:length(prcList)
            val  = prctile(abs(Granger_Rand_Dist_SingleSS),prcList(iPrc));
            RandPrc{prcList(iPrc)} = val;
            val  = prctile(imag(Granger_Rand_Dist_SingleSS),prcList(iPrc));
            RandPrc_Im{prcList(iPrc)} = val;
        end
        
        %% Min max
        RandStats.Min = min(abs(Granger_Rand_Dist_SingleSS));
        RandStats.Max = max(abs(Granger_Rand_Dist_SingleSS));
        RandStats.Mean = mean(abs(Granger_Rand_Dist_SingleSS));
        RandStats.Std = std(abs(Granger_Rand_Dist_SingleSS));
        RandStats.NoPermutations = size(Granger_Rand_Dist_SingleSS,1);
        
        RandStats_Im.Min = min(imag(Granger_Rand_Dist_SingleSS));
        RandStats_Im.Max = max(imag(Granger_Rand_Dist_SingleSS));
        RandStats_Im.Mean = mean(imag(Granger_Rand_Dist_SingleSS));
        RandStats_Im.Std = std(imag(Granger_Rand_Dist_SingleSS));
        RandStats_Im.NoPermutations = size(Granger_Rand_Dist_SingleSS,1);
        
        %% Save results
        if Task_Period == 1
            strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Granger_significance\Encoding\']
        elseif Task_Period == 2
            strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Granger_significance\Maintenance\']
        end
        mkdir(strPaths.stats)
        save([strPaths.stats,'Rand_Stats_Patient_',num2str(38,'%.2d'),'Reref_Rand_Granger_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'],...
            'Real_Granger','Real_Granger2','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','cfgFreq','-v7.3')
        clear Real_Granger Real_Granger2 RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis cfgFreq  PLV_Rand_SS
%         clear gdata_Pairs_ss

    end
    tStopRand = cputime-tStartRand;
    
    
    
    %% Load results / Granger for all pairs / Randomized trials /
    %% encoding
    % strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Stats 180405\'];
    strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Granger_significance\Encoding\']
    Granger_Rand_Ret_Variables_Pair_SS = cell(1,size(dataBipolar_Ret_SS,2));
    for iSS = iSS_ToCompute%1:size(dataBipolar_Ret_SS,2)
        Granger_Rand_Ret_Variables_Pair_SS{iSS} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            strVariablePath = [strPaths.stats,'Rand_Stats_Patient_',num2str(38,'%.2d'),'Reref_Rand_Granger_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'];
            try
                Vars = load(strVariablePath);
                Granger_Rand_Ret_Variables_Pair_SS{iSS}{nPair} = Vars;
                clear Vars
            catch
            end
        end
    end
    
    %% maintenance
    % strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Stats 180405\'];
    strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Granger_significance\Maintenance\']
    Granger_Rand_Ret_Variables_Pair_SS_maint = cell(1,size(dataBipolar_Ret_SS,2));
    for iSS = iSS_ToCompute%1:size(dataBipolar_Ret_SS,2)
        Granger_Rand_Ret_Variables_Pair_SS_maint{iSS} = cell(1,size(nChannelPairs,1));
        for nPair = 1:size(nChannelPairs,1)
            nChannel_1 = nChannelPairs(nPair,1);
            nChannel_2 = nChannelPairs(nPair,2);
            strVariablePath = [strPaths.stats,'Rand_Stats_Patient_',num2str(38,'%.2d'),'Reref_Rand_Granger_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'];
            try
                Vars = load(strVariablePath);
                Granger_Rand_Ret_Variables_Pair_SS_maint{iSS}{nPair} = Vars;
                clear Vars
            catch
            end
        end
    end
    
        %%
%     Panel_l_Dat = FigureMaint_Data.Granger;
    Grngr_Enc1 = Granger_Rand_Ret_Variables_Pair_SS{4}{nPairs_ToRun}.Real_Granger;
    Grngr_Enc2 = Granger_Rand_Ret_Variables_Pair_SS{4}{nPairs_ToRun}.Real_Granger2;%Panel_l_Dat.Grng_Encoding{1};
    Grngr_Maint1 =Granger_Rand_Ret_Variables_Pair_SS_maint{4}{nPairs_ToRun}.Real_Granger;
    Grngr_Maint2 =Granger_Rand_Ret_Variables_Pair_SS_maint{4}{nPairs_ToRun}.Real_Granger2;
    freq = [4:1:100];%Panel_l_Dat.freq;
    light_blue = [0.30,0.75,0.93];
    orange       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',orange};
    
    figure;
    hold on;
    semilogx(freq,Grngr_Enc2*100,'Color',Colors{2},'LineWidth',3);
    semilogx(freq,Grngr_Maint1*100,'Color',Colors{3},'LineWidth',3);
    semilogx(freq,Grngr_Maint2*100,'Color',Colors{4},'LineWidth',3);
    semilogx(freq,Grngr_Enc1*100,'Color',Colors{1},'LineWidth',3);
    xlim([4 30]);
    ylim([-0.01 0.2]*100);
    set(gca,'XScale','log','TickDir','out','XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.025]);
    
    Vars = Granger_Rand_Ret_Variables_Pair_SS{iSS}{nPairs_ToRun};
    Vars2 = Granger_Rand_Ret_Variables_Pair_SS_maint{iSS}{nPairs_ToRun};

    indMaxBandCortexHipp_Enc = Difference_Bar((Grngr_Enc1-Grngr_Enc2)*100,Vars.RandPrc{99}*100,freq,[0 0.2], -0.005*100,Colors{1});
    indMaxBandHippCortex_Maint = Difference_Bar((Grngr_Maint1-Grngr_Maint2)*100,Vars2.RandPrc{99}*100,freq,[0 0.2], -0.007*100,Colors{3});

%     indMaxBandHippCortex_Maint = Difference_Bar(Grngr_Maint.grangerspctrm(1,:)-Vars.RandPrc{99},Vars.RandPrc{99},freq,[0 0.2], 0.005,Colors{3});
%%  
    k =1;
    for i = ind(1):ind(end)
        Vars = Granger_Rand_Ret_Variables_Pair_SS{iSS}{i}
        surrogate =  (Granger_Rand_Ret_Variables_Pair_SS{iSS}{i}.Real_Granger) %-(Vars.RandPrc{95}) %StatsVars.RandPrc{95};
        FreqBand = [8,18];
        freqAxis = Vars.freqAxis%[1:100];
        [~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
        [~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
        Granger_temp = surrogate([indFreq1:indFreq2]);
        Granger_In_Band_SS(k) = min(zscore(Granger_temp));
        k = k +1;
    end
    Granger_In_Band_ToPlot = reshape(Granger_In_Band_SS,8,8);
    %%
    figure;
    imagesc(Granger_In_Band_ToPlot,[-5 -1])
    
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
    
    