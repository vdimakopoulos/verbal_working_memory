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



%% Load data for 2 sessions - Macro Data
strMacroDataDir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\12 BP\';
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
 data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end


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
    
    end
end
mkdir(Path_for_PLV_sgnf)


%% Merge sessions
dataBipolar = data_ses(1).data;
for nSes = 2:length(data_ses)
    cfg = [];
    dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
    
end

%% Rereferencing
cfg               =   [];
cfg.reref         =  'yes' 
cfg.refchannel    =  [1 2 3] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 9;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [65 67 67 68 68 69 69 70 70 71 71 72 72 73 73 74 74 66]%[65 67 65 68 67 68 65 69 67 69 68 69] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
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
cfg.refchannel = [75 83]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;


%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);


%% macro data
% Select only correct trials
strChannelNameList = dataBipolar.label;

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


%% Channel Pairs
if Hippocampus_Grid_Coupling
    PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
    clear nChannelPairs;
    
    nGridChans = length(find(contains(dataBipolar.label,'GR')));
    nChannelPairs = [];%[[ones(nStripChans,1),(1:nStripChans)'+AHL_chans];[nStripChans+AHL_chans+ones(nStripChans,1),(1:nStripChans)'+AHL_chans]];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChans,1),(1:nGridChans)']];
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
    save([strVariableFolder,'Patient_',num2str(12,'%.2d '),'_',coupling_type,'_Randomized_Trial_Numbers_200_perms_maint_vs_encoding.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')
    
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
        ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','-')));
        ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','GR')));
        
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
                    cfg.foi         = [1:1:29, 30:5:100];
                    
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
            
            
            if PLV_significance_maint_vs_fix
                strVariableFolder = [strPLVResultFolder,'Grid_Hippocampus\Fix_vs_maint\'];
                mkdir(strVariableFolder)
                save([strVariableFolder,'Patient_',num2str(12),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'_Hippocampus_',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'Fix_vs_maint_mt_tsf_2_Set_Size_freq_1_30_5_100_',num2str(2*iSS+2,'%.1d'),'.mat'],...
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
           
            if Hippocampus_Grid_Coupling
                if PLV_significance_maint_vs_fix
                    strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
                    mkdir(strPaths.stats)
                    save([strPaths.stats,'Patient_',num2str(12),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat'],...
                        'Real_PLV_AcrossPeriods','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
                    clear Real_PLV_AcrossPeriods RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
                end   
                        
            end
        end
        clear PLV_SS;
        tStopRand = cputime-tStartRand;
        %% Load the results
        
        if PLV_significance_maint_vs_fix
            strPaths.stats = [strPaths.Results,'Statistics\Rand_Stats\Maintenance_vs_Encoding\Grid_Hippocampus\Fix_vs_Maint\']
            
        end
        PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,size(Set_Sizes,1));
        for iSS = 1:2%size(Set_Sizes,1)
            PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)} = cell(1,size(nChannelPairs,1));
            for nPair = 1:size(nChannelPairs,1)
                nChannel_1 = nChannelPairs(nPair,1);
                nChannel_2 = nChannelPairs(nPair,2);
                if PLV_significance_maint_vs_fix
                    strVariablePath = [strPaths.stats,'Patient_',num2str(12),'_Pair_',num2str(nPair,'%.4d'),'_Grid_',dataBipolar_Ret_SS{iSS}.label{nChannel_2},'Hippocampus',dataBipolar_Ret_SS{iSS}.label{nChannel_1},'_mt_tsf_2_Stats_Summary_Across_Periods_Maint_vs_Fixation_freq_1_30','.mat']
                end
                try
                    Vars = load(strVariablePath);
                    PLV_Rand_Enc_Maint_Variables_Pair_SS{nSS_Ran(iSS)}{nPair} = Vars;
                    clear Vars
                catch
                end
            end
        end
        
        
        %% Run Once
        strPaths.Figures = [strPaths.Results,'12 BP Results\PSD\Strip Chans\'];
        mkdir(strPaths.Figures)
        strChannelNameList = dataBipolar.label;
        ind2 = find(~cellfun(@isempty,strfind(strChannelNameList()','GR')));
        
        for i =1:length(ind2)
            
            nGridChanNum =i;
            strChanName = strcat('GR',int2str(nGridChanNum));
            gridGR_plot_order(i) = min(find(~cellfun(@isempty,strfind(strChannelNameList()', strChanName))));
            
        end
        nGridChans = length(find(contains(dataBipolar.label,'GR')));
        k = 0;
        strGridRows = {'A','B','C','D','E','F','G','H'};
        for i = 1:nGridChans
            
            nGridNum{i} = str2double(strsplit(dataBipolar.label{i},'GR'))
            
            nChanNum(i) = nGridNum{i}(2);
            nRowLetter(i) = ceil(nChanNum(i)/8);
            dataBipolar.label{i}  = strcat(strGridRows{nRowLetter(i)},int2str((nChanNum(i)-(nRowLetter(i)-1)*8)));
            
        end

        %% Plot the PLV results for every bipolar channel

        iSS_to_Plot = 4;
        k =0;
        % Load the PLV results
        strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\'
        PLV_Vars = load([strVariableFolder,'Patient_',num2str(12,'%.2d '),'_','Grid_Depth','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'])
        PLV_Pairs_maint = PLV_Vars.PLV_Pairs;
        freqAxis = PLV_Vars.freqAxis;
        
        PLV_Vars = load([strVariableFolder,'Patient_',num2str(12,'%.2d '),'_','Grid_Depth','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'])
        PLV_Pairs_fix = PLV_Vars.PLV_Pairs;
        nChannelPairs = PLV_Vars.nChannelPairs;
        freqAxis_fix = PLV_Vars.freqAxis;
        
        for i=1:nPairs_ToRun(end)
            
            if i==1 || ~mod(i-1,nGridChans)
                fig = figure;
                
                ha = tight_subplot(8,8,[.03 .03],[.1 .1],[.1 .1])
                strTitle =['PLV of ',dataBipolar.label(nChannelPairs(i,1)),' and grid GR, ss ',num2str(Set_Sizes{iSS_to_Plot})]
                suptitle(strTitle);
                if i~=1
                    k = floor(i/nGridChans);
                end
            end
            ind  = find(gridGR_plot_order == (i-(k*64)))
            axes(ha(ind));
            semilogx(freqAxis_fix,abs(PLV_Pairs_fix{i,iSS_to_Plot}),'k','LineWidth',3);
            hold on;
            semilogx(freqAxis,abs(PLV_Pairs_maint{i,iSS_to_Plot}),'r','LineWidth',3);
            xlim([4 100]);
            ylim([0 1]);
            Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS_to_Plot}{i};
            Fix_Maint_Diff = abs(PLV_Pairs_fix{i,iSS_to_Plot})-abs(PLV_Pairs_maint{i,iSS_to_Plot});
            indMaxBand = Difference_Bar(Fix_Maint_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'k');
            Maint_Fix_Diff = abs(PLV_Pairs_maint{i,iSS_to_Plot})-abs(PLV_Pairs_fix{i,iSS_to_Plot});
            indMaxBand = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis,[0 1],0.03,'r');
            ylabel(dataBipolar.label((i-(k*64))));
            
        end
end