%% Close all figures, clear variables and command window
close all
clear
clc

%% Paths
% Folder that contains the folder 'Sternberg Task'
strPaths.Main = 'I:\MATLAB Codes\';
strPaths.GeneralFunctions = [strPaths.Main,'General Functions\'];

% Main folder of task and subfolders
strPaths.Sternberg                                  = [strPaths.Main,'Sternberg Task\'];
strPaths.EEGLABDataClean.Main                             = [strPaths.Sternberg,'EEGLAB Data Clean 1\'];
strPaths.EEGLABDataClean.Patients.Main                    = [strPaths.EEGLABDataClean.Main,'Patients\'];
strPaths.EEGLABDataClean.Patients.Macro                   = [strPaths.EEGLABDataClean.Patients.Main,'Macro Data\'];
strPaths.EEGLABDataClean.Patients.Macro_Bipolar           = [strPaths.EEGLABDataClean.Patients.Main,'Macro Bipolar Data\'];
strPaths.EEGLABDataClean.Patients.Scalp                   = [strPaths.EEGLABDataClean.Patients.Main,'Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Macro_Bipolar_Scalp     = [strPaths.EEGLABDataClean.Patients.Main,'Macro Bipolar Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Micro                   = [strPaths.EEGLABDataClean.Patients.Main,'Micro Data\'];
strPaths.EEGLABDataClean.Patients.First_3_Bipolar_Scalp   = [strPaths.EEGLABDataClean.Patients.Main,'Macro First 3 Bipolar Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Params                  = [strPaths.EEGLABDataClean.Patients.Main,'Artifact Rejection Params\'];

% Recording information
strPaths.PatientRecordingInformation = [strPaths.Sternberg,'Patient Information\'];

% Results
strPaths.Results = [strPaths.Sternberg,'Results\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'I:\MATLAB Codes\Toolboxes\fieldtrip-20170925\';
% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = 'I:\MATLAB Codes\Toolboxes\eeglab14_1_1b\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(genpath(strPaths.Sternberg))
addpath(genpath(strPaths.GeneralFunctions))

% Add EEGLAB to path
addpath(strPaths.Toolboxes.EEGLAB)
eeglab
close all

% Remove FieldTrip toolbox from path
rmpath(strPaths.Toolboxes.FieldTrip)

% Plot colors
strPlotColors = {'b','g','r','k'};

%% Check the contents of the folder for patient data with macro bipolar and scalp channels
% Get list of patients with data in this folder
strFieldName_Main = 'First_3_Bipolar_Scalp';
strPatientDataPath_Main = strPaths.EEGLABDataClean.Patients.(strFieldName_Main);
[nPatientList,nPatientNumberFile,strPatientFileList] = Get_Patient_List_In_Folder(strPatientDataPath_Main);
nNumberOfPatients = length(nPatientList);

%% Select patient
nPatient = 28;
iPatient = find(nPatientList==nPatient);
strPatientFileName_Main = strPatientFileList{iPatient};

%% Run certain parts of code
flag_Run_Load.Load_Data = 1;
flag_Run_Load.Run_Real_PLV_Ret = 1;
flag_Run_Load.Load_Real_PLV_Ret = 0;
flag_Run_Load.Plot_Real_PLV_Ret = 0;
flag_Run_Load.Rand_PLV_TrialNumberList = 0;

flag_Run_Load.Run_Rand_PLV_Ret = 0;

flag_Run_Load.Run_Real_PLV_Stim = 0;
flag_Run_Load.Load_Real_PLV_Stim = 0;

flag_Run_Load.Run_Real_PLV_Ret_Depth = 1;

flag_Run_Load.Run_Real_PLV_CrossSpectrum = 1;

%% Load patient data
if flag_Run_Load.Load_Data % Load already saved variables if 0
    %% Load clean macro bipolar data
    % Macro bipolar scalp
    strFieldName = 'First_3_Bipolar_Scalp';
    strPatientDataPath = strrep(strPatientDataPath_Main,strrep(strFieldName_Main,'_',' '),strrep(strFieldName,'_',' '));
    strPatientFileName = strrep(strPatientFileName_Main,strFieldName_Main,strFieldName);
    EEG_Data = load([strPatientDataPath,strPatientFileName]);
    
    EEGBipolarScalp         = EEG_Data.EEGBipolarScalp;
    TrialInformationTable   = EEG_Data.TrialInformationTable;
    EEG_Data = rmfield(EEG_Data,'EEGBipolarScalp');
    
    %% Clean artifacts if there are any remaining
    nTrialsRejected_3 = Get_PLV_FT_TrialsRejected_3_Patient_180105(nPatient);    
    EEGBipolarScalp = pop_rejepoch(EEGBipolarScalp,nTrialsRejected_3,0);
    
    %% Update trial information table and artifact rejection parameters
    TrialInformationTable(nTrialsRejected_3,:) = [];
    TrialInformationTable.TrialNumber(:) = 1:size(TrialInformationTable,1);
    EEG_Data.ArtifactRejectionParams.nTrialsRejected_Steps{3} = nTrialsRejected_3;
    
    %% Convert into FieldTrip data
    dataBipolarScalp = eeglab2fieldtrip(EEGBipolarScalp,'preprocessing');
    dataBipolarScalp = Get_FieldTrip_Events_From_EEGLAB_Events(dataBipolarScalp,EEGBipolarScalp);
    clear EEG_Data EEGBipolarScalp ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG eeglabUpdater LASTCOM PLUGINLIST STUDY
    
    %% Remove EEGLAB from path and add FieldTrip
    rmpath(genpath(strPaths.Toolboxes.EEGLAB))
    addpath(strPaths.Toolboxes.FieldTrip)
    ft_defaults
    
    %% Plot data
    cfg = [];
    cfg.viewmode = 'vertical';
    cfg.continuous = 'no';
    % ft_databrowser(cfg,dataBipolarScalp);
    
    %% Resample at 200 Hz
    fs = 200;
    cfg = [];
    cfg.resamplefs = fs;
    cfg.detrend = 'no';
    dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);
    
    %% Add channels
    strChannelNamesAdded = Get_PLV_FT_Channel_Names_Added_Patient_180105(nPatient);
    dataBipolarScalp = Get_Added_Channels_Dataset_From_Channel_Names_FieldTrip(dataBipolarScalp,strChannelNamesAdded);
    
    %% Select only correct trials
    [dataBipolarScalp,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable);
    
    %% Divide data into set sizes
    % 4 6 8 6&8
    [dataBipolarScalp_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolarScalp,TrialInformationTable);
    
    %% Scalp and depth channels for phase coherence analysis
    strChannelNameList = dataBipolarScalp.label;
    % Scalp EEG channels
    nChannelList_Scalp = find(cellfun(@isempty,strfind(strChannelNameList,'-')));
    % Bipolar EEG channels / channel 1-2
    nChannelList_Depth = find(~cellfun(@isempty,strfind(strChannelNameList,'-')));
    % nChannelList_Depth = find(~cellfun(@isempty,strfind(strChannelNameList,'1-2')));
    % nChannelList_Depth = [nChannelList_Depth,find(~cellfun(@isempty,strfind(strChannelNameList,'2-4')))];
    
    %% Pairs of channels for phase coherence analysis
    nChannelPairs = [];
    for iChannel1 = 1:length(nChannelList_Scalp)
        for iChannel2 = 1:length(nChannelList_Depth)
            nChannelPairs = [nChannelPairs;[nChannelList_Scalp(iChannel1),nChannelList_Depth(iChannel2)]];
        end
    end
    
    %% Number of trials for each set size
    nNumberOfTrialsSetSizes = zeros(1,4);
    for iSS = 1:4
        nNumberOfTrialsSetSizes(iSS) = length(dataBipolarScalp_SS{iSS}.trial);
    end
    
    %% Extract 2 seconds of retention
    for iSS = 1:4
        cfg = [];
        cfg.latency = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        dataBipolarScalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolarScalp_SS{iSS});
    end
    
    %% Save variables
    strVariablePath = [strPaths.Results,'PLV FieldTrip\Load Patient Data Variables 180405\'];
    mkdir(strVariablePath)
    save([strVariablePath,sprintf('Patient_%.2d_Load_Patient_Data_Variables_180405.mat',nPatient)],'-regexp','^(?!(strPaths|strPatientDataPath_Main|strPatientDataPath|strVariablePath|flag_Run_Load|ans)$).')
    
else
    %% Remove EEGLAB from path and add FieldTrip
    rmpath(genpath(strPaths.Toolboxes.EEGLAB))
    addpath(strPaths.Toolboxes.FieldTrip)
    ft_defaults
    
    %% Load variables
    strVariablePath = [strPaths.Results,'PLV FieldTrip\Load Patient Data Variables 180405\'];
    load([strVariablePath,sprintf('Patient_%.2d_Load_Patient_Data_Variables_180405.mat',nPatient)])
end

%% Run PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Run_Real_PLV_Ret
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs
    % ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','1-2')));
    % nPairs_ToRun = intersect(nPairs_ToRun,ind);
    
    %% Analysis for each pair of electrodes - Real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Ret Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.foi         = 0.5:0.5:30;
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1:2;
            Freq            = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            %% Coherence
            cfg             = [];
            cfg.method      = 'coh';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
        end
        clear Freq cfg
        
        %% Store results for the channel pair
        strPairLabels = PLV_SS{iSS}.label';
        freqAxis = PLV_SS{iSS}.freq;
        for iSS = 1:4
            PLV_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
            Coh_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
        end
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'],'PLV_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','Coh_SS','-v7.3')
        clear PLV_SS Coh_SS strPairLabels freqAxis cfgFreq
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
    end
end

%% Load results / PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Load_Real_PLV_Ret % Run results for patient from folder if 0    
    %% Load results for all pairs
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Ret Tapsmofrq 2 Values 180405\'];
    PLV_Real_Ret_SS_Pair = cell(1,size(nChannelPairs,1));
    Coh_Real_Ret_SS_Pair = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Real_Ret_SS_Pair{nPair} = Vars.PLV_SS;
            Coh_Real_Ret_SS_Pair{nPair} = Vars.Coh_SS;
            freqAxis = Vars.freqAxis;
            clear Vars nChannel_1 nChannel_2 strVariablePath
        catch
            clear Vars nChannel_1 nChannel_2 strVariablePath
        end
    end
    strPairLabels_Pair = dataBipolarScalp_Ret_SS{iSS}.label(nChannelPairs);
    
    %% Save real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Ret Tapsmofrq 2 Values 180405\'];
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_mt_tsf_2.mat'],'PLV_Real_Ret_SS_Pair','strPairLabels_Pair','freqAxis','Coh_Real_Ret_SS_Pair','-v7.3')
    
else
    %% Load real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Ret Tapsmofrq 2 Values 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_mt_tsf_2.mat'])
    
end

%% Plot the real values on for each depth electrode - real values
if flag_Run_Load.Plot_Real_PLV_Ret
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Abs Ret Tapsmofrq 2 Scalp Layouts 180405\'];
mkdir(strImageFolder)
figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = 1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs(:,2)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs(nPairs_WithDepthChannel,1);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs(:,1)==nChannelScalp&nChannelPairs(:,2)==nChannel_Depth);
        
        %% Subplot number for scalp channel
        nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
        nSubplot = figLayout{nEntryInTable,2};
        subplot(5,5,nSubplot)
        
        %% Plot
        for iSS = 1:3 %[1,4,2,3]
            plot(freqAxis,abs(PLV_Real_Ret_SS_Pair{nPairs_WithDepthScalpChannel}{iSS}),strPlotColors{iSS},'LineWidth',2)
            hold on
            xlim([1,30])
            ylim([0,0.8])
            set(gca,'XTick',[1,5:5:30])
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
        end
        %%
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.fig'],'fig')
    close all
    
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList
end
end

%% Random trial numbers for each set size separately
if flag_Run_Load.Rand_PLV_TrialNumberList
    nNumberOfPermutations = 1000;
    nRandomTrialNumberList_Pair_SS = cell(1,4);
    for iSS = 1:4
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
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Trial Numbers 180405\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_Pair_SS','nNumberOfPermutations','-v7.3')
else
    %% Load variables
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Trial Numbers 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'])
end

%% Run PLV for all pairs / Randomized trials / 2 seconds of retention
if flag_Run_Load.Run_Rand_PLV_Ret
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Values 180405\'];
    strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Stats 180405\'];
    mkdir(strVariableFolder)
    mkdir(strStatsFolder)

    %% Select pairs for analysis
    thrPLV = 0.2;
    iSS = 3;
    flagRunAnalysis_Pair = zeros(size(nChannelPairs,1),1);
    for iPair = 1:size(nChannelPairs,1)
        flagRunAnalysis_Pair(iPair) = ~isempty(find(abs(PLV_Real_Ret_SS_Pair{iPair}{iSS})>thrPLV,1));
    end
    nPairs_ToRun = find(flagRunAnalysis_Pair);
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    
    %% Analysis for each pair of electrodes - Randomization
    iSS_ToCompute = 3;
    tStartRand = cputime;
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = iSS_ToCompute
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
        end
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = cell(1,4);
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
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Convert cell to double
        for iSS = 1:4
            try
                PLV_Rand_SS{iSS} = cell2mat(PLV_Rand_SS{iSS}');
            catch
                PLV_Rand_SS{iSS} = [];
            end
        end
        
        %% Save values
        iSS = iSS_ToCompute;
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_ss',num2str(2*iSS+2),'.mat'],...
            'PLV_Rand_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
                
        %% Real values
        nRand = 1;
        Real_PLV = PLV_Rand_SS{iSS}(1,:);
        
        %% Null distribution
        PLV_Rand_Dist_SingleSS = PLV_Rand_SS{iSS}(2:end,:);
        
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
        save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'],...
            'Real_PLV','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
        clear Real_PLV RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis cfgFreq TimeInterval PLV_Rand_SS
        
    end
    
    tStopRand = cputime-tStartRand;
    
end

%% Load results / PLV for all pairs / Randomized trials / 2 seconds of retention
strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Ret Tapsmofrq 2 Stats 180405\'];
PLV_Rand_Ret_Variables_Pair_SS = cell(1,4);
for iSS = 1:4
    PLV_Rand_Ret_Variables_Pair_SS{iSS} = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Rand_Ret_Variables_Pair_SS{iSS}{nPair} = Vars;
            clear Vars
        catch
        end
    end
end

%% Plot results for each depth electrode for randomized values
iSS = 3;
strImageFolder = [strPaths.Sternberg,'Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Scalp Layout Plots 180320\'];
figLayout = Get_Figure_Scalp_Layout_Information();
% Initialize plot
Max_PLV_Pairs_Table = zeros(size(nChannelPairs,1),6);
Max_PLV_Pairs_Table = array2table(Max_PLV_Pairs_Table);
Max_PLV_Pairs_Table.Properties.VariableNames = {'Scalp','Depth','Max_PLV','freq_Max_PLV','AreaUnder_PLV','freqBand_AreaUnder_PLV'};
Max_PLV_Pairs_Table.Scalp = cell(size(nChannelPairs,1),1);
Max_PLV_Pairs_Table.Depth = cell(size(nChannelPairs,1),1);
Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV = cell(size(nChannelPairs,1),1);
% Plot layout for each depth electrode
for iChannel_Depth = 1:3:length(nChannelList_Depth) % 13
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    % Depth channel
    nChannel_Depth = nChannelList_Depth(iChannel_Depth);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    % Scalp channels
    indDepthChannelPairs = find(nChannelPairs(:,2)==nChannel_Depth);
    nChannelList_Scalp = nChannelPairs(indDepthChannelPairs,1);
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp);
    % Plot results on layout
    for iScalp = 1:length(strScalpChannelNameList)
        %% Scalp channel number and label
        nChannel_Scalp = nChannelList_Scalp(iScalp);
        strScalpChannelName = strScalpChannelNameList{nChannel_Scalp};
        %% Pair
        nPair = find(nChannelPairs(:,1)==nChannel_Scalp&nChannelPairs(:,2)==nChannel_Depth);
        %% Subplot in layout
        indSubplot = find(strcmpi(figLayout(:,1),strScalpChannelName));
        nSubplot = figLayout{indSubplot,2};
        subplot(5,5,nSubplot)
        %% Plot results
        Vars = PLV_Rand_Ret_Variables_Pair_SS{iSS}{nPair};
        if(~isempty(Vars))
            % Frequency axis
            freqAxis = Vars.freqAxis;
            % Real PLV
            Real_PLV = Vars.Real_PLV;
            
            %% Percentiles
            prcList = [1,5,50,95,99];
            for iPrc = 1:length(prcList)
                Vars.RandPrc{prcList(iPrc)};
                Vars.RandPrc{prcList(iPrc)};
            end
            
            %% Plot
            plot(freqAxis,abs(Real_PLV),strPlotColors{iSS})
            hold on
            shadedErrorBar(freqAxis,Vars.RandPrc{50},[Vars.RandPrc{95}-Vars.RandPrc{50};Vars.RandPrc{50}-Vars.RandPrc{5}],{'Color',0.5*[1,1,1]});
            plot(freqAxis,Vars.RandStats.Mean)
            plot(freqAxis,Vars.RandStats.Mean+Vars.RandStats.Mean,'m')
            plot(freqAxis,Vars.RandStats.Mean-Vars.RandStats.Mean,'m')
            plot(freqAxis,Vars.RandPrc{99},'g')
            plot(freqAxis,abs(Real_PLV),strPlotColors{iSS})
            
            if(~(strcmpi(Vars.strPairLabels{1},strScalpChannelName)&&strcmpi(Vars.strPairLabels{2},strDepthChannelName)))
                error('Channel labels do not match')
            end
            
            %% Axis limits
            ylim([0,1])
            xlim([1,30])
            set(gca,'XTick',[1,5:5:30])
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
            
            %% Chanel names
            Max_PLV_Pairs_Table.Scalp{nPair} = strScalpChannelName;
            Max_PLV_Pairs_Table.Depth{nPair} = strDepthChannelName;
            
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
                    stem(freqAxis(indMaxBand),Real_PLV(indMaxBand))
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
                %                 pl = area(freqAxis(indMaxBand),x(indMaxBand));
                title(sprintf('%s - %s - %.2f',strScalpChannelName,strDepthChannelName,AreaUnder_PLV_Prc))
                flagSig_SS6(nPair) = length(cc.PixelIdxList);
            end
            
            %%
        end
    end
    
    %%
    if 0
        %% Select band
        Max_PLV_Pairs_Table_SingleDepthChannel = Max_PLV_Pairs_Table(strcmpi(Max_PLV_Pairs_Table.Depth,strDepthChannelName),:);%&(~(strcmpi(Max_PLV_Pairs_Table.Scalp,'A1')|strcmpi(Max_PLV_Pairs_Table.Scalp,'A2'))),:);
        Max_PLV_Pairs_Table_SingleDepthChannel((strcmpi(Max_PLV_Pairs_Table_SingleDepthChannel.Scalp,'A1')|strcmpi(Max_PLV_Pairs_Table_SingleDepthChannel.Scalp,'A2')),:).AreaUnder_PLV = zeros(2,1);
        [~,indMaxArea] = max(Max_PLV_Pairs_Table_SingleDepthChannel.AreaUnder_PLV);
        freqBandSelected = Max_PLV_Pairs_Table_SingleDepthChannel.freqBand_AreaUnder_PLV{indMaxArea};
        
        
        %% Area under curve for the selected band
        for iScalp = 1:length(strScalpChannelNameList)
            %% Scalp channel number and label
            nChannel_Scalp = nChannelList_Scalp(iScalp);
            strScalpChannelName = strScalpChannelNameList{nChannel_Scalp};
            %% Pair
            nPair = find(nChannelPairs(:,1)==nChannel_Scalp&nChannelPairs(:,2)==nChannel_Depth);
            %% Plot results
            Vars = PLV_Rand_Variables_Pairs{iSS}{nPair};
            if(~isempty(Vars))
                % Frequency axis
                freqAxis = Vars.freqAxis;
                % Real PLV
                Real_PLV = Vars.Real_PLV;
                Real_PLV = PLV_Real_SS_Pair{nPair}{1};
                
                %% Area under curve for selected band
                thrPrc = 99;
                thrNeighborPts = 3;
                thrFreq = 2;
                Real_PLV = abs(Real_PLV);
                Thr_PLV = Vars.RandPrc{thrPrc};
                indFreqBandSelected = find(ismember(freqAxis,freqBandSelected));
                Real_PLV(indFreqBandSelected)-Thr_PLV(indFreqBandSelected)
                
                indAboveThr = Real_PLV>=Thr_PLV;
                cc = bwconncomp(indAboveThr);
                if(cc.NumObjects~=0)
                    % Threshold number of points
                    numPixels = cellfun(@numel,cc.PixelIdxList);
                    cc.PixelIdxList = cc.PixelIdxList(numPixels>=thrNeighborPts);
                    % Threshold lower bound for frequency
                    cc.PixelIdxList = cc.PixelIdxList(freqAxis(cellfun(@min,cc.PixelIdxList))>=thrFreq);
                    % Largest cluster
                    numPixels = cellfun(@numel,cc.PixelIdxList);
                    [~,indMaxCluster] = max(numPixels);
                    % Cluster with the highest area under
                    plv_all = cellfun(@(x) Real_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    thr_plv_all = cellfun(@(x) Thr_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    
                    if(~isempty(indMaxCluster))
                        indMaxBand = cc.PixelIdxList{indMaxCluster};
                        freqBand = freqAxis(indMaxBand);
                        AreaUnder_PLV_Prc = sum(Real_PLV(indMaxBand)-Thr_PLV(indMaxBand));
                        %                     stem(freqAxis(indMaxBand),Real_PLV(indMaxBand))
                    else
                        freqBand = [];
                        AreaUnder_PLV_Prc = 0;
                    end
                    Max_PLV_Pairs_Table.AreaUnder_PLV(nPair) = AreaUnder_PLV_Prc;
                    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV{nPair} = freqBand;
                    %                 pl = area(freqAxis(indMaxBand),x(indMaxBand));
                end
                
                %%
            end
        end
    end
    %% Save image
%     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Depth_Channel_',num2str(iChannel_Depth,'%.3d'),'_',strDepthChannelName,'_mt_tsf_2_Rand_Set_Size_',num2str(2*iSS+2,'%.1d'),'.png'],'png')
%     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Depth_Channel_',num2str(iChannel_Depth,'%.3d'),'_',strDepthChannelName,'_mt_tsf_2_Rand_Set_Size_',num2str(2*iSS+2,'%.1d'),'.fig'],'fig')
    close all
end

%%
Max_PLV_Pairs_Table_2 = Max_PLV_Pairs_Table(Max_PLV_Pairs_Table.Max_PLV~=0,:);
Max_PLV_Pairs_Table_2 = Max_PLV_Pairs_Table_2(Max_PLV_Pairs_Table_2.AreaUnder_PLV~=0,:);
[~,indMax] = max(Max_PLV_Pairs_Table_2.AreaUnder_PLV)
Max_PLV_Pairs_Table_2(indMax,:)

[~,indSort] = sort(Max_PLV_Pairs_Table_2.AreaUnder_PLV,'descend');
Max_PLV_Pairs_Table_2(indSort,:)

%%
addpath(strPaths.Toolboxes.EEGLAB)
eeglab

%%
% load(['I:\MATLAB Codes\Sternberg Task\EEGLAB Data Clean 1\Patients\Macro First 3 Bipolar Scalp Data\First_3_Bipolar_Scalp_EEGLAB_Data_Patient_28.mat'],'EEGBipolarScalp')
load(['I:\MATLAB Codes\Sternberg Task\EEGLAB Data Raw\Patients\Scalp Data\','Scalp_EEGLAB_Data_Patient_28.mat'],'EEGScalp')
load(['I:\MATLAB Codes\Sternberg Task\EEGLAB Data Clean 1\Patients\Artifact Rejection Params\','Artifact_Rejection_Params_EEGLAB_Data_Patient_28.mat'],'ArtifactRejectionParams')
EEGBipolarScalp

%%
Max_PLV_Pairs_Table_SingleDepth = Max_PLV_Pairs_Table_2(ismember(Max_PLV_Pairs_Table_2.Depth,'PHL1-2'),:);
Max_PLV_Pairs_Table_SingleDepth = Max_PLV_Pairs_Table_SingleDepth(~ismember(Max_PLV_Pairs_Table_SingleDepth.Scalp,{'A1','A2'}),:);

% Max_PLV_Pairs_Table_SingleDepth_Temp =
sortrows(Max_PLV_Pairs_Table_SingleDepth,'AreaUnder_PLV_Prc_SelectedBand','descend')


datavectorToPlot = Max_PLV_Pairs_Table_SingleDepth.AreaUnder_PLV_Prc_SelectedBand;
datavectorToPlot(datavectorToPlot<0) = 0;

ChanlocsToPlot = EEGBipolarScalp.chanlocs(1:17);

figure(101)
topoplot(datavectorToPlot,ChanlocsToPlot,'maplimits',[0,5])
colorbar
title(sprintf('Mean PLV in %.1f-%.1f Hz during retention',FreqBand(1),FreqBand(2)))

%%
if 0
    
    %% Max values for each pair
    iSS = 3;
    PLV_Max_Pairs = zeros(size(nChannelPairs,1),1);
    for nPair = 1:length(PLV_Real_Ret_SS_Pair)
        try
            vals = PLV_Real_Ret_SS_Pair{nPair}{iSS};
            PLV_Max_Pairs(nPair) = max(vals);
        catch
        end
    end
    
    figure,plot(sort(PLV_Max_Pairs))
    
    %% Table for max PLV values
    Max_PLV_Pairs_Table = zeros(length(PLV_Real_Ret_SS_Pair),3);
    Max_PLV_Pairs_Table = array2table(Max_PLV_Pairs_Table);
    Max_PLV_Pairs_Table.Properties.VariableNames = {'Scalp','Depth','Max_PLV'};
    Max_PLV_Pairs_Table.Scalp = cell(length(PLV_Real_Ret_SS_Pair),1);
    Max_PLV_Pairs_Table.Depth = cell(length(PLV_Real_Ret_SS_Pair),1);
    for nPair = 1:length(PLV_Real_Ret_SS_Pair)
        Max_PLV_Pairs_Table.Scalp(nPair) = dataBipolarScalp_Ret_SS{3}.label(nChannelPairs(nPair,1));
        Max_PLV_Pairs_Table.Depth(nPair) = dataBipolarScalp_Ret_SS{3}.label(nChannelPairs(nPair,2));
        Max_PLV_Pairs_Table.Max_PLV(nPair) = PLV_Max_Pairs(nPair);
    end
    
    %%
    Max_PLV_Pairs_Table_Sorted = sortrows(Max_PLV_Pairs_Table,'Max_PLV','descend');
    Max_PLV_Pairs_Table_Sorted = Max_PLV_Pairs_Table_Sorted(cellfun(@isempty,strfind(Max_PLV_Pairs_Table_Sorted.Scalp,'A')),:);
    cond1 = ~cellfun(@isempty,strfind(Max_PLV_Pairs_Table_Sorted.Depth,'1-2'));
    % cond2 = ~cellfun(@isempty,strfind(Max_PLV_Pairs_Table_Sorted.Depth,'PHL2-4'));
    cond2 = ismember(Max_PLV_Pairs_Table_Sorted.Depth,{'PHL2-4','PHL2-3'});
    % cond3 = ~cellfun(@isempty,strfind(Max_PLV_Pairs_Table_Sorted.Scalp,'A1'));
    % cond4 = ~cellfun(@isempty,strfind(Max_PLV_Pairs_Table_Sorted.Scalp,'A2'));
    cond = cond1|cond2;
    Max_PLV_Pairs_Table_Sorted = Max_PLV_Pairs_Table_Sorted(cond,:);
    Max_PLV_Pairs_Table_Sorted(1:20,:)
    
    % figure,plot(Max_PLV_Pairs_Table_Sorted.Max_PLV)
    
    %%
    strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 180315 Variables\';
    % PLV_Dist_SinglePair = load([strVariableFolder,'Patient_28_Pair_0265_Scalp_O2_Depth_Channel_PHL1-2_mt_tsf_2']);
    % PLV_Dist_SinglePair = PLV_Dist_SinglePair.PLV_SS;
    iSS = 3;
    nRand = 1;
    RealPLV = abs(squeeze(PLV_Dist_SinglePair{iSS}{nRand}.plvspctrm(1,2,:)))';
    RandDist = [];
    for nRand = 1:length(PLV_Dist_SinglePair{iSS})-1
        RandDist(nRand,:) = abs(squeeze(PLV_Dist_SinglePair{iSS}{nRand+1}.plvspctrm(1,2,:)))';
    end
    freqAxis = 2:0.5:30;
    
    prcList = [1,5,50,95,99];
    RandPrc = {};
    for iPrc = 1:length(prcList)
        val  = prctile(RandDist,prcList(iPrc));
        RandPrc{prcList(iPrc)} = val;
    end
    
    figure
    plot(freqAxis,RealPLV,'r','LineWidth',1)
    hold on
    plot(freqAxis,RandPrc{95},'k')
    
    x = RealPLV;
    y = RandPrc{95};
    ind = x>y;
    cc = bwconncomp(ind);
    
    numPixels = cellfun(@numel,cc.PixelIdxList);
    cc.PixelIdxList = cc.PixelIdxList(numPixels>=3);
    numPixels = cellfun(@numel,cc.PixelIdxList);
    [biggest,idx] = max(numPixels);
    
    indMaxBand = cc.PixelIdxList{idx};
    freqBandWidth = freqAxis(indMaxBand(end))-freqAxis(indMaxBand(1));
    
    sum(x(indMaxBand)-y(indMaxBand))
    
    plot(freqAxis(indMaxBand),x(indMaxBand),'r','LineWidth',2)
    plot(freqAxis(indMaxBand),y(indMaxBand),'k','LineWidth',2)
    
end

%%
%% Set size comparison
nSS_Ran = [3,1];
nNumberOfPermutations = 1000;
nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
nTrialNumbersAll_AcrossSS = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
nRandomTrialNumberList_AcrossSS_Pairs = cell(1,size(nChannelPairs,1));
for iPair = 1:size(nChannelPairs,1) % length(nPairs_ToRun)
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
    nRandomTrialNumberList_AcrossSS_Pairs{nPair} = nRandomTrialNumberList;
end

strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Across SS 8 vs 4 Trial Numbers 180320\';
mkdir(strVariableFolder)
save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_AcrossSS_Pairs','-v7.3')

%% Channel pairs to run analysis for
thrPLV = 0.2;
iSS = 3;
flagRunAnalysis_Pair = zeros(size(nChannelPairs,1),1);
for iPair = 1:size(nChannelPairs,1)
    flagRunAnalysis_Pair(iPair) = length(find(abs(PLV_Real_SS_Pair{iPair}{iSS})>thrPLV))>3;
end
nPairs_ToRun = find(flagRunAnalysis_Pair);
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','1-2')));
nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,find(flagSig_SS6));

%% Analysis for each pair of electrodes - Randomization
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Values Across SS 8 vs 4 180320\';
tStartRand = cputime;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,2);
    for iSS_AcrossSS = 1:length(nSS_Ran)
        iSS = nSS_Ran(iSS_AcrossSS);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
    end
    
    %% Randomly mix trials between set sizes
    for nRand = 1:nNumberOfPermutations+1
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        iTrialList = nRandomTrialNumberList_AcrossSS_Pairs{nPair}(nRand,:);
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
        for iSS = 1:2
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.tapsmofrq   = 2;
            cfg.pad         = 2;
            cfg.foi         = 0.5:0.5:30;
            cfg.channel     = 1:2;
            Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            %% Coherence
            % cfg             = [];
            % cfg.method      = 'coh';
            % cfg.complex     = 'complex';
            % Coh_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS}{nRand});
            
        end
    end
    clear Freq_SS
    
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}{1}.label';
    freqAxis = PLV_SS{iSS}{1}.freq;
    for iSS = 1:2
        for nRand = 1:5%nNumberOfPermutations+1
            PLV_AcrossSS{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
        end
    end
    clear PLV_SS
    
    %% Save results
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
        'PLV_AcrossSS','strPairLabels','freqAxis','-v7.3')
    
    %% Save statistical summary
    %% Real values
    nRand = 1;
    for iSS = 1:2
        Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
    end
    
    %% Null distribution
    for iSS = 1:2
        RandDist_PLV_All_AcrossSS{iSS} = [];
        for nRand = 1:length(PLV_AcrossSS{iSS})-1
            RandDist_PLV_All_AcrossSS{iSS}(nRand,:) = PLV_AcrossSS{iSS}{nRand+1};
        end
    end
    clear PLV_AcrossSS
    
    PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
    RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
    
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
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Stats Summary Across SS 8 vs 4 180320\';
    mkdir(strStatsFolder)
    save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'],...
        'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
    clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end

tStopRand = cputime-tStartRand;

%%
if 0
    %% Analysis for each pair of electrodes - Randomization
    strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Values Across SS 8 vs 4 180320\';
    tStartRand = cputime;
    for iPair = 1:length(nPairs_ToRun) % TEMP
        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Save results
        load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
            'PLV_AcrossSS','strPairLabels','freqAxis')
        
        %% Save statistical summary
        %% Real values
        nRand = 1;
        for iSS = 1:2
            Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
        end
        
        %% Null distribution
        for iSS = 1:2
            RandDist_PLV_All_AcrossSS{iSS} = [];
            for nRand = 1:length(PLV_AcrossSS{iSS})-1
                RandDist_PLV_All_AcrossSS{iSS}(nRand,:) = PLV_AcrossSS{iSS}{nRand+1};
            end
        end
        clear PLV_AcrossSS
        
        RandDist_PLV_All = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
        RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
        
        %% Percentiles
        prcList = [1,5,50,95,99];
        RandPrc = {};
        RandPrc_Im = {};
        for iPrc = 1:length(prcList)
            val  = prctile(RandDist_PLV_All,prcList(iPrc));
            RandPrc{prcList(iPrc)} = val;
            val  = prctile(RandDist_PLV_All_Im,prcList(iPrc));
            RandPrc_Im{prcList(iPrc)} = val;
        end
        
        %% Min max
        RandStats.Min = min(RandDist_PLV_All);
        RandStats.Max = max(RandDist_PLV_All);
        RandStats.Mean = mean(RandDist_PLV_All);
        RandStats.Std = std(RandDist_PLV_All);
        RandStats.NoPermutations = size(RandDist_PLV_All,1);
        
        RandStats_Im.Min = min(RandDist_PLV_All_Im);
        RandStats_Im.Max = max(RandDist_PLV_All_Im);
        RandStats_Im.Mean = mean(RandDist_PLV_All_Im);
        RandStats_Im.Std = std(RandDist_PLV_All_Im);
        RandStats_Im.NoPermutations = size(RandDist_PLV_All_Im,1);
        
        %% Save results
        strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Stats Summary Across SS 8 vs 4 180320\';
        mkdir(strStatsFolder)
        save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'],...
            'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
        clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
        
    end
    
    tStopRand = cputime-tStartRand;
end

%% Load results for all pairs
PLV_Rand_Variables_Pairs = cell(1,3);
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Stats Summary Across SS 8 vs 4 180320\';
    strVariablePath = [strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'];
    try
        Vars = load(strVariablePath);
        PLV_Rand_Variables_AcrossSS_Pairs{nPair} = Vars;
        clear Vars
    catch
    end
end

%% Plot results for each depth electrode for randomized values
strImageFolder = [strPaths.Sternberg,'Results\Multitaper Tapsmofrq 2 Rand PLV Across SS 8 vs 4 Scalp Layout Plots 180320\'];
mkdir(strImageFolder)
figLayout = Get_Figure_Scalp_Layout_Information();
% Initialize plot
Max_PLV_Pairs_AcrossSS_Table = zeros(size(nChannelPairs,1),6);
Max_PLV_Pairs_AcrossSS_Table = array2table(Max_PLV_Pairs_AcrossSS_Table);
Max_PLV_Pairs_AcrossSS_Table.Properties.VariableNames = {'Scalp','Depth','Max_PLV','freq_Max_PLV','AreaUnder_PLV','freqBand_AreaUnder_PLV'};
Max_PLV_Pairs_AcrossSS_Table.Scalp = cell(size(nChannelPairs,1),1);
Max_PLV_Pairs_AcrossSS_Table.Depth = cell(size(nChannelPairs,1),1);
Max_PLV_Pairs_AcrossSS_Table.freqBand_AreaUnder_PLV = cell(size(nChannelPairs,1),1);
freqAxis = 0.5:0.5:30;
% Plot layout for each depth electrode
for iChannel_Depth = 1:length(nChannelList_Depth)
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    % Depth channel
    nChannel_Depth = nChannelList_Depth(iChannel_Depth);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    % Scalp channels
    indDepthChannelPairs = find(nChannelPairs(:,2)==nChannel_Depth);
    nChannelList_Scalp = nChannelPairs(indDepthChannelPairs,1);
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp);
    % Plot results on layout
    for iScalp = 1:length(strScalpChannelNameList)
        %% Scalp channel number and label
        nChannel_Scalp = nChannelList_Scalp(iScalp);
        strScalpChannelName = strScalpChannelNameList{nChannel_Scalp};
        %% Pair
        nPair = find(nChannelPairs(:,1)==nChannel_Scalp&nChannelPairs(:,2)==nChannel_Depth);
        %% Subplot in layout
        indSubplot = find(strcmpi(figLayout(:,1),strScalpChannelName));
        nSubplot = figLayout{indSubplot,2};
        subplot(5,5,nSubplot)
        %% Plot results
        Vars = PLV_Rand_Variables_AcrossSS_Pairs{nPair};
        if(~isempty(Vars))
            % Frequency axis
            %             freqAxis = Vars.freqAxis;
            % Real PLV
            for iSS = 1:2
                Real_PLV_SS{iSS} = Vars.Real_PLV_AcrossSS{iSS};
                %             Real_PLV = PLV_Real_SS_Pair{nPair}{1};
            end
            
            %% Plot
            %             plot(freqAxis,abs(Real_PLV),strPlotColors{iSS})
            hold on
            shadedErrorBar(freqAxis,Vars.RandPrc{50},[Vars.RandPrc{95}-Vars.RandPrc{50};Vars.RandPrc{50}-Vars.RandPrc{5}],{'Color','b'});%0.5*[1,1,1]});
            %             plot(freqAxis,Vars.RandStats.Mean,'m')
            %             plot(freqAxis,Vars.RandStats.Mean+Vars.RandStats.Mean,'m')
            %             plot(freqAxis,Vars.RandStats.Mean-Vars.RandStats.Mean,'m')
            plot(freqAxis,Vars.RandPrc{95},'g')
            plot(freqAxis,abs(Real_PLV_SS{1})-abs(Real_PLV_SS{2}),strPlotColors{3})
            
            if(~(strcmpi(Vars.strPairLabels{1},strScalpChannelName)&&strcmpi(Vars.strPairLabels{2},strDepthChannelName)))
                error('Channel labels do not match')
            end
            
            %% Axis limits
            ylim([-0.5,0.5])
            xlim([1,30])
            set(gca,'XTick',[1,5:5:30])
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
            
            %% Chanel names
            Max_PLV_Pairs_AcrossSS_Table.Scalp{nPair} = strScalpChannelName;
            Max_PLV_Pairs_AcrossSS_Table.Depth{nPair} = strDepthChannelName;
            
            if 1
                %% Max PLV
                [valMax,indMax] = max(abs(abs(Real_PLV_SS{1})-abs(Real_PLV_SS{2})));
                Max_PLV_Pairs_AcrossSS_Table.Max_PLV(nPair) = valMax;
                Max_PLV_Pairs_AcrossSS_Table.freq_Max_PLV(nPair) = Vars.freqAxis(indMax);
                
                %% Area under curve
                thrPrc = 95;
                thrNeighborPts = 3;
                thrFreq = 0;
                Real_PLV = abs(abs(Real_PLV_SS{1})-abs(Real_PLV_SS{2})); % abs(Real_PLV);
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
                        stem(freqAxis(indMaxBand),Real_PLV(indMaxBand))
                        FreqBand = [6.5,18.5];%[6.5,18.5];
                        indFreqBand = freqAxis>=FreqBand(1)&freqAxis<=FreqBand(2);
                        AreaUnder_8_vs_4_Prc_SelectedBand = sum(Real_PLV(indFreqBand));
                    else
                        freqBand = [];
                        AreaUnder_PLV_Prc = 0;
                        AreaUnder_8_vs_4_Prc_SelectedBand = 0;
                    end
                    Max_PLV_Pairs_AcrossSS_Table.AreaUnder_PLV(nPair) = AreaUnder_PLV_Prc;
                    Max_PLV_Pairs_AcrossSS_Table.freqBand_AreaUnder_PLV{nPair} = freqBand;
                    Max_PLV_Pairs_AcrossSS_Table.freqBand_AreaUnder_PLV_Min{nPair} = min(freqBand);
                    Max_PLV_Pairs_AcrossSS_Table.freqBand_AreaUnder_PLV_Max{nPair} = max(freqBand);
                    Max_PLV_Pairs_AcrossSS_Table.AreaUnder_8_vs_4_Prc_SelectedBand(nPair) = AreaUnder_8_vs_4_Prc_SelectedBand;
                    
                    %                 pl = area(freqAxis(indMaxBand),x(indMaxBand));
                    title(sprintf('%s - %s - %.2f',strScalpChannelName,strDepthChannelName,AreaUnder_PLV_Prc))
                end
                
            end
            %%
        end
    end
    
    %%
    if 0
        %% Select band
        Max_PLV_Pairs_Table_SingleDepthChannel = Max_PLV_Pairs_Table(strcmpi(Max_PLV_Pairs_Table.Depth,strDepthChannelName),:);%&(~(strcmpi(Max_PLV_Pairs_Table.Scalp,'A1')|strcmpi(Max_PLV_Pairs_Table.Scalp,'A2'))),:);
        Max_PLV_Pairs_Table_SingleDepthChannel((strcmpi(Max_PLV_Pairs_Table_SingleDepthChannel.Scalp,'A1')|strcmpi(Max_PLV_Pairs_Table_SingleDepthChannel.Scalp,'A2')),:).AreaUnder_PLV = zeros(2,1)
        [~,indMaxArea] = max(Max_PLV_Pairs_Table_SingleDepthChannel.AreaUnder_PLV);
        freqBandSelected = Max_PLV_Pairs_Table_SingleDepthChannel.freqBand_AreaUnder_PLV{indMaxArea};
        
        %% Area under curve for the selected band
        for iScalp = 1:length(strScalpChannelNameList)
            %% Scalp channel number and label
            nChannel_Scalp = nChannelList_Scalp(iScalp);
            strScalpChannelName = strScalpChannelNameList{nChannel_Scalp};
            %% Pair
            nPair = find(nChannelPairs(:,1)==nChannel_Scalp&nChannelPairs(:,2)==nChannel_Depth);
            %% Plot results
            Vars = PLV_Rand_Variables_Pairs{iSS}{nPair};
            if(~isempty(Vars))
                % Frequency axis
                freqAxis = Vars.freqAxis;
                % Real PLV
                Real_PLV = Vars.Real_PLV;
                %             Real_PLV = PLV_Real_SS_Pair{nPair}{1};
                
                %% Area under curve for selected band
                thrPrc = 99;
                thrNeighborPts = 3;
                thrFreq = 2;
                Real_PLV = abs(Real_PLV);
                Thr_PLV = Vars.RandPrc{thrPrc};
                indFreqBandSelected = find(ismember(freqAxis,freqBandSelected));
                Real_PLV(indFreqBandSelected)-Thr_PLV(indFreqBandSelected)
                
                indAboveThr = Real_PLV>=Thr_PLV;
                cc = bwconncomp(indAboveThr);
                if(cc.NumObjects~=0)
                    % Threshold number of points
                    numPixels = cellfun(@numel,cc.PixelIdxList);
                    cc.PixelIdxList = cc.PixelIdxList(numPixels>=thrNeighborPts);
                    % Threshold lower bound for frequency
                    cc.PixelIdxList = cc.PixelIdxList(freqAxis(cellfun(@min,cc.PixelIdxList))>=thrFreq);
                    % Largest cluster
                    numPixels = cellfun(@numel,cc.PixelIdxList);
                    [~,indMaxCluster] = max(numPixels);
                    % Cluster with the highest area under
                    plv_all = cellfun(@(x) Real_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    thr_plv_all = cellfun(@(x) Thr_PLV(x),cc.PixelIdxList,'UniformOutput',0);
                    
                    if(~isempty(indMaxCluster))
                        indMaxBand = cc.PixelIdxList{indMaxCluster};
                        freqBand = freqAxis(indMaxBand);
                        AreaUnder_PLV_Prc = sum(Real_PLV(indMaxBand)-Thr_PLV(indMaxBand));
                        %                     stem(freqAxis(indMaxBand),Real_PLV(indMaxBand))
                    else
                        freqBand = [];
                        AreaUnder_PLV_Prc = 0;
                    end
                    Max_PLV_Pairs_Table.AreaUnder_PLV(nPair) = AreaUnder_PLV_Prc;
                    Max_PLV_Pairs_Table.freqBand_AreaUnder_PLV{nPair} = freqBand;
                    %                 pl = area(freqAxis(indMaxBand),x(indMaxBand));
                end
                
                %%
            end
        end
        
    end
    
    %% Save image
    %     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Depth_Channel_',num2str(iChannel_Depth,'%.3d'),'_',strDepthChannelName,'_mt_tsf_2_Rand_AcrossSS_8_vs_4','.png'],'png')
    %     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Depth_Channel_',num2str(iChannel_Depth,'%.3d'),'_',strDepthChannelName,'_mt_tsf_2_Rand_AcrossSS_8_vs_4','.fig'],'fig')
    close all
    
end

%%
Max_PLV_Pairs_AcrossSS_Table

%%
Max_PLV_Pairs_AcrossSS_Table_2 = Max_PLV_Pairs_AcrossSS_Table(Max_PLV_Pairs_AcrossSS_Table.Max_PLV~=0,:);
Max_PLV_Pairs_AcrossSS_Table_2 = Max_PLV_Pairs_AcrossSS_Table_2(Max_PLV_Pairs_AcrossSS_Table_2.AreaUnder_PLV~=0,:);
[~,indMax] = max(Max_PLV_Pairs_AcrossSS_Table_2.AreaUnder_PLV)
Max_PLV_Pairs_AcrossSS_Table_2(indMax,:)

[~,indSort] = sort(Max_PLV_Pairs_AcrossSS_Table_2.AreaUnder_PLV,'descend');
Max_PLV_Pairs_AcrossSS_Table_2(indSort,:)

%%
Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2 = Max_PLV_Pairs_AcrossSS_Table_2;%(ismember(Max_PLV_Pairs_AcrossSS_Table_2.Depth,'PHL1-2'),:);
Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2 = Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2(~ismember(Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2.Scalp,{'A1','A2'}),:);

% Max_PLV_Pairs_Table_SingleDepth_Temp =
sortrows(Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2,'AreaUnder_8_vs_4_Prc_SelectedBand','descend')
sortrows(Max_PLV_Pairs_Table_SingleDepth,'AreaUnder_PLV_Prc_SelectedBand','descend')

%%
datavectorToPlot = Max_PLV_Pairs_AcrossSS_Table_SingleDepth_2.AreaUnder_8_vs_4_Prc_SelectedBand;
datavectorToPlot(datavectorToPlot<0) = 0;

ChanlocsToPlot = EEGBipolarScalp.chanlocs(1:17);

figure(101)
topoplot(datavectorToPlot,ChanlocsToPlot,'maplimits',[0,2])
colorbar
title(sprintf('Mean PLV in %.1f-%.1f Hz during retention',FreqBand(1),FreqBand(2)))

%% Encoding maintenance comparison 
%% Set size comparison
nSS_Ran = [3,3];
nNumberOfPermutations = 1000;
nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
nTrialNumbersAll_AcrossSS = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
nRandomTrialNumberList_AcrossTime_Pairs = cell(1,size(nChannelPairs,1));
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
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

strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint ss8 Trial Numbers 180320\';
mkdir(strVariableFolder)
save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')

%% Extract 2 seconds of stimulus
dataBipolarScalp_Stim_SS = cell(1,4);
for iSS = 1:4
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolarScalp_SS{iSS}.fsample];
    dataBipolarScalp_Stim_SS{iSS} = ft_selectdata(cfg,dataBipolarScalp_SS{iSS});
end
%% Extract 1 second of stimulus
dataBipolarScalp_Stim_2_SS = cell(1,4);
for iSS = 1:4
    cfg = [];
    cfg.latency = [-4,-3-1/dataBipolarScalp_SS{iSS}.fsample];
    dataBipolarScalp_Stim_2_SS{iSS} = ft_selectdata(cfg,dataBipolarScalp_SS{iSS});
end

%% Pairs to run analysis for
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','Pz')));

nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);

%%
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint 180320\';
tStartRand = cputime;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,2);
    iSS_AcrossSS = 1;
    iSS = nSS_Ran(iSS_AcrossSS);
    cfg = [];
    cfg.channel = [nChannel_1,nChannel_2];
    data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
    
    iSS_AcrossSS = 2;
    iSS = nSS_Ran(iSS_AcrossSS);
    cfg = [];
    cfg.channel = [nChannel_1,nChannel_2];
    data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Stim_SS{iSS});
    
    %% Randomly mix trials between set sizes
    for nRand = 1:201 % nNumberOfPermutations+1
        
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
        for iSS = 1:2
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.tapsmofrq   = 2;
            cfg.pad         = 2;
            cfg.foi         = 0.5:0.5:30;
            cfg.channel     = 1:2;
            Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            %% Coherence
            % cfg             = [];
            % cfg.method      = 'coh';
            % cfg.complex     = 'complex';
            % Coh_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS}{nRand});
            
        end
    end
    clear Freq_SS
    
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}{1}.label';
    freqAxis = PLV_SS{iSS}{1}.freq;
    for iSS = 1:2
        for nRand = 1:201%nNumberOfPermutations+1
            PLV_AcrossSS{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
        end
    end
    clear PLV_SS
    
    %% Save results
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
        'PLV_AcrossSS','strPairLabels','freqAxis','-v7.3')
    
    %% Save statistical summary
    %% Real values
    nRand = 1;
    for iSS = 1:2
        Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
    end
    
    %% Null distribution
    for iSS = 1:2
        RandDist_PLV_All_AcrossSS{iSS} = [];
        for nRand = 1:length(PLV_AcrossSS{iSS})-1
            RandDist_PLV_All_AcrossSS{iSS}(nRand,:) = PLV_AcrossSS{iSS}{nRand+1};
        end
    end
    clear PLV_AcrossSS
    
    PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
    RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
    
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
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint Stats 180320\';
    mkdir(strStatsFolder)
    save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'],...
        'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
    clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end

tStopRand = cputime-tStartRand;

%% Load results / PLV for all pairs / Randomized trials / Maintenance vs Encoding
strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint Stats 180320\';
PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,4);
for iSS = 1:4
    PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS} = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair} = Vars;
            clear Vars
        catch
        end
    end
end

%% Display results for a single pair
nPair = 589;
iSS = 3;
Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair};

fig = figure;
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1})-abs(Vars.Real_PLV_AcrossSS{2}),'LineWidth',2)
hold on
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1}))
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{2}))

% figure
% plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1})-abs(Vars.Real_PLV_AcrossSS{2}),'LineWidth',2)
% hold on
plot(Vars.freqAxis,Vars.RandPrc{95})
plot(Vars.freqAxis,Vars.RandPrc{5})

% legend('Difference','95%')
legend('Difference','Maint','Enc','95%','5%')

strImageFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Enc vs Maint Single Pair Images\';
saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,1)},'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,2)},'_mt_tsf_2_Enc_Maint','.png'],'png')
saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,1)},'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,2)},'_mt_tsf_2_Enc_Maint','.fig'],'fig')
% Vars.RandStats

%%
if 0
%% Set size comparison / high vs low
nSS_Ran = [4,1];
nNumberOfPermutations = 1000;
nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
nTrialNumbersAll_AcrossSS_HL = [1:nNumberOfSetSizes_ToCompare(1),-1*(1:nNumberOfSetSizes_ToCompare(2))];
nRandomTrialNumberList_AcrossSS_Pairs_HL = cell(1,size(nChannelPairs,1));
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    nRandomTrialNumberList = zeros(nNumberOfPermutations,sum(nNumberOfSetSizes_ToCompare));
    for nRand = 1:nNumberOfPermutations
        cond = 1;
        while(cond)
            indRand = randperm(sum(nNumberOfSetSizes_ToCompare));
            cond = ~isempty(find((indRand-(1:length(indRand)))==0,1));
        end
        nRandomTrialNumberList(nRand,:) = nTrialNumbersAll_AcrossSS_HL(indRand);
    end
    nRandomTrialNumberList = [nTrialNumbersAll_AcrossSS_HL;nRandomTrialNumberList];
    nRandomTrialNumberList_AcrossSS_Pairs_HL{nPair} = nRandomTrialNumberList;
end

strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Across SS 68 vs 4 Trial Numbers 180320\';
mkdir(strVariableFolder)
save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_AcrossSS_Pairs_HL','-v7.3')

%% Channel pairs to run analysis for
thrPLV = 0.2;
iSS = 3;
flagRunAnalysis_Pair = zeros(size(nChannelPairs,1),1);
for iPair = 1:size(nChannelPairs,1)
    flagRunAnalysis_Pair(iPair) = length(find(abs(PLV_Real_SS_Pair{iPair}{iSS})>thrPLV))>3;
end
nPairs_ToRun = find(flagRunAnalysis_Pair);
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
nPairs_ToRun = intersect(nPairs_ToRun,ind);

%% Analysis for each pair of electrodes - Randomization
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Values Across SS 68 vs 4 180320\';
tStartRand = cputime;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,2);
    for iSS_AcrossSS = 1:length(nSS_Ran)
        iSS = nSS_Ran(iSS_AcrossSS);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
    end
    
    %% Randomly mix trials between set sizes
    for nRand = 1:nNumberOfPermutations+1
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        iTrialList = nRandomTrialNumberList_AcrossSS_Pairs_HL{nPair}(nRand,:);
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
        for iSS = 1:2
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.tapsmofrq   = 2;
            cfg.pad         = 2;
            cfg.foi         = 0.5:0.5:30;
            cfg.channel     = 1:2;
            Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            %% Coherence
            % cfg             = [];
            % cfg.method      = 'coh';
            % cfg.complex     = 'complex';
            % Coh_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS}{nRand});
            
        end
    end
    clear Freq_SS
    
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}{1}.label';
    freqAxis = PLV_SS{iSS}{1}.freq;
    for iSS = 1:2
        for nRand = 1:5%nNumberOfPermutations+1
            PLV_AcrossSS{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
        end
    end
    clear PLV_SS
    
    %% Save results
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_68_vs_4','.mat'],...
        'PLV_AcrossSS','strPairLabels','freqAxis','-v7.3')
    
    %% Save statistical summary
    %% Real values
    nRand = 1;
    for iSS = 1:2
        Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
    end
    
    %% Null distribution
    for iSS = 1:2
        RandDist_PLV_All_AcrossSS{iSS} = [];
        for nRand = 1:length(PLV_AcrossSS{iSS})-1
            RandDist_PLV_All_AcrossSS{iSS}(nRand,:) = PLV_AcrossSS{iSS}{nRand+1};
        end
    end
    clear PLV_AcrossSS
    
    PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
    RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
    
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
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Stats Summary Across SS 68 vs 4 180320\';
    mkdir(strStatsFolder)
    save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_68_vs_4','.mat'],...
        'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
    clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end

tStopRand = cputime-tStartRand;
%193


end
%% Extract 2 seconds of stimulus
dataBipolarScalp_Stim_SS = cell(1,4);
for iSS = 1:4
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolarScalp_SS{iSS}.fsample];
    dataBipolarScalp_Stim_SS{iSS} = ft_selectdata(cfg,dataBipolarScalp_SS{iSS});
end

%% Run PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Run_Real_PLV_Stim
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    % ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','1-2')));
    % nPairs_ToRun = intersect(nPairs_ToRun,ind);
    
    %% Analysis for each pair of electrodes - Real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Stim Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Stim_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.tapsmofrq   = 2;
            cfg.foi         = 0.5:0.5:30;
            cfg.channel     = 1:2;
            cfg.feedback    = 'no';
            Freq            = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            %% Coherence
            cfg             = [];
            cfg.method      = 'coh';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
        end
        clear Freq cfg
        
        %% Store results for the channel pair
        strPairLabels = PLV_SS{iSS}.label';
        freqAxis = PLV_SS{iSS}.freq;
        for iSS = 1:4
            PLV_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
            Coh_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
        end
        TimeInterval = [-5,-3-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'],'PLV_SS','strPairLabels','freqAxis','Coh_SS','TimeInterval','-v7.3')
        clear PLV_SS Coh_SS strPairLabels freqAxis TimeInterval
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
        
    end
    
end

%% Load results / for PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Load_Real_PLV_Stim % Run results for patient from folder if 0    
    %% Load results for all pairs
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Stim Tapsmofrq 2 Values 180405\'];
    PLV_Real_Stim_SS_Pair = cell(1,size(nChannelPairs,1));
    Coh_Real_Stim_SS_Pair = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Real_Stim_SS_Pair{nPair} = Vars.PLV_SS;
            Coh_Real_Stim_SS_Pair{nPair} = Vars.Coh_SS;
            freqAxis = Vars.freqAxis;
            clear Vars nChannel_1 nChannel_2 strVariablePath
        catch
            clear Vars nChannel_1 nChannel_2 strVariablePath
        end
    end
    strPairLabels_Pair = dataBipolarScalp_Ret_SS{iSS}.label(nChannelPairs);
    
    %% Save real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Stim Tapsmofrq 2 Values 180405\'];
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Stim_mt_tsf_2.mat'],'PLV_Real_Stim_SS_Pair','strPairLabels_Pair','freqAxis','Coh_Real_Stim_SS_Pair','-v7.3')
    
else
    %% Load real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Stim Tapsmofrq 2 Values 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Stim_mt_tsf_2.mat'])
    
end

%% Plot the real values on for each depth electrode - real values
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Stim Tapsmofrq 2 Scalp Layouts 180405\'];
mkdir(strImageFolder)
figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = 1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs(:,2)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs(nPairs_WithDepthChannel,1);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs(:,1)==nChannelScalp&nChannelPairs(:,2)==nChannel_Depth);
        
        %% Subplot number for scalp channel
        nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
        nSubplot = figLayout{nEntryInTable,2};
        subplot(5,5,nSubplot)
        
        %% Plot
        for iSS = 1:3 %[1,4,2,3]
            plot(freqAxis,abs(PLV_Real_Stim_SS_Pair{nPairs_WithDepthScalpChannel}{iSS}),strPlotColors{iSS},'LineWidth',2)
            hold on
            xlim([1,30])
            ylim([0,0.8])
            set(gca,'XTick',[1,5:5:30])
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
        end
        %%
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Stim_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Stim_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.fig'],'fig')
    close all
    
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList
end
%

%%
if 1
   %% Encoding maintenance comparison 
%% Set size comparison
nSS_Ran = [3,3];
nNumberOfPermutations = 1000;
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

strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint ss8 Trial Numbers 180320\';
mkdir(strVariableFolder)
save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_AcrossTime_Pairs','-v7.3')

%% Extract 2 seconds of stimulus
dataBipolarScalp_Stim_SS = cell(1,4);
for iSS = 1:4
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolarScalp_SS{iSS}.fsample];
    dataBipolarScalp_Stim_SS{iSS} = ft_selectdata(cfg,dataBipolarScalp_SS{iSS});
end

%% Pairs to run analysis for
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','Pz')));

nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);

%%
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint 180320\';
tStartRand = cputime;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,2);
    iSS_AcrossSS = 1;
    iSS = nSS_Ran(iSS_AcrossSS);
    cfg = [];
    cfg.channel = [nChannel_1,nChannel_2];
    data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
    
    iSS_AcrossSS = 2;
    iSS = nSS_Ran(iSS_AcrossSS);
    cfg = [];
    cfg.channel = [nChannel_1,nChannel_2];
    data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_Stim_SS{iSS});
    
    %% Randomly mix trials between set sizes
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
        for iSS = 1:2
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.tapsmofrq   = 2;
            cfg.pad         = 2;
            cfg.foi         = 0.5:0.5:30;
            cfg.channel     = 1:2;
            Freq_SS{iSS}    = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            PLV_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS});
            
            %% Coherence
            % cfg             = [];
            % cfg.method      = 'coh';
            % cfg.complex     = 'complex';
            % Coh_SS{iSS}{nRand} = ft_connectivityanalysis(cfg,Freq_SS{iSS}{nRand});
            
        end
    end
    clear Freq_SS
    
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}{1}.label';
    freqAxis = PLV_SS{iSS}{1}.freq;
    for iSS = 1:2
        for nRand = 1:nNumberOfPermutations+1
            PLV_AcrossSS{iSS}{nRand} = squeeze(PLV_SS{iSS}{nRand}.plvspctrm(1,2,:))';
        end
    end
    clear PLV_SS
    
    %% Save results
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
        'PLV_AcrossSS','strPairLabels','freqAxis','-v7.3')
    
    %% Save statistical summary
    %% Real values
    nRand = 1;
    for iSS = 1:2
        Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
    end
    
    %% Null distribution
    for iSS = 1:2
        RandDist_PLV_All_AcrossSS{iSS} = [];
        for nRand = 1:length(PLV_AcrossSS{iSS})-1
            RandDist_PLV_All_AcrossSS{iSS}(nRand,:) = PLV_AcrossSS{iSS}{nRand+1};
        end
    end
    clear PLV_AcrossSS
    
    PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
    RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
    
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
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint Stats 180320\';
    mkdir(strStatsFolder)
    save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'],...
        'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','-v7.3')
    clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end

tStopRand = cputime-tStartRand;

%% Load results / PLV for all pairs / Randomized trials / Maintenance vs Encoding
strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint Stats 180320\';
PLV_Rand_Enc_Maint_Variables_Pair_SS = cell(1,4);
for iSS = 1:4
    PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS} = cell(1,size(nChannelPairs,1));
    for nPair = 1:size(nChannelPairs,1)
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariablePath = [strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair} = Vars;
            clear Vars
        catch
        end
    end
end

%% Display results for a single pair
nPair = nPairs_ToRun;
iSS = 3;
Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS}{nPair};

fig = figure;
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1})-abs(Vars.Real_PLV_AcrossSS{2}),'LineWidth',2)
hold on
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1}))
plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{2}))

% figure
% plot(Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1})-abs(Vars.Real_PLV_AcrossSS{2}),'LineWidth',2)
% hold on
plot(Vars.freqAxis,Vars.RandPrc{95})
plot(Vars.freqAxis,Vars.RandPrc{5})

% legend('Difference','95%')
legend('Difference','Maint','Enc','95%','5%')
strImageFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Enc vs Maint Single Pair Images\';
saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,1)},'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,2)},'_mt_tsf_2_Enc_Maint','.png'],'png')
saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,1)},'_',dataBipolarScalp_SS{1}.label{nChannelPairs(nPair,2)},'_mt_tsf_2_Enc_Maint','.fig'],'fig')

% Vars.RandStats 

%%
nPair = nPairs_ToRun;
iSS = 3;
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Enc Maint 180320\';
iSS_ToLoad = 2;
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        strVariableFolder = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS_ToLoad+2,'%.1d'),'.mat'];

Vars = load(strVariableFolder);

PLV_Maint = cell2mat(Vars.PLV_AcrossSS{1}');
PLV_Maint_Real = PLV_Maint(1,:);
PLV_Maint = PLV_Maint(2:end,:);
PLV_Enc = cell2mat(Vars.PLV_AcrossSS{2}');
PLV_Enc_Real = PLV_Enc(1,:);
PLV_Enc = PLV_Enc(2:end,:);

PLV_Diff_Rand = PLV_Maint-PLV_Enc;
PLV_Diff_Real = PLV_Maint_Real-PLV_Enc_Real;

%%
clear timelockSS1 timelockSS2
iSS = 1;
freqAxis = 0.5:0.5:30;
for tr = 1:size(PLV_Diff_Real,1)*1000
    timelockSS1.trial(tr,1,:)   = abs(PLV_Diff_Real(1,:));
    timelockSS1.time            = freqAxis;
end
timelockSS1.label       = '1';
timelockSS1.fsample     = 1;
timelockSS1.dimord      = 'rpt_chan_time';
iSS = 2;
for tr = 1:size(PLV_Diff_Rand,1)
    timelockSS2.trial(tr,1,:)   = abs(PLV_Diff_Rand(tr,:));
    timelockSS2.time            = freqAxis;
end
timelockSS2.label       = '2';
timelockSS2.fsample     = 1;
timelockSS2.dimord      = 'rpt_chan_time';

%%
cfg                             = [];
cfg.method                      = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic                   = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
% evaluate the effect at the sample level
cfg.correctm                    = 'cluster';
cfg.clusteralpha                = 0.05;         % alpha level of the sample-specific test statistic that
% will be used for thresholding
cfg.clusterstatistic            = 'maxsum'; % test statistic that will be evaluated under the
% permutation distribution.
cfg.minnbchan                   = 0;               % minimum number of neighborhood channels that is
cfg.tail                        = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail                 = 0;
cfg.alpha                       = 0.05;               % alpha level of the permutation test
cfg.numrandomization            = 100;      % number of draws from the permutation distribution
% Design matrix
design = zeros(1,size(timelockSS1.trial,1) + size(timelockSS2.trial,1));
design(1,1:size(timelockSS1.trial,1)) = 1;
design(1,(size(timelockSS1.trial,1)+1):(size(timelockSS1.trial,1) + size(timelockSS2.trial,1)))= 2;
cfg.design = design;             % design matrix
cfg.ivar  = 1;                   % number or list with indices indicating the independent variable(s)
cfg.channel       = {'1','2'};     % cell-array with selected channel labels
cfg.latency       = [0.5 30];       % time interval over which the experimental
% conditions must be compared (in seconds)
for ch = 1
    neighbours(ch).label = '';
    neighbours(ch).neighblabel = {''};
end
cfg.neighbours    = neighbours;  % the neighbours specify for each sensor with
[stat] = ft_timelockstatistics(cfg,timelockSS1,timelockSS2);

end

%%

if 1
    strPaths.Results = [strPaths.Sternberg,'Results\'];

%% Run PLV using cross spectrum 
if flag_Run_Load.Run_Real_PLV_CrossSpectrum
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','Pz')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    %% Analysis for each pair of electrodes - Real values
%     strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f Values 180405\'];
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f LowRe Values 180503\']; % Low Re
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Time-frequency analysis
            cfg = [];
            cfg.output     = 'powandcsd'; % fourier pow
            cfg.method     = 'mtmconvol';
            cfg.foi        = 5:1:30; % 5:0.5:30
            cfg.t_ftimwin  = 10./cfg.foi;
            cfg.tapsmofrq  = 0.2*cfg.foi;
            cfg.toi        = -6:0.25:2; % -6:0.1:2
            cfg.keeptrials = 'yes';
            CrossFreq = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV using cross spectrum
            PLV_SS{iSS} = squeeze(CrossFreq.crsspctrm);
            PLV_SS{iSS} = squeeze(nanmean(PLV_SS{iSS}./abs(PLV_SS{iSS}),1));
            
        end
        
        %% Store results for the channel pair
        strPairLabels = CrossFreq.labelcmb';
        freqAxis = CrossFreq.freq;
        timeAxis = CrossFreq.time;
        TimeInterval = [-6,2];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_0_2f.mat'],'PLV_SS','strPairLabels','freqAxis','timeAxis','cfgFreq','TimeInterval','-v7.3')
        clear PLV_SS strPairLabels freqAxis timeAxis TimeInterval CrossFreq
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
    end
end

%% Plot the real values on for each depth electrode - real values
iSS_ToPlot = 3;
% strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f Values 180405\'];
strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f LowRe Values 180503\']; % Low Re
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f Scalp Layouts 180405\'];
mkdir(strImageFolder)
figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelList_Depth),'PHL1-2')))
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs(:,2)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs(nPairs_WithDepthChannel,1);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs(:,1)==nChannelScalp&nChannelPairs(:,2)==nChannel_Depth);
        nPair = nPairs_WithDepthScalpChannel;
        
        %% Plot for pairs with results
        try
            %% Load variables
            Vars = load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannelScalp},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_Depth},'_mt_tsf_0_2f.mat']);
            PLV = Vars.PLV_SS{iSS_ToPlot};
            freqAxis = Vars.freqAxis;
            timeAxis = Vars.timeAxis;
            
            %% Subplot number for scalp channel
            nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
            nSubplot = figLayout{nEntryInTable,2};
            subplot(5,5,nSubplot)
            
            %% Plot
            imagesc(timeAxis,freqAxis,abs(PLV),[0,0.6])
            hold on
            set(gca,'YTick',[1,5:5:30])
            set(gca,'YDir','normal')
            title(sprintf('%s - %s / ss%d',strScalpChannelName,strDepthChannelName,2*iSS_ToPlot+2))
            colorbar
            colormap jet
        catch
        end
    end
    %% Save image
%     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.png'],'png')
%     saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.fig'],'fig')
%     close all
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList Vars
end

%% Plot the real values on for each depth electrode - real values
iSS_ToPlot = 3;
strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum Tapsmofrq 0_2f Values 180405\'];
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Cross Spectrum In Time 0_2f Scalp Layouts 180405\'];
figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelList_Depth),'PHL1-2')))
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs(:,2)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs(nPairs_WithDepthChannel,1);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs(:,1)==nChannelScalp&nChannelPairs(:,2)==nChannel_Depth);
        nPair = nPairs_WithDepthScalpChannel;
        
        %% Plot for pairs with results
        try
            %% Load variables
            Vars = load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannelScalp},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_Depth},'_mt_tsf_0_2f.mat']);
            
            freqAxis = Vars.freqAxis;
            timeAxis = Vars.timeAxis;
            
            %%
            FreqBand = [9,12];
            [~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
            [~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
            
            %% Subplot number for scalp channel
            nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
            nSubplot = figLayout{nEntryInTable,2};
            subplot(5,5,nSubplot)
            
            %% Plot
            for iSS = 1:3
                PLV = Vars.PLV_SS{iSS};
                PLV = abs(PLV);
                PLV_In_Band = mean(PLV(indFreq1:indFreq2,:));
                plot(timeAxis,PLV_In_Band,strPlotColors{iSS})
                hold on
                title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
            end
        catch
        end
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.fig'],'fig')
%     close all
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList Vars
end

%% Set size comparison
nSS_Ran = [3,1];
nNumberOfPermutations = 1000;
nNumberOfSetSizes_ToCompare = nNumberOfTrialsSetSizes(nSS_Ran);
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Multitaper Tapsmofrq 2 Rand PLV Across SS 8 vs 4 Trial Numbers 180320\';
load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_AcrossSS_Pairs')

%% Pairs to run analysis for
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','Pz')));

nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);

%% Analysis for each pair of electrodes - Randomization
% strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f Rand Values Across SS 8 vs 4 180427\';
strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f Rand LowRe Values Across SS 8 vs 4 180503\';  % Low Re
tStartRand = cputime;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,2);
    for iSS_AcrossSS = 1:length(nSS_Ran)
        iSS = nSS_Ran(iSS_AcrossSS);
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS_AcrossSS} = ft_preprocessing(cfg,dataBipolarScalp_SS{iSS});
    end        
    
    %% Randomly mix trials between set sizes
    clear PLV_SS
    for nRand = 1:nNumberOfPermutations+1
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        iTrialList = nRandomTrialNumberList_AcrossSS_Pairs{nPair}(nRand,:);
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
        for iSS = 1:2
            %% Time-frequency analysis
            cfg = [];
            cfg.output     = 'powandcsd'; % fourier pow
            cfg.method     = 'mtmconvol';
            cfg.foi        = 5:1:30; % 5:0.5:30;
            cfg.t_ftimwin  = 10./cfg.foi;
            cfg.tapsmofrq  = 0.2*cfg.foi;
            cfg.toi        = -6:0.25:2; % -6:0.1:2;
            cfg.keeptrials = 'yes';
            CrossFreq = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV using cross spectrum
            PLV_SS{iSS}{nRand} = squeeze(CrossFreq.crsspctrm);
            PLV_SS{iSS}{nRand} = squeeze(nanmean(PLV_SS{iSS}{nRand}./abs(PLV_SS{iSS}{nRand}),1));
            
        end
    end
    
    %% Store results for the channel pair
    strPairLabels = CrossFreq.labelcmb';
    freqAxis = CrossFreq.freq;
    timeAxis = CrossFreq.time;
    TimeInterval = [-6,2];
    PLV_AcrossSS = PLV_SS;
    
    %% Save results
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat'],...
        'PLV_AcrossSS','strPairLabels','freqAxis','timeAxis','TimeInterval','-v7.3')
    
    %% Save statistical summary
    %% Real values
    nRand = 1;
    for iSS = 1:2
        Real_PLV_AcrossSS{iSS} = PLV_AcrossSS{iSS}{nRand};
    end
    
    %% Null distribution
    for iSS = 1:2
        RandDist_PLV_All_AcrossSS{iSS} = [];
        for nRand = 1:length(PLV_AcrossSS{iSS})-1
            RandDist_PLV_All_AcrossSS{iSS}(nRand,:,:) = PLV_AcrossSS{iSS}{nRand+1};
        end
    end
%     clear PLV_AcrossSS
    
    PLV_Rand_Dist_SingleSS = abs(RandDist_PLV_All_AcrossSS{1})-abs(RandDist_PLV_All_AcrossSS{2});
    RandDist_PLV_All_Im = imag(RandDist_PLV_All_AcrossSS{1})-imag(RandDist_PLV_All_AcrossSS{2});
    
    %% Percentiles
    prcList = [1,5,50,95,99];
    RandPrc = {};
    RandPrc_Im = {};
    for iPrc = 1:length(prcList)
        for ii = 1:size(PLV_Rand_Dist_SingleSS,2)
            for jj = 1:size(PLV_Rand_Dist_SingleSS,3)
                temp = PLV_Rand_Dist_SingleSS(:,ii,jj);
                val  = prctile(temp,prcList(iPrc));
                RandPrc{prcList(iPrc)}(ii,jj) = val;
                temp = RandDist_PLV_All_Im(:,ii,jj);
                val  = prctile(temp,prcList(iPrc));
                RandPrc_Im{prcList(iPrc)}(ii,jj) = val;
            end
        end
    end
    
    %% Min max
    for ii = 1:size(PLV_Rand_Dist_SingleSS,2)
        for jj = 1:size(PLV_Rand_Dist_SingleSS,3)
    RandStats.Min(ii,jj) = min(PLV_Rand_Dist_SingleSS(:,ii,jj));
    RandStats.Max(ii,jj) = max(PLV_Rand_Dist_SingleSS(:,ii,jj));
    RandStats.Mean(ii,jj) = mean(PLV_Rand_Dist_SingleSS(:,ii,jj));
    RandStats.Std(ii,jj) = std(PLV_Rand_Dist_SingleSS(:,ii,jj));
    RandStats.NoPermutations(ii,jj) = size(PLV_Rand_Dist_SingleSS(:,ii,jj),1);
    
    RandStats_Im.Min(ii,jj) = min(RandDist_PLV_All_Im(:,ii,jj));
    RandStats_Im.Max(ii,jj) = max(RandDist_PLV_All_Im(:,ii,jj));
    RandStats_Im.Mean(ii,jj) = mean(RandDist_PLV_All_Im(:,ii,jj));
    RandStats_Im.Std(ii,jj) = std(RandDist_PLV_All_Im(:,ii,jj));
    RandStats_Im.NoPermutations(ii,jj) = size(RandDist_PLV_All_Im(:,ii,jj),1);
        end
    end
    %% Save results
%     strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f Stats Summary Across SS 8 vs 4 180320\';
    strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f LowRe tats Summary Across SS 8 vs 4 180503\'; % LowRe
    mkdir(strStatsFolder)
    save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat'],...
        'Real_PLV_AcrossSS','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','timeAxis','TimeInterval','-v7.3')
    clear Real_PLV_AcrossSS RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis
    
end

tStopRand = cputime-tStartRand;

%% Plot results
nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PHL1-2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','Pz')));

nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);

%% 
% strVariableFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f Rand Values Across SS 8 vs 4 180427\';
strStatsFolder = 'I:\MATLAB Codes\Sternberg Task\Results\Phase Coherence FieldTrip\Rand PLV Cross Spectrum Tapsmofrq 0_2f Stats Summary Across SS 8 vs 4 180320\';
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    %% Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    %% Load results
% Vars = load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Set_Size_',num2str(2*iSS+2,'%.1d'),'.mat']);
Vars = load([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Pair_',num2str(nPair,'%.4d'),'_Scalp_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_Depth_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_Summary_Across_SS_8_vs_4','.mat']);
%%
thrPrc = 95;
Thr_PLV = Vars.RandPrc{thrPrc};
%
figure
imagesc(Vars.timeAxis,Vars.freqAxis,abs(Vars.Real_PLV_AcrossSS{1}))
set(gca,'YDir','normal')
hold on
temp = (Vars.Real_PLV_AcrossSS{1}-Vars.Real_PLV_AcrossSS{2})>Thr_PLV;
for ii = 1:size(temp,1)
    for jj = 1:size(temp,2)
        if(temp(ii,jj)>1)
            plot(Vars.freqAxis(ii),Vars.timeAxis(jj),'ro','MarkerSize',10)
        end
    end
end
figure,imagesc(temp)
set(gca,'YDir','normal')


end




end

%%

%% PLV filter in time
FrequencyBandLimits = [9,12];
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Filt In Time 0_2f Scalp Layouts 180405\'];
mkdir(strImageFolder)
figLayout = Get_Figure_Scalp_Layout_Information();

for iChannelDepth = 13 % 1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs(:,2)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs(nPairs_WithDepthChannel,1);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs(:,1)==nChannelScalp&nChannelPairs(:,2)==nChannel_Depth);
        nPair = nPairs_WithDepthScalpChannel;
        
        %%
        Freq1 = FrequencyBandLimits(1);
        Freq2 = FrequencyBandLimits(2);
        
        cfg             = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [Freq1,Freq2];
        cfg.channel     = [nChannelScalp,nChannel_Depth];
        
        dataFilteredSS = cell(4,1);
        for iSS = 1:3
            dataFilteredSS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_SS{iSS});
        end
        
        %%
        plv_SS = cell(4,1);
        for iSS = 1:3
            eegData = NaN(2,size(dataFilteredSS{iSS}.trial{1},2),length(dataFilteredSS{iSS}.trial));
            for nTrial = 1:length(dataFilteredSS{iSS}.trial)
                for nChannel = 1:2
                    Temp = dataFilteredSS{iSS}.trial{nTrial}(nChannel,:);
                    eegData(nChannel,:,nTrial) = Temp;
                end
            end
            fs = dataFilteredSS{iSS}.fsample;
            time_axis_PLV = dataFilteredSS{iSS}.time{1};
            plv = pn_eegPLV_mod_170301(eegData);
            
            plv_SS{iSS} = plv;
            clear plv
        end
        
        %% Subplot number for scalp channel
            nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
            nSubplot = figLayout{nEntryInTable,2};
            subplot(5,5,nSubplot)
            
            %% Plot
            tSmooth = 0.5;
            nSmooth = tSmooth/0.005+1;
            for iSS = 1:3
                plot(time_axis_PLV,smooth(plv_SS{iSS}(:,1,2),nSmooth),strPlotColors{iSS})
                hold on
            end
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
    end
     %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_CS_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_0_2f_ss',num2str(2*iSS_ToPlot+2),'.fig'],'fig')
%     close all
end

%% PLV between depth-depth contacts
flag_Run_Load.Run_Real_PLV_Ret_Depth = 1;
flag_Run_Load.Rand_PLV_Depth_TrialNumberList = 1;
flag_Run_Load.Run_Rand_PLV_Depth_Ret = 1;

%% Pairs of depth channels for phase coherence analysis
nChannelPairs_Depth = [];
for iChannel1 = 1:length(nChannelList_Depth)
    for iChannel2 = 1:length(nChannelList_Depth)
        nChannelPairs_Depth = [nChannelPairs_Depth;[nChannelList_Depth(iChannel1),nChannelList_Depth(iChannel2)]];
    end
end

%% Run PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Run_Real_PLV_Ret_Depth 
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs_Depth,1); % all pairs
    ind = find(nChannelPairs_Depth(:,1)<nChannelPairs_Depth(:,2));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Depth(:,1))','1-2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    %% Analysis for each pair of electrodes - Real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Ret Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.foi         = 0.5:0.5:30;
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1:2;
            Freq            = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            %% Coherence
            cfg             = [];
            cfg.method      = 'coh';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
        end
        clear Freq cfg
        
        %% Store results for the channel pair
        strPairLabels = PLV_SS{iSS}.label';
        freqAxis = PLV_SS{iSS}.freq;
        for iSS = 1:4
            PLV_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
            Coh_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
        end
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'],'PLV_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','Coh_SS','-v7.3')
        clear PLV_SS Coh_SS strPairLabels freqAxis cfgFreq
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
    end
end

%% Load results / PLV for all pairs / 2 seconds of retention
if flag_Run_Load.Load_Real_PLV_Ret_Depth % Run results for patient from folder if 0    
    %% Load results for all pairs
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Ret Tapsmofrq 2 Values 180405\'];
    PLV_Real_Depth_Ret_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    Coh_Real_Depth_Ret_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    PLV_Direction = zeros(1,size(nChannelPairs_Depth,1));
    for nPair = 1:size(nChannelPairs_Depth,1)
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
        if(exist(strVariablePath,'file')==2) % File for pair exists
            try
                Vars = load(strVariablePath);
                PLV_Real_Depth_Ret_SS_Pair{nPair} = Vars.PLV_SS;
                Coh_Real_Depth_Ret_SS_Pair{nPair} = Vars.Coh_SS;
                freqAxis = Vars.freqAxis;
                PLV_Direction(nPair) = 1;
                clear Vars nChannel_1 nChannel_2 strVariablePath
            catch
                clear Vars nChannel_1 nChannel_2 strVariablePath
            end
        else
            nPair_Reverse = find((nChannelPairs_Depth(:,1)==nChannel_2&nChannelPairs_Depth(:,2)==nChannel_1));
            nChannel_1 = nChannelPairs_Depth(nPair_Reverse,1);
            nChannel_2 = nChannelPairs_Depth(nPair_Reverse,2);
            strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair_Reverse,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
            if(exist(strVariablePath,'file')==2) % File for reverse pair exists
                try
                    Vars = load(strVariablePath);
                    PLV_Real_Depth_Ret_SS_Pair{nPair} = Vars.PLV_SS;
                    Coh_Real_Depth_Ret_SS_Pair{nPair} = Vars.Coh_SS;
                    freqAxis = Vars.freqAxis;
                    PLV_Direction(nPair) = -1;
                    clear Vars nChannel_1 nChannel_2 strVariablePath
                catch
                    clear Vars nChannel_1 nChannel_2 strVariablePath
                end
            end
        end
    end
    strPairLabels_Pair = dataBipolarScalp_Ret_SS{iSS}.label(nChannelPairs_Depth);
    
    %% Save real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Ret Tapsmofrq 2 Values 180405\'];
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'],'PLV_Real_Depth_Ret_SS_Pair','PLV_Direction','strPairLabels_Pair','freqAxis','Coh_Real_Depth_Ret_SS_Pair','-v7.3')
    
else
    %% Load real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Ret Tapsmofrq 2 Values 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'])
    
end

%% Plot the real values on for each depth electrode - real values
if flag_Run_Load.Plot_Real_PLV_Ret_Depth
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Abs Ret Tapsmofrq 2 Scalp Layouts 180405\'];
mkdir(strImageFolder)
% figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = 1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs_Depth(:,1)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs_Depth(nPairs_WithDepthChannel,2);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs_Depth(:,2)==nChannelScalp&nChannelPairs_Depth(:,1)==nChannel_Depth);
        
        %% Subplot number for scalp channel
        nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
%         nSubplot = figLayout{nEntryInTable,2};
        nSubplot = iScalp;
        subplot(8,6,nSubplot)
        
        %% Plot
        for iSS = 1:3 %[1,4,2,3]
            try
            plot(freqAxis,abs(PLV_Real_Depth_Ret_SS_Pair{nPairs_WithDepthScalpChannel}{iSS}),strPlotColors{iSS},'LineWidth',2)
            hold on
            xlim([1,30])
            ylim([0,1])
            set(gca,'XTick',[1,5:5:30])
            catch
            end
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
        end
        %%
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.fig'],'fig')
    close all
    
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList
end
end

%% Random trial numbers for each set size separately
if flag_Run_Load.Rand_PLV_Depth_TrialNumberList
    nNumberOfPermutations = 1000;
    nRandomTrialNumberList_Depth_Pair_SS = cell(1,4);
    for iSS = 1:4
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
            nRandomTrialNumberList_Depth_Pair_SS{iSS}{nPair} = nRandomTrialNumberList;
        end
    end
    %% Save variables
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Depth Trial Numbers 180405\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'],'nRandomTrialNumberList_Depth_Pair_SS','nNumberOfPermutations','-v7.3')
else
    %% Load variables
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Depth Trial Numbers 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Randomized_Trial_Numbers.mat'])
end

%% Run PLV for all pairs / Randomized trials / 2 seconds of retention
if flag_Run_Load.Run_Rand_PLV_Depth_Ret
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Depth Ret Tapsmofrq 2 Values 180405\'];
    strStatsFolder = [strPaths.Results,'PLV FieldTrip\Rand PLV Depth Ret Tapsmofrq 2 Stats 180405\'];
    mkdir(strVariableFolder)
    mkdir(strStatsFolder)

    %% Select pairs for analysis
    thrPLV = 0.2;
    iSS = 3;
    flagRunAnalysis_Pair = zeros(size(nChannelPairs,1),1);
    for iPair = 1:size(nChannelPairs,1)
        flagRunAnalysis_Pair(iPair) = ~isempty(find(abs(PLV_Real_Depth_Ret_SS_Pair{iPair}{iSS})>thrPLV,1));
    end
    nPairs_ToRun = find(flagRunAnalysis_Pair);
    ind = find(nChannelPairs_Depth(:,1)<nChannelPairs_Depth(:,2));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','1-2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    
    %% Analysis for each pair of electrodes - Randomization
    iSS_ToCompute = 3;
    tStartRand = cputime;
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = iSS_ToCompute
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Ret_SS{iSS});
        end
        
        %% Create datasets with randomized trials
        data_SinglePair_Rand_SS = cell(1,4);
        nChannel_ToChange = 2;
        for iSS = iSS_ToCompute
            for nRand = 1:nNumberOfPermutations+1
                data_SinglePair_Rand_SS{iSS}{nRand} = data_SinglePair_SS{iSS};
                for nTrial = 1:length(data_SinglePair_SS{iSS}.trial)
                    nTrialInRand = nRandomTrialNumberList_Depth_Pair_SS{iSS}{nPair}(nRand,nTrial);
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
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Convert cell to double
        for iSS = 1:4
            try
                PLV_Rand_SS{iSS} = cell2mat(PLV_Rand_SS{iSS}');
            catch
                PLV_Rand_SS{iSS} = [];
            end
        end
        
        %% Save values
        iSS = iSS_ToCompute;
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_ss',num2str(2*iSS+2),'.mat'],...
            'PLV_Rand_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
        
        %% Save statistical summary
        
        %% Real values
        nRand = 1;
        Real_PLV = PLV_Rand_SS{iSS}(1,:);
        
        %% Null distribution
        PLV_Rand_Dist_SingleSS = PLV_Rand_SS{iSS}(2:end,:);
        
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
        save([strStatsFolder,'Patient_',num2str(nPatient,'%.2d'),'_Rand_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2_Stats_ss',num2str(2*iSS+2),'.mat'],...
            'Real_PLV','RandPrc','RandPrc_Im','RandStats','RandStats_Im','strPairLabels','freqAxis','cfgFreq','TimeInterval','-v7.3')
        clear Real_PLV RandPrc RandPrc_Im RandStats RandStats_Im strPairLabels freqAxis cfgFreq TimeInterval PLV_Rand_SS
        
    end
    
    tStopRand = cputime-tStartRand;
    
end






%% Run PLV for all pairs / 2 seconds of stimulus
flag_Run_Load.Run_Real_PLV_Stim_Depth = 1;
if flag_Run_Load.Run_Real_PLV_Stim_Depth
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs_Depth,1); % all pairs
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Depth(:,1))','PHL1-2')));
    ind2 = find(cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Depth(:,2))','PHL1-2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    %% Analysis for each pair of electrodes - Real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Stim_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.foi         = 0.5:0.5:30;
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1:2;
            Freq            = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            %% Coherence
            cfg             = [];
            cfg.method      = 'coh';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
        end
        clear Freq cfg
        
        %% Store results for the channel pair
        strPairLabels = PLV_SS{iSS}.label';
        freqAxis = PLV_SS{iSS}.freq;
        for iSS = 1:4
            PLV_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
            Coh_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
        end
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'],'PLV_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','Coh_SS','-v7.3')
        clear PLV_SS Coh_SS strPairLabels freqAxis cfgFreq
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
    end
end

%% Load results / PLV for all pairs / 2 seconds of stimulus
flag_Run_Load.Load_Real_PLV_Stim_Depth = 1;
if flag_Run_Load.Load_Real_PLV_Stim_Depth % Run results for patient from folder if 0    
    %% Load results for all pairs
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim Tapsmofrq 2 Values 180405\'];
    PLV_Real_Depth_Stim_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    Coh_Real_Depth_Stim_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    for nPair = 1:size(nChannelPairs_Depth,1)
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Real_Depth_Stim_SS_Pair{nPair} = Vars.PLV_SS;
            Coh_Real_Depth_Stim_SS_Pair{nPair} = Vars.Coh_SS;
            freqAxis = Vars.freqAxis;
            clear Vars nChannel_1 nChannel_2 strVariablePath
        catch
            clear Vars nChannel_1 nChannel_2 strVariablePath
        end
    end
    strPairLabels_Pair = dataBipolarScalp_Ret_SS{iSS}.label(nChannelPairs_Depth);
    
    %% Save real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim Tapsmofrq 2 Values 180405\'];
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'],'PLV_Real_Depth_Stim_SS_Pair','strPairLabels_Pair','freqAxis','Coh_Real_Depth_Stim_SS_Pair','-v7.3')
    
else
    %% Load real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim Tapsmofrq 2 Values 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'])
    
end

%% Plot the real values on for each depth electrode - real values
flag_Run_Load.Plot_Real_PLV_Stim_Depth = 1;
if flag_Run_Load.Plot_Real_PLV_Stim_Depth
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Abs Stim Tapsmofrq 2 Scalp Layouts 180405\'];
mkdir(strImageFolder)
% figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = 13%1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs_Depth(:,1)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs_Depth(nPairs_WithDepthChannel,2);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs_Depth(:,2)==nChannelScalp&nChannelPairs_Depth(:,1)==nChannel_Depth);
        
        %% Subplot number for scalp channel
        nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
%         nSubplot = figLayout{nEntryInTable,2};
        nSubplot = iScalp;
        subplot(7,6,nSubplot)
        
        %% Plot
        for iSS = 1:3 %[1,4,2,3]
            try
            plot(freqAxis,abs(PLV_Real_Depth_Stim_SS_Pair{nPairs_WithDepthScalpChannel}{iSS}),strPlotColors{iSS},'LineWidth',2)
            hold on
            xlim([1,30])
            ylim([0,1])
            set(gca,'XTick',[1,5:5:30])
            catch
            end
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
        end
        %%
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.fig'],'fig')
    close all
    
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList
end
end


%% Run PLV for all pairs / 2 seconds of stimulus
flag_Run_Load.Run_Real_PLV_Stim_2_Depth = 1;
if flag_Run_Load.Run_Real_PLV_Stim_2_Depth
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs_Depth,1); % all pairs
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Depth(:,1))','PHL1-2')));
    ind2 = find(cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Depth(:,2))','PHL1-2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
    
    %% Analysis for each pair of electrodes - Real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim 2 Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
        
        %% Channel numbers
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,4);
        for iSS = 1:4
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolarScalp_Stim_2_SS{iSS});
        end
        
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
        
        %% Phase coherence
        for iSS = 1:4
            %% Frequency analysis
            cfg             = [];
            cfg.method      = 'mtmfft';
            cfg.taper       = 'dpss';
            cfg.output      = 'fourier';
            cfg.foi         = 1:1:30;
            cfg.tapsmofrq   = 2;
            cfg.channel     = 1:2;
            Freq            = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
            
            cfgFreq = cfg;
            
            %% PLV analysis
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
            %% Coherence
            cfg             = [];
            cfg.method      = 'coh';
            cfg.complex     = 'complex';
            cfg.feedback    = 'no';
            Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
            
        end
        clear Freq cfg
        
        %% Store results for the channel pair
        strPairLabels = PLV_SS{iSS}.label';
        freqAxis = PLV_SS{iSS}.freq;
        for iSS = 1:4
            PLV_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
            Coh_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
        end
        TimeInterval = [-2,-1/dataBipolarScalp_SS{iSS}.fsample];
        
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'],'PLV_SS','strPairLabels','freqAxis','cfgFreq','TimeInterval','Coh_SS','-v7.3')
        clear PLV_SS Coh_SS strPairLabels freqAxis cfgFreq
        
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
    end
end

%% Load results / PLV for all pairs / 2 seconds of stimulus
flag_Run_Load.Load_Real_PLV_Stim_2_Depth = 1;
if flag_Run_Load.Load_Real_PLV_Stim_2_Depth % Run results for patient from folder if 0    
    %% Load results for all pairs
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim 2 Tapsmofrq 2 Values 180405\'];
    PLV_Real_Depth_Stim_2_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    Coh_Real_Depth_Stim_2_SS_Pair = cell(1,size(nChannelPairs_Depth,1));
    for nPair = 1:size(nChannelPairs_Depth,1)
        nChannel_1 = nChannelPairs_Depth(nPair,1);
        nChannel_2 = nChannelPairs_Depth(nPair,2);
        strVariablePath = [strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_1},'_',dataBipolarScalp_Ret_SS{iSS}.label{nChannel_2},'_mt_tsf_2.mat'];
        try
            Vars = load(strVariablePath);
            PLV_Real_Depth_Stim_2_SS_Pair{nPair} = Vars.PLV_SS;
            Coh_Real_Depth_Stim_2_SS_Pair{nPair} = Vars.Coh_SS;
            freqAxis = Vars.freqAxis;
            clear Vars nChannel_1 nChannel_2 strVariablePath
        catch
            clear Vars nChannel_1 nChannel_2 strVariablePath
        end
    end
    strPairLabels_Pair = dataBipolarScalp_Ret_SS{iSS}.label(nChannelPairs_Depth);
    
    %% Save real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim 2 Tapsmofrq 2 Values 180405\'];
    mkdir(strVariableFolder)
    save([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'],'PLV_Real_Depth_Stim_2_SS_Pair','strPairLabels_Pair','freqAxis','Coh_Real_Depth_Stim_2_SS_Pair','-v7.3')
    
else
    %% Load real values
    strVariableFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Stim Tapsmofrq 2 Values 180405\'];
    load([strVariableFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Depth_Ret_mt_tsf_2.mat'])
    
end

%% Plot the real values on for each depth electrode - real values
flag_Run_Load.Plot_Real_PLV_Stim_2_Depth = 1;
if flag_Run_Load.Plot_Real_PLV_Stim_2_Depth
strImageFolder = [strPaths.Results,'PLV FieldTrip\Real PLV Depth Abs Stim 2 Tapsmofrq 2 Scalp Layouts 180405\'];
mkdir(strImageFolder)
% figLayout = Get_Figure_Scalp_Layout_Information();
for iChannelDepth = 13%1:length(nChannelList_Depth)
    %% Figure
    fig = figure;
    Make_Plot_Secondary_Display_Fullscreen_Gca;
    %% Channel pairs with the selected depth channel
    nChannel_Depth = nChannelList_Depth(iChannelDepth);
    nPairs_WithDepthChannel = nChannelPairs_Depth(:,1)==nChannel_Depth;
    nChannelList_Scalp_WithDepthChannel = nChannelPairs_Depth(nPairs_WithDepthChannel,2);
    strDepthChannelName = strChannelNameList{nChannel_Depth};
    strScalpChannelNameList = strChannelNameList(nChannelList_Scalp_WithDepthChannel);
    
    %% For each scalp channel
    for iScalp = 1:length(strScalpChannelNameList)
        %% Pair number
        nChannelScalp = nChannelList_Scalp_WithDepthChannel(iScalp);
        strScalpChannelName = strScalpChannelNameList{iScalp};
        nPairs_WithDepthScalpChannel = find(nChannelPairs_Depth(:,2)==nChannelScalp&nChannelPairs_Depth(:,1)==nChannel_Depth);
        
        %% Subplot number for scalp channel
        nEntryInTable = find(strcmpi(figLayout(:,1),strScalpChannelNameList{iScalp}));
%         nSubplot = figLayout{nEntryInTable,2};
        nSubplot = iScalp;
        subplot(7,6,nSubplot)
        
        %% Plot
        for iSS = 1:3 %[1,4,2,3]
            try
            plot(freqAxis,abs(PLV_Real_Depth_Stim_2_SS_Pair{nPairs_WithDepthScalpChannel}{iSS}),strPlotColors{iSS},'LineWidth',2)
            hold on
            xlim([1,30])
            ylim([0,1])
            set(gca,'XTick',[1,5:5:30])
            catch
            end
            title(sprintf('%s - %s',strScalpChannelName,strDepthChannelName))
        end
        %%
    end
    %% Save image
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.png'],'png')
    saveas(fig,[strImageFolder,'Patient_',num2str(nPatient,'%.2d'),'_Real_PLV_Ret_Layout_Depth_',num2str(iChannelDepth,'%.2d'),'_',strDepthChannelName,'_mt_tsf_2.fig'],'fig')
    close all
    
    %%
    clear nChannel_Depth nPairs_WithDepthChannel strDepthChannelName strScalpChannelNameList
end
end


%%
%%
%%
%%
%%





