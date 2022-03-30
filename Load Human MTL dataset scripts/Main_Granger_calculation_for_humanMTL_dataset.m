%% Close all figures, clear variables and command window
close all
clear
clc

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
% strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20191126\';
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


% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Parameters
analysis_type = 'correct trials'; %'incorrect trials'; %or 'correct trials'
Extra_Analysis = 0;
% Granger Causality calculation between a hippocampal contact of a depth 
%electrode and a contact in the cortex of a depth electrode (outer contacts)
flag_Hipp_GC_Ctx_Depth = 1; 

%% Load the data
strDataPath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Dataset_2020_MTL_iEEG_scalpEEG\';
strDataPath_Macro = [strDataPath,'Macro\'];
strDataPath_Scalp = [strDataPath,'Scalp\'];
nSubjects = [1:9];
strPat_IDs = {'WC','SW','IA','HJ','PL','AM','HF','MA','MS'};
nSubject_ID = [28 22 19 30 33 13 23 29 16];
% specify the patient's data directory
selected_subj =  9;
strPatientDataPath_Macro = sprintf('%s%d %s\\',strDataPath_Macro,nSubject_ID(selected_subj),strPat_IDs{selected_subj});
strPatientDataPath_Scalp = sprintf('%s%d %s\\',strDataPath_Scalp,nSubject_ID(selected_subj),strPat_IDs{selected_subj});

[dataBipolar, dataBipolar_Scalp, TrialInformationTable] = ...
    Load_Dataset_from_human_MTL(strPatientDataPath_Macro,strPatientDataPath_Scalp);

cd(strPaths.Main)

%% Reject Artifactual trials
[dataBipolar_clean, TrialInformationTable_iEEG_clean ] = ...
    Get_Only_Clean_Trials_FieldTrip(dataBipolar, TrialInformationTable.Macro);

[dataBipolar_Scalp_clean, TrialInformationTable_scalp_EEG_clean ] = ...
    Get_Only_Clean_Trials_FieldTrip(dataBipolar_Scalp, TrialInformationTable.Scalp);


%% Rereference
% Load anatomical locations for the current patient
strAnatLocationsPath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Dataset_2020_MTL_iEEG_scalpEEG\Anatomical Locations\';
strFilePrefix = 'Anatomical_Locations_Patient_';
strSubjAnatLocationFile = sprintf('%s%s%d.mat',strAnatLocationsPath,strFilePrefix,nSubject_ID(selected_subj));
load(strSubjAnatLocationFile);

% Find white matter contacts
nWM_contacts = find(strcmp(AnatomicalLocations.Anatomical_Location,'no_label_found'))
strWM_labels = AnatomicalLocations.Labels(nWM_contacts);
hipp_chans_WM_contacts = find(contains(strWM_labels,'H'));
refChannel_iEEG = nWM_contacts(hipp_chans_WM_contacts(1));
% Rereferencing iEEG
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =   refChannel_iEEG;%white matter contacts referencing
cfg.refmethod     =  'avg';
dataBipolar_reref       = ft_preprocessing(cfg,dataBipolar_clean);
dataBipolar_reref.label = strrep(dataBipolar.label,'m','');

% Rereferencing scalp EEG
A1_chan = find(contains(dataBipolar_Scalp_clean.label,'A1'));
A2_chan = find(contains(dataBipolar_Scalp_clean.label,'A2'));


cfg               = []
cfg.channel       = {'all'};
cfg.reref         =  'yes';
cfg.refchannel    =  [A1_chan, A2_chan];
cfg.refmethod     =  'avg';
dataBipolarScalp_reref  =  ft_preprocessing(cfg,dataBipolar_Scalp_clean);


%% Apply montage
[montage,chans_to_be_bipolar,dataBip,Bip_chans] = GetDataInBipolar(dataBipolar_reref,selected_subj,flag_Hipp_GC_Ctx_Depth);

%% Downsample to 200 Hz
cfg = [];
cfg.resamplefs = 80;%200;
dataBipolarResampled = ft_resampledata(cfg,dataBip);
dataBipolarScalpResampled = ft_resampledata(cfg,dataBipolarScalp_reref);

%% Append scalp and iEEG data
if ~flag_Hipp_GC_Ctx_Depth
    cfg = [];
    dataBipolarResampled.time = dataBipolarScalpResampled.time;
    data_bipolar_appended = ft_appenddata(cfg,dataBipolarResampled,dataBipolarScalpResampled)
else
    data_bipolar_appended = dataBipolarResampled;
end

%% macro data
% Select only correct trials
switch analysis_type
    case 'correct trials'
        [dataBipolar_SS,TrialInformationTableCorrect] = Get_Only_Correct_Trials_FieldTrip(data_bipolar_appended,TrialInformationTable_iEEG_clean);
    case 'incorrect trials'
        [dataBipolar_SS, TrialInformationTableCorrect] = Get_Only_Incorrect_Trials_FieldTrip(data_bipolar_appended, TrialInformationTable_iEEG_clean )
end
%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
TrialInformationTableCorrect.SetSize = TrialInformationTableCorrect.Setsize;
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTableCorrect);

%% Select the latencies for every task period for each task period
[dataBipolar_Ret_SS{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
[dataBipolar_Ret_SS{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
[dataBipolar_Ret_SS{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance

%% Channel pairs to run
if ~flag_Hipp_GC_Ctx_Depth
    %To be configured for  every combination of electrode pairs you have
    clear nChannelPairs;
    
    nScalpChans = length(dataBipolar_Scalp.label);
    
    nChannelPairs = [];
    
    for i=1:length(Bip_chans)
        nChannelPairs = [nChannelPairs;[i+size(montage.labelold,2)+zeros(nScalpChans,1),(1:nScalpChans)'+Bip_chans(end) ]];
    end
    strChannelNameList = dataBipolar_Ret_SS{1}{1}.label;
    
    hipp_chans= find(contains(strChannelNameList(nChannelPairs(:,1)),'H'));
    hipp_chanPairs = nChannelPairs(hipp_chans,:);
    
    %Validate that the channel pairs are correct
    disp(strChannelNameList(hipp_chanPairs))
else
    % Channel pairs to run for iEEG Granger (e.g.: AHL2-3 - AHL6-7)
    %To be configured for  every combination of electrode pairs you have
    clear nChannelPairs;
    
    nScalpChans = length(dataBipolar_Scalp.label);
    
    nChannelPairs = [];
    
    for i=1:length(Bip_chans)/2
        nChannelPairs = [nChannelPairs; [i+size(montage.labelold,2)+zeros(floor(length(Bip_chans)/2),1),(1:length(Bip_chans)/2)'+floor(length(Bip_chans)/2) + size(montage.labelold,2)]];
    end
    strChannelNameList = dataBipolar_Ret_SS{1}{1}.label;
    
    hipp_chans= find(contains(strChannelNameList(nChannelPairs(:,1)),'H'));
    hipp_chanPairs = nChannelPairs(hipp_chans,:);
    
    %Validate that the channel pairs are correct
    disp(strChannelNameList(hipp_chanPairs))
end
%%
if Extra_Analysis
%% Granger session Wise
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\SubjectDeltaGranger.mat');
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\SubjectAccuracy.mat');
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\SignificantBarMaint_all_subjects.mat');

currentSubject = 44;%nSubject_ID(selected_subj);
indDeltaGranger = find(SubjectDeltaGrangerMaint(:,1) == currentSubject);
DeltaGranger_currentSubj = SubjectDeltaGrangerMaint(indDeltaGranger,2);
SubjectAccuracy_currentSubj = SubjectAccuracy(indDeltaGranger,2);

[TrialInfo, dataBipolar_Session, CorrectResponseRate_Session, dataBipolar_Ret_SS_session ] = ...
    getDataInSessions(dataBipolar_SS,TrialInformationTable_SS,TrialInformationTable_iEEG_clean)
nTrial_per_Session = cellfun(@numel,TrialInfo)./size(TrialInfo{1},2);


[gdata_session, DeltaGranger_Maint] = getGranger_session_Wise(dataBipolar_Ret_SS_session,hipp_chanPairs,strChannelNameList)

[maxGranger_Maint_sessions{selected_subj} ] = ...
    Plot_Granger_sessionWise(gdata_session,DeltaGranger_Maint,CorrectResponseRate_Session,...
    SignificantBars_GrangerScalp_Maint,selected_subj,DeltaGranger_currentSubj,SubjectAccuracy_currentSubj,nTrial_per_Session)


Granger_SessionWise{selected_subj}.maxGranger_MaintSession = maxGranger_Maint_sessions{selected_subj};
Granger_SessionWise{selected_subj}.CorrectReponseSession = CorrectResponseRate_Session;
Granger_SessionWise{selected_subj}.SignificantBarsMaint = SignificantBars_GrangerScalp_Maint{selected_subj};
Granger_SessionWise{selected_subj}.DeltaGranger_Maint = DeltaGranger_Maint;
Granger_SessionWise{selected_subj}.DeltaGranger_allSessions_x_Accuracy = [DeltaGranger_currentSubj SubjectAccuracy_currentSubj];

%%
figure;
plot(SubjectDeltaGrangerMaint(:,2),SubjectAccuracy(:,2),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor','k');
    hold on;
set(gca,'box','off','FontSize',16);

model = fitlm(SubjectDeltaGrangerMaint(3:15,2),SubjectAccuracy(3:15,2),'RobustOpts','on')
h = plot(model);
hLegend = findobj(gcf, 'Type', 'Legend');
title('');

set(hLegend,'Visible','off')
% set('legend','visible','off')
xlabel('?Granger (%)');
ylabel('Correct Response Accuracy (%)')
[rho, pval] = corr(SubjectAccuracy(3:15,2),SubjectDeltaGrangerMaint(3:15,2));

annotation('textbox',...
    [0.564285714285714 0.789476193882171 0.34196427605514 0.138095234689259],...
    'String',{'R =0.35, p = 0.03','rho = 0.4, p = 0.15'},...
    'FontSize',14,...
    'EdgeColor','none');

ylim([75 100])
%% Power calculation - trial wise during maintenance

        %% Correct Trials
%         TrialInformationTable_iEEG_clean.SetSize = TrialInformationTable_iEEG_clean.Setsize; %only for scalp Patients
        pID = 14%selected_subj;
        chanPair = 19;
        
        [dataBipolar_Maint TrialInfoTable_Temp] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable_iEEG_clean);
        [dataBipolar_Maint_SS TrialInfoTable_Temp] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_Maint,TrialInfoTable_Temp);
        dataBipolar_Ret_SS{1} = Extract_task_period_Data('fix',dataBipolar_Maint_SS);
        dataBipolar_Ret_SS{2} = Extract_task_period_Data('maint',dataBipolar_Maint_SS);
        switch subjectID
            case 37
                channel_cmb = {'PHL1-PHL2','PL1'};        % S37 grid patient
            case 42
                channel_cmb = {'_AHL2-_AHL3','_GL_C2'};   % S42 grid patient
            case 38
                channel_cmb = {'PHR1-PHR3','TLSL1'};      % S38 grid patient
        end
%         cmb1 = strChannelNameList{hipp_chanPairs(chanPair,1)};
%         cmb2 = strChannelNameList{hipp_chanPairs(chanPair,2)};
%         channel_cmb ={cmb1,cmb2}
        iSS = 4;
        grangerChannelPair = channel_cmb;
        foi_1 = 9:18;
        foi_2 = 9:18;
        % ?
        [Fr_Maint_Hipp_Correct_scalp{pID}.theta, Fr_Fix_Hipp_Correct_scalp{pID}.theta]= ...
            getPower_trialWise(dataBipolar_Ret_SS{2},dataBipolar_Ret_SS{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[4:8])
        % ?
         [Fr_Maint_Hipp_Correct_scalp{pID}.alpha, Fr_Fix_Hipp_Correct_scalp{pID}.alpha]= ...
            getPower_trialWise(dataBipolar_Ret_SS{2},dataBipolar_Ret_SS{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[8:12])
        % ?-?
         [Fr_Maint_Hipp_Correct_scalp{pID}.thetaalpha, Fr_Fix_Hipp_Correct_scalp{pID}.thetaalpha]= ...
            getPower_trialWise(dataBipolar_Ret_SS{2},dataBipolar_Ret_SS{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[4:12])
        
        [Fr_Maint_Hipp_Correct{pID}.beta, Fr_Fix_Hipp_Correct{pID}.beta]= ...
            getPower_trialWise(dataBipolar_Ret_SS{2},dataBipolar_Ret_SS{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[9:18])
        
        %% Incorrect Trials
        foi_1 = 9:18%1:100;
        foi_2 = 9:18%1:100;
        [dataBipolar_Incor TrialInfoTable_Temp] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable_iEEG_clean);
        [dataBipolar_Maint_SS TrialInfoTable_Temp] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_Incor,TrialInfoTable_Temp);
        dataBipolar_Fix_Incor{1} = Extract_task_period_Data('fix',dataBipolar_Maint_SS);
        dataBipolar_Maint_Incor{2} = Extract_task_period_Data('maint',dataBipolar_Maint_SS);

        Fr_Maint_Hipp_Incorrect_Scalp{pID}.theta = getPower_trialWise(dataBipolar_Maint_Incor{2},dataBipolar_Fix_Incor{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[4:8])
        
         Fr_Maint_Hipp_Incorrect_Scalp{pID}.alpha = getPower_trialWise(dataBipolar_Maint_Incor{2},dataBipolar_Fix_Incor{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[8:12])
        
         Fr_Maint_Hipp_Incorrect_Scalp{pID}.thetaalpha = getPower_trialWise(dataBipolar_Maint_Incor{2},dataBipolar_Fix_Incor{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[4:12])
        
         Fr_Maint_Hipp_Incorrect{pID}.beta = getPower_trialWise(dataBipolar_Maint_Incor{2},dataBipolar_Fix_Incor{1},...
            iSS,grangerChannelPair,foi_1,foi_2,[9:18])

        %% Visualize results for scalp patients
%         figure;
        subplot(3,3,selected_subj)
        for i = 129:length(grangerChannelPair)
            [modelCorr{i,selected_subj} modelIncor{i,selected_subj}] = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{selected_subj}.beta{i},Fr_Maint_Hipp_Incorrect{selected_subj}.beta{i},'beta',grangerChannelPair(i,:))
        end
end
%% Calculate Granger - Power Correlation for all channel pairs
clear gdata PowCorr;
tStart = tic;
nSet_Size = size(Set_Sizes,1);
Fix_Data = dataBipolar_Ret_SS{1};
Enc_Data = dataBipolar_Ret_SS{2};
Maint_Data = dataBipolar_Ret_SS{3};
for iPair = 1:size(hipp_chanPairs,1)
    fprintf('Calculating Granger for pair %d for all set sizes\n',iPair)
    for iSS = 4%:nSet_Size%1:nSet_Size
        freq                   = [];
        freq.freqcfg           = [];
        freq.freqcfg.method    = 'mtmfft';
%         freq.freqcfg.foi       = [1:1:100];%[4:1:100];
        freq.freqcfg.output    = 'fourier';
        freq.freqcfg.taper     = 'hanning';
        freq.freqcfg.pad       = 20;
        freq.freqcfg.tapsmofrq = 2;
        Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS});
        Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS});
        Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS});
        
        grangercfg = [];
        grangercfg.method  = 'granger';
        grangercfg.granger.conditional = 'no';
        grangercfg.granger.sfmethod = 'bivariate';
        grngChan1 = strChannelNameList{hipp_chanPairs(iPair,1)};
        grngChan2 = strChannelNameList{hipp_chanPairs(iPair,2)};
        grangercfg.channelcmb = {grngChan1,grngChan2};
        gdata.Maint{iPair,iSS}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
        gdata.Enc{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Encoding_freq);
        gdata.Fix{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Fixation_freq);
        
        
%         cfg = [];
%         cfg.method = 'powcorr';
%         cfg.channelcmb = grangercfg.channelcmb;
%         PowCorr.Fix{iPair,iSS}   =   ft_connectivityanalysis(cfg,Fixation_freq);
%         PowCorr.Enc{iPair,iSS}   =   ft_connectivityanalysis(cfg,Encoding_freq);
%         PowCorr.Maint{iPair,iSS} =   ft_connectivityanalysis(cfg,Maintenance_freq);
        
        %% Delta Granger Incorrect Trials
        Dgranger_enc = gdata.Enc{iPair,iSS}.grangerspctrm(2,:)-gdata.Enc{iPair,iSS}.grangerspctrm(1,:);
        Dgranger_maint = gdata.Maint{iPair,iSS}.grangerspctrm(1,:)-gdata.Maint{iPair,iSS}.grangerspctrm(2,:);
        
    end
    tStop = toc(tStart);
    fprintf('\n\n\n\n\n\n\n\n\n\n Pair %d %f seconds elapsed \n',iPair,tStop)
end
%% Save results
strSavePath = 'F:/Vasileios/Task Analysis/Data/Analysis Data/Granger for human MTL dataset/';
strPatSavePath = [strSavePath,'Patient ', num2str(nSubject_ID(selected_subj)),' ' strPat_IDs{selected_subj},'\'];
mkdir(strPatSavePath)
save([strPatSavePath,'Granger_depth_scalp_all_channel_pairs_all_set_sizes'],'gdata','-v7.3');


%% Save results for Power correlation
strSavePath = 'F:/Vasileios/Task Analysis/Data/Analysis Data/Power Correlation for human MTL dataset/';
mkdir(strSavePath);
strPatSavePath = [strSavePath,'Patient ', num2str(nSubject_ID(selected_subj)),' ' strPat_IDs{selected_subj},'\'];
mkdir(strPatSavePath)
switch analysis_type
    case 'correct trials'
      save([strPatSavePath,'Power_correlation_EEG_hipp_all_pairs_correct_trials'],'gdata','-v7.3');
     case 'incorrect trials'
        save([strPatSavePath,'Power_correlation_EEG_hipp_all_pairs_inccorrect_trials'],'gdata','-v7.3');
end

%%
Granger_spectra_incorrectTrials{selected_subj}.Enc.grangerspctrm = gdata.Enc{iPair,iSS}.grangerspctrm
Granger_spectra_incorrectTrials{selected_subj}.Maint.grangerspctrm = gdata.Maint{iPair,iSS}.grangerspctrm;
%% DGranger save values for selected subj
frequency_point = [17];
Dgranger_IncorrectTrials_EEG{selected_subj}.Enc = Dgranger_enc(frequency_point);
Dgranger_IncorrectTrials_EEG{selected_subj}.Maint = Dgranger_maint(frequency_point)
disp(Dgranger_IncorrectTrials_EEG{selected_subj})



%% Visualize the results (granger spectra in scalp eeg arrangement)
%Select Bipolar channel to plot
hipp_chan_to_plot = find(contains(strChannelNameList(Bip_chans),'H'));
Bip_chans_to_plot = Bip_chans(hipp_chan_to_plot);
%Select set size to plot
iSS = 4;
% colors
dark_Green = [0.02 0.47 0.04];
orange       = [1 0.45 0.45];
Colors     = {dark_Green,'g','r',orange};

freq_ax    =  gdata.Enc{1,4}.freq%gdata.Fix{1,1}.freq;

for i = 1:numel(Bip_chans_to_plot)
    indBipChan = find(hipp_chanPairs==Bip_chans_to_plot(i))
    figure('units','normalized','outerposition',[0 0 1 1])
    ha = tight_subplot(5,9,[.07 .03],[.1 .01],[.01 .01])
    set(gcf,'color','white')
    figLayout = Get_Figure_Scalp_Layout_Information();
    for j = 1:numel(indBipChan)
        chanPair = indBipChan(j);
        Scalp_channel = strChannelNameList{hipp_chanPairs(chanPair,2)};
        indSubplot = find(strcmp(figLayout(:,1),Scalp_channel));
        if ~isempty(indSubplot)
            nSubplot(j) = figLayout{indSubplot,2};
            axes(ha(nSubplot(j)));
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
            hold on;
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
            xlim([4 30])
            ylabel(Scalp_channel)
            set(gca,'box','off')
            %             ylim([0 0.2])
        end
    end
    for k = 1:size(ha,1)
        
        if ~ismember(k,nSubplot)
            set(ha(k),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    suptitle(sprintf('Depth Electrode %s - Scalp Granger',dataBipolar_Ret_SS{1}{iSS}.label{Bip_chans_to_plot(i)}))
end


%% Visualize the results (iEEG granger spectra)
if flag_Hipp_GC_Ctx_Depth
    %Select Bipolar channel to plot
    hipp_chan_to_plot = find(contains(strChannelNameList(Bip_chans),'H'));
    Bip_chans_to_plot = Bip_chans(hipp_chan_to_plot(1:end/2));
    
    %Select set size to plot
    iSS = 4;
    % colors
    
    
    light_blue = [0.30,0.75,0.93];
    light_red       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',light_red};
    freq_ax    =  gdata.Enc{17,4}.freq;%gdata.Fix{1,1}.freq;
    
    
    for i = 1:numel(Bip_chans_to_plot)
        indBipChan = find(hipp_chanPairs==Bip_chans_to_plot(i))
        figure('units','normalized','outerposition',[0 0 1 1])
        ha = tight_subplot(4,4,[.07 .05],[.1 .01],[.03 .01])
        set(gcf,'color','white')
        %     figLayout = Get_Figure_Scalp_Layout_Information();
        for j = 1:numel(indBipChan)
            chanPair = indBipChan(j);
            Cx_iEEG_chan = strChannelNameList{hipp_chanPairs(chanPair,2)};
            %         indSubplot = find(strcmp(figLayout(:,1),Scalp_channel));
            axes(ha(j));
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
            hold on;
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
            xlim([4 20])
            ylabel(Cx_iEEG_chan)
            set(gca,'box','off')
            %             ylim([0 0.2])
            
        end
        %     for k = 1:size(ha,1)
        %
        %         if ~ismember(k,nSubplot)
        %             set(ha(k),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        %         end
        %
        %     end
        suptitle(sprintf('Depth Electrode %s - iEEG Granger',dataBipolar_Ret_SS{1}{iSS}.label{Bip_chans_to_plot(i)}))
    end
end

%% Visualize results for Power Correlation
chanPair = 1;
iSS = 4;
dark_Green = [0.02 0.47 0.04];
Colors = {'k',dark_Green,'r'};
freq_ax    =  PowCorr.Fix{1,4}.freq;%gdata.Fix{1,1}.freq;
figure;
semilogx(freq_ax,PowCorr.Fix{chanPair,iSS}.powcorrspctrm(1,:),'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(freq_ax,PowCorr.Enc{chanPair,iSS}.powcorrspctrm(1,:),'Color',Colors{2},'LineWidth',3);
semilogx(freq_ax,PowCorr.Maint{chanPair,iSS}.powcorrspctrm(1,:),'Color',Colors{3},'LineWidth',3);
xlim([4 30]);
set(gca,'FontSize',16,'box','off');
ylabel('Power Correlation');
xlabel('Frequency (Hz)');
SubjectPowCorr{selected_subj}.Fix = PowCorr.Fix{chanPair,iSS};
SubjectPowCorr{selected_subj}.Enc = PowCorr.Enc{chanPair,iSS};
SubjectPowCorr{selected_subj}.Maint = PowCorr.Maint{chanPair,iSS};


%% Time Frequency Granger
freq = [4:100];
[Granger_struct, Granger_TFR_ss] = ...
    get_TFR_Granger(freq,selected_subj,hipp_chanPairs,strChannelNameList,dataBipolar_SS);

%% Plot TFR granger for all set sizes
for iSS = 4:length(Granger_TFR_ss)
    
    grangerTimeAxis = [-6:0.25:2];
    grangerFreqAxis = [1:100];
    figure;
    clim = [-0.1 0.15];%[-0.05 0.05];
    contourf(grangerTimeAxis,grangerFreqAxis,-Granger_TFR_ss{iSS},100,'LineColor','none');
    %             contourf(grangerTimeAxis,grangerFreqAxis,temp,100,'LineColor','none');
    
    set(gca,'clim',clim,'yscale','log');
    set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
    colormap(bluewhitered(128));
    set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
    ylim([4 100]);
    ylim([5 30])
    colorbar;
    ylabel('Frequency (Hz)');
    xlabel('Time (s)');
    set(gca,'FontSize',16,'box','off');
    set(gcf,'Color','white')
end