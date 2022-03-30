function [dataBipolarScalpResampled,dataBipolar_Ret_SS,TrialInformationTable] = LoadHumanMTL_EEG(pID,analysis_type)

%% Parameters
% analysis_type = 'incorrect trials'; %'incorrect trials'; %or 'correct trials'
Extra_Analysis = 0;
%% Load the data
strDataPath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Dataset_2020_MTL_iEEG_scalpEEG\';
strDataPath_Macro = [strDataPath,'Macro\'];
strDataPath_Scalp = [strDataPath,'Scalp\'];
nSubjects = [1:9];
strPat_IDs = {'WC','SW','IA','HJ','PL','AM','HF','MA','MS'};
nSubject_ID = [28 22 19 30 33 13 23 29 16];
% specify the patient's data directory
selected_subj =  pID;
strPatientDataPath_Macro = sprintf('%s%d %s\\',strDataPath_Macro,nSubject_ID(selected_subj),strPat_IDs{selected_subj});
strPatientDataPath_Scalp = sprintf('%s%d %s\\',strDataPath_Scalp,nSubject_ID(selected_subj),strPat_IDs{selected_subj});

[dataBipolar, dataBipolar_Scalp, TrialInformationTable] = ...
    Load_Dataset_from_human_MTL(strPatientDataPath_Macro,strPatientDataPath_Scalp);

% cd(strPaths.Main)

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
cfg.reref         =  'yes';
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
[montage,chans_to_be_bipolar,dataBip,Bip_chans] = GetDataInBipolar(dataBipolar_reref,selected_subj,0);

%% Downsample to 200 Hz
cfg = [];
cfg.resamplefs = 80;%200;
dataBipolarResampled = ft_resampledata(cfg,dataBip);
dataBipolarScalpResampled = ft_resampledata(cfg,dataBipolarScalp_reref);

%% Append scalp and iEEG data
cfg = [];
dataBipolarResampled.time = dataBipolarScalpResampled.time;
data_bipolar_appended = ft_appenddata(cfg,dataBipolarResampled,dataBipolarScalpResampled);

%% macro data
% Select only correct trials
switch analysis_type
    case 'correct trials'
        [dataBipolar_SS,TrialInformationTableCorrect] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalpResampled,TrialInformationTable_iEEG_clean);
    case 'incorrect trials'
        [dataBipolar_SS, TrialInformationTableCorrect] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolarScalpResampled, TrialInformationTable_iEEG_clean )
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
