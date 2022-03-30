function [macro_data, scalp_data, TrialInfoTable]= ...
    Load_Dataset_from_human_MTL(strDataPath_Macro,strDataPath_Scalp)

files_to_Load_Macro = dir([strDataPath_Macro,'*.mat']);
files_to_Load_Scalp = dir([strDataPath_Scalp,'*.mat']);

%% Check whether scalp EEG and iEEG data have the same number of sessions
assert(isequal(length(files_to_Load_Macro),length(files_to_Load_Scalp)),...
    'Number of sessions do not match between scalp EEG and iEEG')

%% Load each session file

%% Macro data
cd(strDataPath_Macro)
for iFile = 1:length(files_to_Load_Macro)
    data_ses(iFile) = load(files_to_Load_Macro(iFile).name); %load the files
end

%% Scalp data
cd(strDataPath_Scalp)
for iFile = 1:length(files_to_Load_Macro)
    data_scalp_ses(iFile) = load(files_to_Load_Scalp(iFile).name);%load the files
end

TrialInformationTable = [];
for iFile=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(iFile).tempMacro.TrialInformationTable];
end

%% Obtain Trial Information table
TrialInformationTable_Scalp = [];
for iFile=1:length(data_scalp_ses)
    TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_scalp_ses(iFile).tempScalp.TrialInformationTable];
end

TrialInfoTable.Macro = TrialInformationTable;
TrialInfoTable.Scalp = TrialInformationTable;

%% Merge_sessions
macro_data = data_ses(1).tempMacro.data;
for nSes = 2:length(data_ses)
    cfg = [];
    cfg.keepsampleinfo = 'no';
    macro_data = ft_appenddata(cfg,macro_data,data_ses(nSes).tempMacro.data);
    
end

scalp_data = data_scalp_ses(1).tempScalp.data
for nSes = 2:length(data_scalp_ses)
    cfg = [];
    scalp_data = ft_appenddata(cfg,scalp_data,data_scalp_ses(nSes).tempScalp.data);
    
end
end