function [ data_SessionRange, TrialInformationTable_SessionRange ] = Get_Selected_Session_Data_FieldTrip( data_All, nSessionRange, TrialInformationTable )

%% Correct trials
TrialInformationTable_SessionRange = TrialInformationTable(ismember(TrialInformationTable.Session,nSessionRange),:);
nTrialList_SessionRange = TrialInformationTable_SessionRange.TrialNumber;

%% Select correct trials
cfg = [];
cfg.trials = nTrialList_SessionRange;

data_SessionRange = ft_selectdata(cfg,data_All);

%% Update trial information table
TrialInformationTable_SessionRange.TrialNumber(:) = (1:length(nTrialList_SessionRange));

end