function [ data_Correct, TrialInformationTable_Clean ] = Get_Only_Clean_Trials_FieldTrip( data_All, TrialInformationTable )

%% Correct trials
nTrialList_Clean = find(~TrialInformationTable.Artifact);
TrialInformationTable_Clean = TrialInformationTable(nTrialList_Clean,:);
% TrialInformationTable_Correct = TrialInformationTable(TrialInformationTable.Correct==1,:);
% nTrialList_Correct = TrialInformationTable_Correct.TrialNumber;

%% Select correct trials
cfg = [];
cfg.trials = nTrialList_Clean;
data_Correct = ft_selectdata(cfg,data_All);

%% Update trial information table
TrialInformationTable_Clean.Trialnumber(:) = (1:length(nTrialList_Clean));

end