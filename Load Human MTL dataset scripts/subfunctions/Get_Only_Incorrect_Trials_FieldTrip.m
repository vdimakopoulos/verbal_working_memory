function [ data_Incorrect, TrialInformationTable_Incorrect ] = Get_Only_Incorrect_Trials_FieldTrip( data_All, TrialInformationTable )

%% Correct trials
nTrialList_Incorrect = find(~(TrialInformationTable.Correct));
TrialInformationTable_Incorrect = TrialInformationTable(nTrialList_Incorrect,:);
% TrialInformationTable_Correct = TrialInformationTable(TrialInformationTable.Correct==1,:);
% nTrialList_Correct = TrialInformationTable_Correct.TrialNumber;

%% Select incorrect trials
cfg = [];
cfg.trials = nTrialList_Incorrect;

data_Incorrect = ft_selectdata(cfg,data_All);

%% Update trial information table
TrialInformationTable_Incorrect.TrialNumber(:) = (1:length(nTrialList_Incorrect));

end