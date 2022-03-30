function [ data_Correct, TrialInformationTable_Correct ] = Get_Only_Correct_Trials_FieldTrip( data_All, TrialInformationTable )

%% Correct trials
nTrialList_Correct = find(TrialInformationTable.Correct);
TrialInformationTable_Correct = TrialInformationTable(nTrialList_Correct,:);
% TrialInformationTable_Correct = TrialInformationTable(TrialInformationTable.Correct==1,:);
% nTrialList_Correct = TrialInformationTable_Correct.TrialNumber;

%% Select correct trials
cfg = [];
cfg.trials = nTrialList_Correct;

data_Correct = ft_selectdata(cfg,data_All);

%% Update trial information table
TrialInformationTable_Correct.TrialNumber(:) = (1:length(nTrialList_Correct));

end