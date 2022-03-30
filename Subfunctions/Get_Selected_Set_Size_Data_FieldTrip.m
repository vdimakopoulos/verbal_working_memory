function [ data_SingleSetSize, TrialInformationTable_SingleSetSize, nTrialList_TT ] = Get_Selected_Set_Size_Data_FieldTrip( data_All, nSetSize, TrialInformationTable )

%% Events to select set size
if(isequal(nSetSize,[6,8]))
    strSetSizeEvent = {['Stim_ss',num2str(nSetSize(1))],['Stim_ss',num2str(nSetSize(2))]};
elseif(ismember(nSetSize,[4,6,8]))
    strSetSizeEvent = {['Stim_ss',num2str(nSetSize)]};
end

%% Trial information table for the set size
TrialInformationTable_SingleSetSize = TrialInformationTable(ismember(TrialInformationTable.SetSize,nSetSize),:);

%% Trial numbers for the set size
% Correct Trials
nTrialList_TT = TrialInformationTable_SingleSetSize.TrialNumber;

%Incorrect trials
% nTrialList_TT = find(TrialInformationTable.SetSize == nSetSize)%TrialInformationTable_SingleSetSize.TrialNumber;
% if (isequal(nSetSize,[6,8]))
%     nTrialList_TT = find(TrialInformationTable.SetSize(ismember(TrialInformationTable.SetSize,nSetSize)))
% elseif(ismember(nSetSize,[4,6,8]))
%        nTrialList_TT = find(TrialInformationTable.SetSize(ismember(TrialInformationTable.SetSize,nSetSize)))
% 
% end
%% Select trials
cfg = [];
cfg.trials = nTrialList_TT;

data_SingleSetSize = ft_selectdata(cfg,data_All);

%% Check number of trials - not implemented
% nNumberOfTrials_EEG = EEG_SingleSetSize.trials;
% nNumberOfTrials_TT = size(TrialInformationTable_SingleSetSize,1);
% 
% if(nNumberOfTrials_EEG~=nNumberOfTrials_TT)
%     error('Different number of trials from event selection and trial information table')
% end

%% Check list of trial numbers - we can't now because the trial numbers in the EEGLAB structure are not updated
% strEvents = {EEG_SingleSetSize.event.type}';
% strEvents = strEvents(~cellfun(@isempty,strfind(strEvents,'Fix_')));
% strEvents = strrep(strEvents,'Fix_s','');
% strEvents = cellfun(@(x) x(4:end),strEvents,'UniformOutput',0);
% % ind = strfind(strEvents,'_ss');
% strEvents = cellfun(@(x,y) x(1:y-1),strEvents,ind,'UniformOutput',0);
% nTrialList_EEG = str2double(strEvents);

% nTrialList_TT = TrialInformationTable_SingleSetSize.TrialNumber;

end