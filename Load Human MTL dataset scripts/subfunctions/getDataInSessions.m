function [TrialInfo, dataBipolar_Session, CorrectResponseRate_Session, dataBipolar_Ret_SS_session ] = ...
    getDataInSessions(dataBipolar_SS,TrialInformationTable_SS,TrialInfo_allTrials)

nSessionTrials = 50;
nTrialsHighWL = TrialInformationTable_SS{4}.TrialNumber
numberOfSessions = ceil(nTrialsHighWL./50);

nTrials_ss_6 = find(TrialInfo_allTrials.SetSize == 6);
nTrials_ss_8 = find(TrialInfo_allTrials.SetSize == 8);
nTrials_ss_68 = union(nTrials_ss_6,nTrials_ss_8);
n_Sessions_allTrials = ceil(nTrials_ss_68./50);

for iSes = 1:numberOfSessions(end)
    
    nTrials{iSes} = find(numberOfSessions == iSes);
    nSes_Trials_CorrIncorr{iSes} = find(n_Sessions_allTrials == iSes);
    
    cfg = [];
    cfg.trials = nTrials{iSes};
    dataBipolar_Session{iSes} = ft_selectdata(cfg,dataBipolar_SS{4});
    
    TrialInfo{iSes} = TrialInformationTable_SS{4}(nTrials{iSes},:);
    TrialSession_CorrectIncorrect = TrialInfo_allTrials(nSes_Trials_CorrIncorr{iSes},:);
    CorrectResponseRate_Session(iSes) = length(find(TrialInfo{iSes}.Correct==1))/length(TrialSession_CorrectIncorrect.Correct);
    
    
end

%% Select the latencies for every task period for each task period
[dataBipolar_Ret_SS_session{1}] = Extract_task_period_Data('fix',dataBipolar_Session); % fixation
[dataBipolar_Ret_SS_session{2}] = Extract_task_period_Data('encod',dataBipolar_Session); % encoding
[dataBipolar_Ret_SS_session{3}] = Extract_task_period_Data('maint',dataBipolar_Session); % maintenance



