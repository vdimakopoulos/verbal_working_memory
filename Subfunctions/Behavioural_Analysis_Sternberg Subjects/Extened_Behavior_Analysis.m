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
%%
strPlotColors = {'b','g','r'};
%% List of patients and sessions
% PatientSessionList = ...
%     [
% %     4,1
% %     4,3
% %     4,4
% %     10,1
% %     10,2
% %     10,3
% %     10,4
% %     10,5
%     13,1
%     13,2
%     13,3
%     13,5
%     13,6
%     13,7
%     13,8
%     16,1
%     16,2
%     19,1
%     19,2
%     19,3
%     22,1
%     22,2
%     22,3
%     22,4
%     22,5
%     22,6
%     22,7
%     23,1
%     23,2
%     23,3
%     23,4
%     28,1
%     28,2
%     28,3
%     28,4
%     29,1
%     29,2
%     29,3
%     29,4
% %     29,5
%     30,1
%     30,2
% %     32,1
% %     32,2
% %     32,3
%     33,1
%     33,2
%     33,3
%     ];

PatientSessionList = ...
    [
    10,1
    10,2
    10,3
    10,4
    10,5
    12,1
    12,2
    14,1
    14,2
    20,1
    20,2
    26,1
    26,2
    37,1
    37,2
    37,3
    37,4
    37,5
    37,6
    38,1
    38,2
    38,3
    38,4
    38,5
    38,6
    38,7
    40,1
    40,2
    40,3
    40,4
    40,5
    40,6
    40,7
    40,8
    42,1
    42,2
    44,1
    44,2
    44,3
    44,4
    44,5
    45,1
    45,2
    45,3
    45,4
    45,5
    45,6
    
    ];
nPart = 1;
nNumberOfPatientsSubjects = size(PatientSessionList,1);
nSubjects  = unique(PatientSessionList(:,1));
%% For each patient/session on list, load data and extract rates
% Initialization for maximum number of trials among each set size / response accuracy pair
nSetSizeValues = [4,6,8];
nCorrIncorrValues = [1,0];
TrialInformationTable_AllSessions = [];

%%
for nPatientSubject = 1:nNumberOfPatientsSubjects
    nPatient =  PatientSessionList(nPatientSubject,1);
    nSession =  PatientSessionList(nPatientSubject,2);
    pID = find(nPatient==nSubjects);
    TrialInformationTable = TrialInformationTable_AllSessions{pID};
    %% Load data
%     strFieldName = 'Macro';
% %     try
%         TrialInformationTable = load([strPaths.NLXDataExtraction.Sessions.(strFieldName),strFieldName,'_Data_Sessions_Patient_',num2str(nPatient,'%.2d'),'_Session_',num2str(nSession,'%.2d'),'_Part_',num2str(nPart,'%.2d')],...
%             'TrialInformationTable');
%         TrialInformationTable = TrialInformationTable.TrialInformationTable;
%     catch
%         strListOfFiles = dir([strPaths.NLXDataExtraction.Sessions.(strFieldName),strFieldName,'_Data_Sessions_Patient_',num2str(nPatient,'%.2d'),'_Session_',num2str(nSession,'%.2d'),'_Part_',num2str(nPart,'%.2d'),'*']);
%         TrialInformationTable = load([strPaths.NLXDataExtraction.Sessions.(strFieldName),strListOfFiles(1).name],'TrialInformationTable');
%         TrialInformationTable = TrialInformationTable.TrialInformationTable;
%     end
    
%% Collect trial information table for all patients 
% TrialInformationTable.Patient = nPatient*ones(size(TrialInformationTable,1),1);
% TrialInformationTable.Session = nSession*ones(size(TrialInformationTable,1),1);
% 
% TrialInformationTable_AllSessions = [TrialInformationTable_AllSessions;TrialInformationTable];

    %% Accuracy
    Accuracy_PatientSubject(nPatientSubject) = length(find(TrialInformationTable.Correct==1))/length(TrialInformationTable.Correct);
    Accuracy_Patient_Session(nPatient,nSession) = length(find(TrialInformationTable.Correct==1))/length(TrialInformationTable.Correct);
    
    %% Accuracy for each set size
    for iSS = 1:3
        TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.SetSize==nSetSizeValues(iSS),:);
        Accuracy_PatientSubject_SS(nPatientSubject,iSS) = length(find(TrialInformationTable_Temp.Correct==1))/length(TrialInformationTable_Temp.Correct);
        Accuracy_Patient_Session_SS(nPatient,nSession,iSS) = length(find(TrialInformationTable_Temp.Correct==1))/length(TrialInformationTable_Temp.Correct);
    end
    clear TrialInformationTable_Temp
    
    %% Accuracy for match/mismatch
    for iMM = 1:2
        TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Match==iMM,:);
        Accuracy_PatientSubject_MM(nPatientSubject,iMM) = length(find(TrialInformationTable_Temp.Correct==1))/length(TrialInformationTable_Temp.Correct);
        Accuracy_Patient_Session_MM(nPatient,nSession,iMM) = length(find(TrialInformationTable_Temp.Correct==1))/length(TrialInformationTable_Temp.Correct);
    end
    clear TrialInformationTable_Temp
    
    %% Reaction time for each set size
    for iSS = 1:3
        ResponseTime_PatientSubject_SS{nPatientSubject}{iSS} = TrialInformationTable(TrialInformationTable.SetSize==nSetSizeValues(iSS),:).ResponseTime;
        ResponseTime_Patient_Session_SS{nPatient}{nSession}{iSS} = TrialInformationTable(TrialInformationTable.SetSize==nSetSizeValues(iSS),:).ResponseTime;
        
        ResponseTime_Mean_PatientSubject_SS(nPatientSubject,iSS) = mean(ResponseTime_PatientSubject_SS{nPatientSubject}{iSS});
        ResponseTime_Median_PatientSubject_SS(nPatientSubject,iSS) = median(ResponseTime_PatientSubject_SS{nPatientSubject}{iSS});
    end
    
    %% Reaction time for correct/incorrect answers
    for iDec = 1:2
        TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Correct==nCorrIncorrValues(iDec),:);
        ResponseTime_PatientSubject_CI{nPatientSubject}{iDec} = TrialInformationTable_Temp.ResponseTime;
        ResponseTime_Patient_Session_CI{nPatient}{nSession}{iDec} = TrialInformationTable_Temp.ResponseTime;
        
        ResponseTime_Mean_PatientSubject_CI(nPatientSubject,iDec) = mean(ResponseTime_PatientSubject_CI{nPatientSubject}{iDec});
        ResponseTime_Median_PatientSubject_CI(nPatientSubject,iDec) = median(ResponseTime_PatientSubject_CI{nPatientSubject}{iDec});
    end
    clear TrialInformationTable_Temp

    %% Reaction time for correct/incorrect answers for each set size
    for iDec = 1:2
        for iSS = 1:3
            TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Correct==nCorrIncorrValues(iDec),:);
            TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.SetSize==nSetSizeValues(iSS),:);
            ResponseTime_PatientSubject_CI_SS{nPatientSubject}{iDec}{iSS} = TrialInformationTable_Temp.ResponseTime;
            ResponseTime_Patient_Session_CI_SS{nPatient}{nSession}{iDec}{iSS} = TrialInformationTable_Temp.ResponseTime;
            
            ResponseTime_Mean_PatientSubject_CI_SS(nPatientSubject,iDec,iSS) = mean(ResponseTime_PatientSubject_CI_SS{nPatientSubject}{iDec}{iSS});
            ResponseTime_Median_PatientSubject_CI_SS(nPatientSubject,iDec,iSS) = median(ResponseTime_PatientSubject_CI_SS{nPatientSubject}{iDec}{iSS});
        end
    end
    clear TrialInformationTable_Temp

    %% Reaction time for correct match/mismatch answers
    for iMM = 1:2
        TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Match==iMM,:);
        TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.Correct==1,:);
        ResponseTime_PatientSubjects_MM{nPatientSubject}{iMM} = TrialInformationTable_Temp.ResponseTime;
        ResponseTime_Patient_Session_MM{nPatient}{nSession}{iMM} = TrialInformationTable_Temp.ResponseTime;
        
        ResponseTime_Mean_PatientSubject_MM(nPatientSubject,iMM) = mean(ResponseTime_PatientSubjects_MM{nPatientSubject}{iMM});
        ResponseTime_Median_PatientSubject_MM(nPatientSubject,iMM) = median(ResponseTime_PatientSubjects_MM{nPatientSubject}{iMM});
    end
    clear TrialInformationTable_Temp

    %% Reaction time for correct match/mismatch answers for each set size
    for iMM = 1:2
        for iSS = 1:3
        TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Match==iMM,:);
        TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.Correct==1,:);
        TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.SetSize==nSetSizeValues(iSS),:);
        ResponseTime_PatientSubjects_MM_SS{nPatientSubject}{iMM}{iSS} = TrialInformationTable_Temp.ResponseTime;
        ResponseTime_Patient_Session_MM_SS{nPatient}{nSession}{iMM}{iSS} = TrialInformationTable_Temp.ResponseTime;
        
        ResponseTime_Mean_PatientSubject_MM_SS(nPatientSubject,iMM,iSS) = mean(ResponseTime_PatientSubjects_MM_SS{nPatientSubject}{iMM}{iSS});
        ResponseTime_Median_PatientSubject_MM_SS(nPatientSubject,iMM,iSS) = median(ResponseTime_PatientSubjects_MM_SS{nPatientSubject}{iMM}{iSS});
        end
    end
    clear TrialInformationTable_Temp
    
    %% Reaction time for correct/incorrect match/mismatch trials
    for iDec = 1:2
        for iMM = 1:2
            TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Correct==nCorrIncorrValues(iDec),:);
            TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.Match==iMM,:);
            
            ResponseTime_PatientSubject_CI_MM{nPatientSubject}{iDec}{iMM} = TrialInformationTable_Temp.ResponseTime;
            ResponseTime_Patient_Session_CI_MM{nPatient}{nSession}{iDec}{iMM} = TrialInformationTable_Temp.ResponseTime;
            
            ResponseTime_Mean_PatientSubject_CI_MM(nPatientSubject,iDec,iMM) = mean(ResponseTime_PatientSubject_CI_MM{nPatientSubject}{iDec}{iMM});
            ResponseTime_Median_PatientSubject_CI_MM(nPatientSubject,iDec,iMM) = median(ResponseTime_PatientSubject_CI_MM{nPatientSubject}{iDec}{iMM});
        end
    end
    clear TrialInformationTable_Temp    
    
    %% Reaction time for correct/incorrect match/mismatch trials for each set size
    for iDec = 1:2
        for iMM = 1:2
            for iSS = 1:3
                TrialInformationTable_Temp = TrialInformationTable(TrialInformationTable.Correct==nCorrIncorrValues(iDec),:);
                TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.Match==iMM,:);
                TrialInformationTable_Temp = TrialInformationTable_Temp(TrialInformationTable_Temp.SetSize==nSetSizeValues(iSS),:);
                
                ResponseTime_PatientSubject_CI_MM_SS{nPatientSubject}{iDec}{iMM}{iSS} = TrialInformationTable_Temp.ResponseTime;
                ResponseTime_Patient_Session_CI_MM_SS{nPatient}{nSession}{iDec}{iMM}{iSS} = TrialInformationTable_Temp.ResponseTime;
                
                ResponseTime_Mean_PatientSubject_CI_MM_SS(nPatientSubject,iDec,iMM,iSS) = mean(ResponseTime_PatientSubject_CI_MM_SS{nPatientSubject}{iDec}{iMM}{iSS});
                ResponseTime_Median_PatientSubject_CI_MM_SS(nPatientSubject,iDec,iMM,iSS) = median(ResponseTime_PatientSubject_CI_MM_SS{nPatientSubject}{iDec}{iMM}{iSS});
            end
        end
    end
    clear TrialInformationTable_Temp        
    
    
end

%% Write results to table
% indNew = [size(TrialInformationTable_AllSessions,2),size(TrialInformationTable_AllSessions,2)-1,1:size(TrialInformationTable_AllSessions,2)-2];
% TrialInformationTable_AllSessions = TrialInformationTable_AllSessions(:,indNew);
% writetable(TrialInformationTable_AllSessions,strRecordingInformationFilePath,'Sheet','Behavioral Analysis')

%% Accuracy for all sessions
Accuracy_PatientSubject_Ord = sort(Accuracy_PatientSubject)*100;
meanAccuracy_PatientSubject_Ord = mean(Accuracy_PatientSubject_Ord);
semAccuracy_PatientSubject_Ord = std(Accuracy_PatientSubject_Ord)/sqrt(length(Accuracy_PatientSubject_Ord));

figure(101);clf
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
plot(Accuracy_PatientSubject_Ord,'k.','MarkerSize',10)
hold on
tempAxis = 1:length(Accuracy_PatientSubject_Ord);
plot(tempAxis,meanAccuracy_PatientSubject_Ord*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
plot(tempAxis,meanAccuracy_PatientSubject_Ord*ones(1,length(tempAxis))+semAccuracy_PatientSubject_Ord,'Color','b')
plot(tempAxis,meanAccuracy_PatientSubject_Ord*ones(1,length(tempAxis))-semAccuracy_PatientSubject_Ord,'Color','b')
% legend('Accuracy','Mean','Mean+s.e.m','Mean-s.e.m')
ylabel('Accuracy (%)')
xlabel('Session no.')
box off
axis square
ylim([75,100])
clear tempAxis

%% Accuracy for set sizes / PUBLICATION
clear indAccuracyOrder
for iSS = 1:3
    temp = Accuracy_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(Accuracy_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
figure(102);clf
clear p1 p2
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(Accuracy_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),Accuracy_PatientSubject_SS(iPatientSubject,:)*100,'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),Accuracy_PatientSubject_SS(iPatientSubject,iSS)*100,[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(Accuracy_PatientSubject_SS(:,iSS)*100);
    semAccTemp = std(Accuracy_PatientSubject_SS(:,iSS)*100)/sqrt(length(Accuracy_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(Accuracy_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Accuracy (%)')
xlabel('Set size')
box off
axis square
xlim([3,9])
ylim([50,100])
set(gca,'XTick',[4,6,8])
set(102,'Position',[680   636   306   342])

%% Accuracy for IN/OUT / PUBLICATION
figure
boxplot(Accuracy_PatientSubject_MM*100)
set(gca,'XTickLabel',{'IN','OUT'})
ylabel('Accuracy (%)')
title('Accuracy for IN/OUT / Sessions')

[h,p] = ttest(Accuracy_PatientSubject_MM(:,1),Accuracy_PatientSubject_MM(:,2))
a = {Accuracy_PatientSubject_MM(:,1)',Accuracy_PatientSubject_MM(:,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000);
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)

%% Response time for for IN/OUT
[h,p] = ttest(Accuracy_PatientSubject_MM(:,1),Accuracy_PatientSubject_MM(:,2))
a = {ResponseTime_Median_PatientSubject_MM(:,1)',ResponseTime_Median_PatientSubject_MM(:,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',20000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)

%% Response times for IN/OUT only correct
% var_PatientSubject_MM = squeeze(ResponseTime_Mean_PatientSubject_CI_MM(:,1,:));
a = {ResponseTime_Mean_PatientSubject_CI_MM(:,1,1)',ResponseTime_Mean_PatientSubject_CI_MM(:,1,2)'};
[h,p] = ttest(ResponseTime_Mean_PatientSubject_CI_MM(:,1,1),ResponseTime_Mean_PatientSubject_CI_MM(:,1,2))
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)

%% Response time for Correct/Incorrect
temp = ResponseTime_Mean_PatientSubject_CI; % ResponseTime_Median_PatientSubject_CI
temp(find(sum(isnan(ResponseTime_Mean_PatientSubject_CI),2)),:) = []; % ResponseTime_Median_PatientSubject_CI
[h,p] = ttest(temp(:,1),temp(:,2));        
fprintf('\nResponse time - Correct %.2f %.2f s - Incorrect %.2f %.2f s - t-test; p = %.6f\n',mean(temp(:,1)),std(temp(:,1)),mean(temp(:,2)),std(temp(:,2)),p)
a = {temp(:,1)',temp(:,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.5f\n',stats,pvals)

%% Accuracy for setsize / PUBLICATION
[p,tbl,stats] = anova1(Accuracy_PatientSubject_SS)
a = {Accuracy_PatientSubject_SS(:,1)',Accuracy_PatientSubject_SS(:,2)',Accuracy_PatientSubject_SS(:,3)'}; % pseudo 'unpaired'
[F df pvals] = statcond(a,'method','perm','naccu',2000,'paired','on') % perform ANOVA
fprintf('\nF %.4f .4f; p = %.4f\n',F,pvals)

%% Response time for correct setsize
var_PatientSubject_SS = squeeze(ResponseTime_Mean_PatientSubject_CI_SS(:,1,:));
[p,tbl,stats] = anova1(var_PatientSubject_SS)
a = {ResponseTime_Mean_PatientSubject_CI_SS(:,1,1)',ResponseTime_Mean_PatientSubject_CI_SS(:,1,2)',ResponseTime_Mean_PatientSubject_CI_SS(:,1,3)'}; % pseudo 'unpaired'
[F df pvals] = statcond(a,'method','perm','naccu',2000,'paired','on') % perform an unpaired ANOVA
    
%% Response times / match mismatch / only correct / setsize
a = {ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,1)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,1)'};
% a = {ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,2)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,2)'};
% a = {ResponseTime_Median_PatientSubject_CI_MM_SS(:,1,1,3)',ResponseTime_Median_PatientSubject_CI_MM_SS(:,1,2,3)'};
% a = {[ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,1)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,2)'],...
%     [ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,1)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,2)']};
% a = {[ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,2)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,1,3)'],...
%     [ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,2)',ResponseTime_Mean_PatientSubject_CI_MM_SS(:,1,2,3)']};
a{1}(isnan(a{2})) = [];
a{2}(isnan(a{2})) = [];
a{2}(isnan(a{1})) = [];
a{1}(isnan(a{1})) = [];
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)
 
%%
%%
%%
%% Statistics at trial level
%% Accuracy for setsizes
Accuracy_TrialLevel_SS = [];
for iSS = 1:length(nSetSizeValues)
Table_SS = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.SetSize==nSetSizeValues(iSS),:);
nAll = size(Table_SS,1);
nCorr = length(find(Table_SS.Correct==1));
Accuracy_TrialLevel_SS(iSS) = nCorr/nAll;
end
Accuracy_TrialLevel_SS

%% Accuracy for IN/OUT
Accuracy_TrialLevel_MM = [];
for iMM = 1:2
Table_MM = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.Match==iMM,:);
nAll = size(Table_MM,1);
nCorr = length(find(Table_MM.Correct==1));
Accuracy_TrialLevel_MM(iMM) = nCorr/nAll;
end
Accuracy_TrialLevel_MM

%% Reaction time for correct set sizes
ReactionTime_TrialLevel_Corr_SS = {};
varSS = [];
varRT = [];
for iSS = 1:length(nSetSizeValues)
Table_SS = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.SetSize==nSetSizeValues(iSS),:);
Table_SS = Table_SS(Table_SS.Correct==1,:);
varSS = [varSS;nSetSizeValues(iSS)*ones(size(Table_SS,1),1)];
varRT = [varRT;Table_SS.ResponseTime];
ReactionTime_TrialLevel_Corr_SS{iSS} = Table_SS.ResponseTime;
end
ReactionTime_TrialLevel_Corr_SS

% varSS = varSS/2-1;
% var_PatientSubject_SS = squeeze(ResponseTime_Mean_PatientSubject_CI_SS(:,1,:));
% [p,tbl,stats] = anova1(var_PatientSubject_SS)
a = {ReactionTime_TrialLevel_Corr_SS{1}',ReactionTime_TrialLevel_Corr_SS{2}',ReactionTime_TrialLevel_Corr_SS{3}'}; % pseudo 'unpaired'
[F df pvals] = statcond(a,'method','perm','naccu',2000,'paired','off') % perform an unpaired ANOVA
    
[mean(ReactionTime_TrialLevel_Corr_SS{1}),mean(ReactionTime_TrialLevel_Corr_SS{2}),mean(ReactionTime_TrialLevel_Corr_SS{3})]
LM = fitlm(varSS,varRT,'linear');

yaxis_temp = 0:5;
figure
boxplot(varRT',varSS')
hold on
plot(yaxis_temp,LM.Coefficients.Estimate(2)*yaxis_temp+LM.Coefficients.Estimate(1),'r')
LM.Coefficients.Estimate(2)

% sprintf('The mean response time (RT) for correct trials increased with workload (Figure TASK c, 48 ms/item, F2,2007 = 21.10, P = 0.0002, permuted-repeated measures ANOVA)',...
%     . 

%% Reaction time for correct IN/OUT
ReactionTime_TrialLevel_Corr_MM = {};
for iMM = 1:2
Table_MM = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.Match==iMM,:);
Table_MM = Table_MM(Table_MM.Correct==1,:);
ReactionTime_TrialLevel_Corr_MM{iMM} = Table_MM.ResponseTime;
end
ReactionTime_TrialLevel_Corr_MM

[h,p] = ttest2(ReactionTime_TrialLevel_Corr_MM{1},ReactionTime_TrialLevel_Corr_MM{2})
fprintf('\nttest2 t-value %.2f; p = %.4f\n',0,p)
a = {ReactionTime_TrialLevel_Corr_MM{1}',ReactionTime_TrialLevel_Corr_MM{2}'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','off','verbose','off');
fprintf('\npermutation t-value %.2f; p = %.4f\n',stats,pvals)

%% Reaction time for correct IN/OUT / each session 
ReactionTime_TrialLevel_Corr_MM = {};
for iMM = 1:2
Table_MM = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.Match==iMM,:);
Table_MM = Table_MM(Table_MM.Correct==1,:);
ReactionTime_TrialLevel_Corr_MM{iMM} = Table_MM.ResponseTime;
end
ReactionTime_TrialLevel_Corr_MM

[h,p] = ttest2(ReactionTime_TrialLevel_Corr_MM{1},ReactionTime_TrialLevel_Corr_MM{2})
fprintf('\nttest2 t-value %.2f; p = %.4f\n',stats,p)
a = {ReactionTime_TrialLevel_Corr_MM{1}',ReactionTime_TrialLevel_Corr_MM{2}'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','off','verbose','off');
fprintf('\npermutation t-value %.2f; p = %.4f\n',stats,pvals)

%% Reaction time for correct incorrect / PUBLICATION
ReactionTime_TrialLevel_Corr_MM = {};
for iC = 1:2
Table_CI = TrialInformationTable_AllSessions(TrialInformationTable_AllSessions.Correct==nCorrIncorrValues(iC),:);
ReactionTime_TrialLevel_CI{iC} = Table_CI.ResponseTime;
end
ReactionTime_TrialLevel_CI

[h,p] = ttest2(ReactionTime_TrialLevel_CI{1},ReactionTime_TrialLevel_CI{2})
fprintf('\nttest2 t-value %.2f; p = %.4f\n',stats,p)
a = {ReactionTime_TrialLevel_CI{1}',ReactionTime_TrialLevel_CI{2}'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',20000,'paired','off','verbose','off');
fprintf('\npermutation t-value %.2f; p = %.4f\n',stats,pvals)

cellfun(@mean,ReactionTime_TrialLevel_CI)
cellfun(@std,ReactionTime_TrialLevel_CI)

% Phrase for publication
sprintf('Correct IN/OUT decisions were made more quickly than incorrect decisions \n%.2f ± %.2f s versus %.2f ± %.2f s, paired permutation t-test; t-value %.2f; P = %.1e ).',...
    mean(ReactionTime_TrialLevel_CI{1}),std(ReactionTime_TrialLevel_CI{1}),mean(ReactionTime_TrialLevel_CI{2}),std(ReactionTime_TrialLevel_CI{2}),...
    -stats,pvals)

%% Reaction time for set sizes / PUBLICATION
var_PatientSubject_SS = ResponseTime_Median_PatientSubject_SS;
clear indAccuracyOrder
for iSS = 1:3
    temp = var_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(var_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
figure(103);clf
clear p1 p2
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(var_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),var_PatientSubject_SS(iPatientSubject,:),'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),var_PatientSubject_SS(iPatientSubject,iSS),[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(var_PatientSubject_SS(:,iSS));
    semAccTemp = std(var_PatientSubject_SS(:,iSS))/sqrt(length(var_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(var_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Reaction time (s)')
xlabel('Session no.')
box off
axis square
xlim([3,9])
% ylim([50,100])
set(gca,'XTick',[4,6,8])

%% Reaction time for set sizes / PUBLICATION / Task b c
fig = figure;
subplot(1,2,1)
clear indAccuracyOrder
for iSS = 1:3
    temp =  Accuracy_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.01;
xoff = -floor(size(Accuracy_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
% figure(102);clf
clear p1 p2
% set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])


for iPatientSubject = 1:size(Accuracy_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),Accuracy_PatientSubject_SS(iPatientSubject,:)*100,'k:');
    hold on
% 
% 
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),Accuracy_PatientSubject_SS(iPatientSubject,iSS)*100,[strPlotColors{iSS},'.']);
    end
end


for iSS = 1:3
    meanAccTemp = mean(Accuracy_PatientSubject_SS(:,iSS)*100);
    semAccTemp = std(Accuracy_PatientSubject_SS(:,iSS)*100)/sqrt(length(Accuracy_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)%+xoff+(0:length(Accuracy_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Accuracy (%)')
xlabel('Set size')
box off
axis square
xlim([3,9])
ylim([60,100])
set(gca,'XTick',[4,6,8])
set(gca,'YTick',60:10:100)
% set(102,'Position',[680   636   306   342])
ttLetter(1) = text(0,0,'b','Units','normalized');


subplot(1,2,2)
var_PatientSubject_SS = ResponseTime_Median_PatientSubject_SS;
clear indAccuracyOrder
for iSS = 1:3
    temp = var_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(var_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
% figure(103);clf
clear p1 p2
% set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(var_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),var_PatientSubject_SS(iPatientSubject,:),'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),var_PatientSubject_SS(iPatientSubject,iSS),[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(var_PatientSubject_SS(:,iSS));
    semAccTemp = std(var_PatientSubject_SS(:,iSS))/sqrt(length(var_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(var_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Reaction time (s)')
xlabel('Set size')
box off
axis square
xlim([3,9])
% ylim([50,100])
set(gca,'XTick',[4,6,8])
set(gca,'YTick',0.8:0.4:2.4)

ttLetter(2) = text(0,0,'c','Units','normalized');

for ii = 1:length(ttLetter)
    set(ttLetter(ii),'HorizontalAlignment','Left')
    set(ttLetter(ii),'FontSize',12)
    set(ttLetter(ii),'FontWeight','bold')
    set(ttLetter(ii),'Position',[-0.1,1.1])
end


figure;
subplot(1,2,1)
plotSpread(Accuracy_PatientSubject_SS*100,'xNames',{'4','6','8'}, 'distributionMarkers', {'.', '.', '.'}, 'distributionColors',{'k','k','k'},'binWidth',0.001)
boxplot(Accuracy_PatientSubject_SS*100);

set(gca,'FontSize',8,'box','off')
xlabel('Set Size')
ylabel('Accuracy (%)')
ylim([60 100]);set(gca,'YTick',[60:10:100],'YTickLabel',[60:10:100],'XTick',[1 2 3],'XTickLabel',[4 6 8]);
axis square


subplot(1,2,2)
plotSpread(var_PatientSubject_SS,'xNames',{'4','6','8'}, 'distributionMarkers', {'.', '.', '.'}, 'distributionColors',{'k','k','k'},'binWidth',0.001)
boxplot(var_PatientSubject_SS);
set(gca,'FontSize',8,'box','off')
ylim([0.8 2.4]);set(gca,'YTick',[0.8:0.4:2.4],'YTickLabel',[0.8:0.4:2.4])
xlabel('Set Size')
ylabel('Reaction Time (s)')
axis square

%% Reaction time for set sizes / Only correct trials
var_PatientSubject_SS = squeeze(ResponseTime_Median_PatientSubject_CI_SS(:,1,:));
clear indAccuracyOrder
for iSS = 1:3
    temp = var_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(var_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
figure(104);clf
clear p1 p2
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(var_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),var_PatientSubject_SS(iPatientSubject,:),'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),var_PatientSubject_SS(iPatientSubject,iSS),[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(var_PatientSubject_SS(:,iSS));
    semAccTemp = std(var_PatientSubject_SS(:,iSS))/sqrt(length(var_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(var_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Reaction time (s)')
xlabel('Set size')
box off
axis square
xlim([3,9])
ylim([0.5,2.5])
set(gca,'XTick',[4,6,8])
set(104,'Position',[680   636   306   342])

%% Reaction time forset sizes only correct trials
var_PatientSubject_SS = squeeze(ResponseTime_Median_PatientSubject_CI_SS(:,1,:));
clear indAccuracyOrder
for iSS = 1:3
    temp = var_PatientSubject_SS(:,3);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

% Plot
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(var_PatientSubject_SS,1)/2)*c1;
c2 = 0.01;
figure(104);clf
clear p1 p2
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(var_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),var_PatientSubject_SS(iPatientSubject,:),'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),var_PatientSubject_SS(iPatientSubject,iSS),[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(var_PatientSubject_SS(:,iSS));
    semAccTemp = std(var_PatientSubject_SS(:,iSS))/sqrt(length(var_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(var_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Reaction time (s)')
xlabel('Session no.')
box off
axis square
xlim([3,9])
% ylim([50,100])

%% Reaction time for set sizes
clear indAccuracyOrder
for iSS = 1:3
    temp = Accuracy_PatientSubject_SS(:,iSS);
    [~,temp2] = sort(temp,'ascend');
    [~,temp3] = sort(temp2,'ascend');
    indAccuracyOrder(:,iSS) = temp3;
end
clear temp temp2 temp3

%% Plot accuracy time for set sizes
strPlotColors = {'b','g','r'};
c1 = 0.02;
xoff = -floor(size(Accuracy_PatientSubject_SS,1)/2)*c1;
c2 = 0.02;
figure(102);clf
clear p1 p2
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for iPatientSubject = 1:size(Accuracy_PatientSubject_SS,1)
    p1 = plot(nSetSizeValues+xoff+c2*indAccuracyOrder(iPatientSubject,:),Accuracy_PatientSubject_SS(iPatientSubject,:),'k:');
    hold on
    for iSS = 1:3
    p2(iSS) = plot(nSetSizeValues(iSS)+xoff+c2*indAccuracyOrder(iPatientSubject,iSS),Accuracy_PatientSubject_SS(iPatientSubject,iSS),[strPlotColors{iSS},'.']);
    end
end

for iSS = 1:3
    meanAccTemp = mean(Accuracy_PatientSubject_SS(:,iSS));
    semAccTemp = std(Accuracy_PatientSubject_SS(:,iSS))/sqrt(length(Accuracy_PatientSubject_SS(:,iSS)));
    tempAxis = nSetSizeValues(iSS)+xoff+(0:length(Accuracy_PatientSubject_SS(:,iSS))-1)*c2;
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis)),'Color','b','LineWidth',1.5)
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))+semAccTemp,'Color','b')
    plot(tempAxis,meanAccTemp*ones(1,length(tempAxis))-semAccTemp,'Color','b')
end
% for iSS = 1:3
% plot(nSetSizeValues(iSS),Accuracy_PatientSubject_SS(:,iSS)*100,[strPlotColors{iSS},'.'],'MarkerSize',10)
% hold on
% end
ylabel('Accuracy (%)')
xlabel('Session no.')
box off
axis square
% xlim([3.5,8.5])
% ylim([60,100])


%% Plot results
figure(102);clf
set(gca,'Position',[0.1586    0.1100    0.7464    0.8150])
for nPatientSubject = 1:nNumberOfPatientsSubjects
    for iSS = 1:3
        ResponseTime_SingleSessionSS = ResponseTime_PatientSubject_SS{nPatientSubject}{iSS};
        meanResponseTime_SingleSession_SS(iSS) = mean(ResponseTime_SingleSessionSS);
    end
    meanResponseTime_SS(nPatientSubject,:) = meanResponseTime_SingleSession_SS;
    plot(nSetSizeValues,meanResponseTime_SingleSession_SS,'k.','MarkerSize',8)
    plot(nSetSizeValues,meanResponseTime_SingleSession_SS,'k:')
    hold on
    xlim([3.5,8.5])
    ylabel('Reaction time (s)')
    xlabel('Set size')
    set(gca,'XTick',nSetSizeValues)
end
box off
axis square

%
for iSS = 1:3
    meanResponseTime_PatientSubject_Ord = mean(meanResponseTime_SS(:,iSS));
    semResponseTime_PatientSubject_Ord = std(meanResponseTime_SS(:,iSS))/sqrt(size(meanResponseTime_SS(:,iSS),1));
    
    ssAxis = [nSetSizeValues(iSS)-0.3,nSetSizeValues(iSS)+0.3];
    plot(ssAxis,meanResponseTime_PatientSubject_Ord.*ones(length(ssAxis),1),'Color','b','LineWidth',1.5)
    plot(ssAxis,meanResponseTime_PatientSubject_Ord.*ones(length(ssAxis),1)+semResponseTime_PatientSubject_Ord,'Color','b')
    plot(ssAxis,meanResponseTime_PatientSubject_Ord.*ones(length(ssAxis),1)-semResponseTime_PatientSubject_Ord,'Color','b')
    %         legend('Accuracy','Mean','Mean+s.e.m','Mean-s.e.m','Location','NorthWest')
    %         title('Response Time')
end

%% Response time for each patient
nPatientList = unique(PatientSessionList(:,1));
for iPatient = 1:length(nPatientList)
   nPatient = nPatientList(iPatient);
   ind = find(cellfun(@isempty,ResponseTime_Patient_Session_CI{nPatient}));
   ResponseTime_Patient_Session_CI{nPatient}(ind) = [];
   nNumberOfSessions = length(ResponseTime_Patient_Session_CI{nPatient});
   for iDec = 1:2
       ResponseTime_SinglePatient_CI{iDec} = [];
       for iSession = 1:nNumberOfSessions
           ResponseTime_SinglePatient_CI{iDec} = [ResponseTime_SinglePatient_CI{iDec};...
               ResponseTime_Patient_Session_CI{nPatient}{iSession}{iDec}];
       end
   end
   for iDec = 1:2
       MedianResponseTime_Patients(iPatient,iDec) = median(ResponseTime_SinglePatient_CI{iDec});
   end
end

%% Paired t-test for median response times / patients
[h,p] = ttest(MedianResponseTime_Patients(:,1),MedianResponseTime_Patients(:,2));

[mean(MedianResponseTime_Patients(:,1)),mean(MedianResponseTime_Patients(:,2))]

%% Response time for each session
for iDec = 1:2
    ResponseTime_Sessions_CI{iDec} = [];
    for nPatientSubject = 1:nNumberOfPatientsSubjects
        ResponseTime_Sessions_CI{iDec} = [ResponseTime_Sessions_CI{iDec};...
            ResponseTime_PatientSubject_CI{nPatientSubject}{iDec}];
        MedianResponseTime_Sessions(nPatientSubject,iDec) = median(ResponseTime_PatientSubject_CI{nPatientSubject}{iDec});
    end
end

%% Paired t-test for median response times / sessions
[h,p] = ttest(MedianResponseTime_Sessions(:,1),MedianResponseTime_Sessions(:,2))

[nanmean(MedianResponseTime_Sessions(:,1)),nanmean(MedianResponseTime_Sessions(:,2))]

%% Response times for set size / correct incorrect / match mismatch
xVal = [];
xPos = [];
for iDec = 1:2
    for iMM = 1:2
        for iSS = 1:3
            ps = (iSS-1)*4+(iDec-1)*2+iMM;
            xVal = [xVal;squeeze(ResponseTime_Median_PatientSubject_CI_MM_SS(:,iDec,iMM,iSS))];
            xPos = [xPos;ps*ones(size(squeeze(ResponseTime_Median_PatientSubject_CI_MM_SS(:,iDec,iMM,iSS))))];
        end
    end
end

figure
boxplot(xVal,xPos)
ylim([-2,10])
set(gca,'XTick',[])
text(2.5,-1.75,'Set size 4','HorizontalAlignment','Center')
text(1.5,-1,'Correct','HorizontalAlignment','Center')
text(3.5,-1,'Incorrect','HorizontalAlignment','Center')
text(1,-0.25,'In','HorizontalAlignment','Center')
text(2,-0.25,'Out','HorizontalAlignment','Center')
text(3,-0.25,'In','HorizontalAlignment','Center')
text(4,-0.25,'Out','HorizontalAlignment','Center')

text(4+2.5,-1.75,'Set size 6','HorizontalAlignment','Center')
text(4+1.5,-1,'Correct','HorizontalAlignment','Center')
text(4+3.5,-1,'Incorrect','HorizontalAlignment','Center')
text(4+1,-0.25,'In','HorizontalAlignment','Center')
text(4+2,-0.25,'Out','HorizontalAlignment','Center')
text(4+3,-0.25,'In','HorizontalAlignment','Center')
text(4+4,-0.25,'Out','HorizontalAlignment','Center')

text(8+2.5,-1.75,'Set size 8','HorizontalAlignment','Center')
text(8+1.5,-1,'Correct','HorizontalAlignment','Center')
text(8+3.5,-1,'Incorrect','HorizontalAlignment','Center')
text(8+1,-0.25,'In','HorizontalAlignment','Center')
text(8+2,-0.25,'Out','HorizontalAlignment','Center')
text(8+3,-0.25,'In','HorizontalAlignment','Center')
text(8+4,-0.25,'Out','HorizontalAlignment','Center')

ylabel('Reaction time (s)')

%% Response times / correct incorrect / match mismatch
xVal = [];
xPos = [];
for iDec = 1:2
    for iMM = 1:2
        ps = (iDec-1)*2+iMM;
        xVal = [xVal;squeeze(ResponseTime_Median_PatientSubject_CI_MM(:,iDec,iMM))];
        xPos = [xPos;ps*ones(size(squeeze(ResponseTime_Median_PatientSubject_CI_MM(:,iDec,iMM))))];
    end
end

figure
boxplot(xVal,xPos)
ylim([-2,10])
set(gca,'XTick',[])
text(1.5,-1,'Correct','HorizontalAlignment','Center')
text(3.5,-1,'Incorrect','HorizontalAlignment','Center')
text(1,-0.25,'In','HorizontalAlignment','Center')
text(2,-0.25,'Out','HorizontalAlignment','Center')
text(3,-0.25,'In','HorizontalAlignment','Center')
text(4,-0.25,'Out','HorizontalAlignment','Center')

ylabel('Reaction time (s)')

%%
[h,p] = ttest(ResponseTime_Median_PatientSubject_CI_MM(:,1,1),ResponseTime_Median_PatientSubject_CI_MM(:,1,2))
% [h,p] = ttest(ResponseTime_Median_PatientSubject_CI_MM(:,2,1),ResponseTime_Median_PatientSubject_CI_MM(:,2,2))

% addpath(genpath('H:\MATLAB Codes\Toolboxes\eeglab14_1_1b\functions\'))

a = {ResponseTime_Median_PatientSubject_CI_MM(:,1,1)',ResponseTime_Median_PatientSubject_CI_MM(:,1,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000);

pvals

mean(ResponseTime_Median_PatientSubject_CI_MM(:,1,1))
mean(ResponseTime_Median_PatientSubject_CI_MM(:,1,2))
mean(ResponseTime_Mean_PatientSubject_CI_MM(:,1,1))
mean(ResponseTime_Mean_PatientSubject_CI_MM(:,1,2))

%% Response times / match mismatch / only correct
a = {ResponseTime_Mean_PatientSubject_CI_MM(:,1,1)',ResponseTime_Mean_PatientSubject_CI_MM(:,1,2)'};
% a = {ResponseTime_Median_PatientSubject_CI_MM(:,1,1)',ResponseTime_Median_PatientSubject_CI_MM(:,1,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)
[h,p] = ttest(ResponseTime_Mean_PatientSubject_CI_MM(:,1,1),ResponseTime_Mean_PatientSubject_CI_MM(:,1,2))


%% Response times for set size / only correct / match mismatch
xVal = [];
xPos = [];
for iDec = 1:1
    for iMM = 1:2
        for iSS = 1:3
            ps = (iSS-1)*4+(iDec-1)*2+iMM;
            xVal = [xVal;squeeze(ResponseTime_Median_PatientSubject_CI_MM_SS(:,iDec,iMM,iSS))];
            xPos = [xPos;ps*ones(size(squeeze(ResponseTime_Median_PatientSubject_CI_MM_SS(:,iDec,iMM,iSS))))];
        end
    end
end

figure
boxplot(xVal,xPos)
ylim([-2,Inf])
set(gca,'XTick',[])
text(1.5,-1.75,'Set size 4','HorizontalAlignment','Center')
text(1.5,-1,'Correct','HorizontalAlignment','Center')
text(1,-0.25,'In','HorizontalAlignment','Center')
text(2,-0.25,'Out','HorizontalAlignment','Center')

text(2+1.5,-1.75,'Set size 6','HorizontalAlignment','Center')
text(2+1.5,-1,'Correct','HorizontalAlignment','Center')
text(2+1,-0.25,'In','HorizontalAlignment','Center')
text(2+2,-0.25,'Out','HorizontalAlignment','Center')

text(4+1.5,-1.75,'Set size 8','HorizontalAlignment','Center')
text(4+1.5,-1,'Correct','HorizontalAlignment','Center')
text(4+1,-0.25,'In','HorizontalAlignment','Center')
text(4+2,-0.25,'Out','HorizontalAlignment','Center')

ylabel('Reaction time (s)')

%%
[h,p] = ttest(ResponseTime_Median_PatientSubject_CI(:,1),ResponseTime_Median_PatientSubject_CI(:,2))


%% Response times for set size / only correct
xVal = [];
xPos = [];
for iDec = 1:1
        for iSS = 1:3
            ps = (iSS-1)*4+(iDec-1)*2+iMM;
            xVal = [xVal;squeeze(ResponseTime_Median_PatientSubject_CI_SS(:,iDec,iSS))];
            xPos = [xPos;ps*ones(size(squeeze(ResponseTime_Median_PatientSubject_CI_SS(:,iDec,iSS))))];
        end
end

figure
boxplot(xVal,xPos)
ylim([-1,Inf])
set(gca,'XTick',[])
text(1,-1.75+1.5,'Set size 4','HorizontalAlignment','Center')
text(1,-1+1.5,'Correct','HorizontalAlignment','Center')

text(1+1,-1.75+1.5,'Set size 6','HorizontalAlignment','Center')
text(1+1,-1+1.5,'Correct','HorizontalAlignment','Center')

text(2+1,-1.75+1.5,'Set size 8','HorizontalAlignment','Center')
text(2+1,-1+1.5,'Correct','HorizontalAlignment','Center')

ylabel('Reaction time (s)')

%%
[median(ResponseTime_Median_PatientSubject_CI_SS(:,1,1));...
    median(ResponseTime_Median_PatientSubject_CI_SS(:,1,2));...
    median(ResponseTime_Median_PatientSubject_CI_SS(:,1,3))]

%%
[h,p] = ttest(ResponseTime_Median_PatientSubject_CI_SS(:,1,3),ResponseTime_Median_PatientSubject_CI_SS(:,1,1))

[p, h,] = ranksum(ResponseTime_Median_PatientSubject_CI_SS(:,1,1),ResponseTime_Median_PatientSubject_CI_SS(:,1,3))

%% Response times for correct incorrect
xVal = [];
xPos = [];
for iDec = 1:2
    ps = (iDec-1);
    xVal = [xVal;squeeze(ResponseTime_Median_PatientSubject_CI(:,iDec))];
    xPos = [xPos;ps*ones(size(squeeze(ResponseTime_Median_PatientSubject_CI(:,iDec))))];
end

figure
boxplot(xVal,xPos)
ylim([-2,10])
set(gca,'XTick',[])
text(1,-1,'Correct','HorizontalAlignment','Center')
text(2,-1,'Incorrect','HorizontalAlignment','Center')

ylabel('Reaction time (s)')


%%
[h,p] = ttest(ResponseTime_Median_PatientSubject_CI(:,1),ResponseTime_Median_PatientSubject_CI(:,2));
[h,p] = ttest(ResponseTime_Mean_PatientSubject_CI(:,1),ResponseTime_Mean_PatientSubject_CI(:,2))
a = {ResponseTime_Mean_PatientSubject_CI(:,1)',ResponseTime_Mean_PatientSubject_CI(:,2)'};
a{1}(isnan(a{2})) = [];
a{2}(isnan(a{2})) = [];
% a = {ResponseTime_Median_PatientSubject_CI_MM(:,1,1)',ResponseTime_Median_PatientSubject_CI_MM(:,1,2)'};
[stats,df,pvals,surrog] = statcond(a,'method','perm','naccu',2000,'paired','on','verbose','off');
fprintf('\nt-value %.2f; p = %.4f\n',stats,pvals)

