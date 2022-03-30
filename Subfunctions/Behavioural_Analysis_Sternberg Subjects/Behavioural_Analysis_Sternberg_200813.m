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
addpath(genpath(strPaths.Project))
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

% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Subject to load
PatientList = [10 12 14 20 26 37 38 40 42 44 45 ]
nPid = PatientList(1:end);
strPaths.Subject_Data = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\';

for i = 1:numel(nPid)
    switch nPid(i)
        case 10
            strPaths.PathToLoad = [strPaths.Subject_Data,'10 SO\'];
        case 12
            strPaths.PathToLoad = [strPaths.Subject_Data,'12 BP\'];
        case 14
            strPaths.PathToLoad = [strPaths.Subject_Data,'14 UC\'];
        case 20
            strPaths.PathToLoad = [strPaths.Subject_Data,'20 SN\'];
        case 26
            strPaths.PathToLoad = [strPaths.Subject_Data,'26 MC\'];
        case 37
            strPaths.PathToLoad = [strPaths.Subject_Data,'37 PN\'];
        case 38
            strPaths.PathToLoad = [strPaths.Subject_Data,'38 NP\'];
        case 40
            strPaths.PathToLoad = [strPaths.Subject_Data,'40 DG\'];
        case 42
            strPaths.PathToLoad = 'F:\Vasileios\Task Analysis\Data\Macro Data\';
        case 44
            strPaths.PathToLoad = [strPaths.Subject_Data,'44 MJ\'];
        case 45
            strPaths.PathToLoad = [strPaths.Subject_Data,'45 SS\'];
    end
    
    
    cd (strPaths.PathToLoad);
    files = dir('*.mat'); %try to look for all the .mat files under the folder
    for ii=1:length(files)
        data_ses(ii) = load(files(ii).name); %load the files
    end
    TrialInformationTable{i} = []
    for ij=1:length(data_ses)
        TrialInformationTable{i} = [TrialInformationTable{i};data_ses(ij).TrialInformationTable];
    end
    
    clear data_ses;
    cd(strPaths.Project);
end


%% Correct Response Rate (Average, IN, Out)
for i = 1:length(TrialInformationTable)
    Correct_Response_Trials{i} = find(TrialInformationTable{i}.Correct == 1);
    nCorrect  = numel(Correct_Response_Trials{i});
    %Average Correct Response rate:
    Correct_Response_Rate{i} =  nCorrect/size(TrialInformationTable{i},1)*100;   
    In_Responses = length(find(TrialInformationTable{i}.Response == 51))+length(find(TrialInformationTable{i}.Response == 1))
    Correct_In = length(find(TrialInformationTable{i}.Response == 51));
    
    Out_Responses = length(find(TrialInformationTable{i}.Response == 52))+length(find(TrialInformationTable{i}.Response == 2));
    Correct_out = length(find(TrialInformationTable{i}.Response == 52));
    
    %Average Correct response rate for IN:
    In_Rate{i} = Correct_In/In_Responses*100; 
    %Average Correct response rate for OUT:
    Out_Rate{i} = Correct_out/Out_Responses*100; 
end

%% Response Rate for Set Sizes
for i = 1:length(TrialInformationTable)
    nSetSize_4{i} = find(TrialInformationTable{i}.SetSize == 4);
    nSetSize_6{i} = find(TrialInformationTable{i}.SetSize == 6);
    nSetSize_8{i} = find(TrialInformationTable{i}.SetSize == 8);
    
    Correct_SS4{i} = intersect(Correct_Response_Trials{i},nSetSize_4{i});
    Correct_SS6{i} = intersect(Correct_Response_Trials{i},nSetSize_6{i});
    Correct_SS8{i} = intersect(Correct_Response_Trials{i},nSetSize_8{i});
    
    %Average Correct response rate for set size 4:
    SS4_rate{i} =  (length(Correct_SS4{i})/length(nSetSize_4{i}))*100; 
    %Average Correct response rate for set size 6:
    SS6_rate{i} =   (length(Correct_SS6{i})/length(nSetSize_6{i}))*100;
    %Average Correct response rate for set size 8:
    SS8_rate{i} =  (length(Correct_SS8{i})/length(nSetSize_8{i}))*100; 
end


%% Memory Capacity for Each Subject
for i = 1:length(TrialInformationTable)
    CorrectIN = find(TrialInformationTable{i}.Response == 51);
    CorrectOut = find(TrialInformationTable{i}.Response == 52);
    
    IncorrectIN = find(TrialInformationTable{i}.Response == 1);
    IncorrectOut = find(TrialInformationTable{i}.Response == 2);
    
    IN = [CorrectIN; IncorrectIN];
    OUT = [CorrectOut; IncorrectOut];
    
    In_SS4 = intersect(IN,find(TrialInformationTable{i}.SetSize == 4));
    In_SS6 = intersect(IN,find(TrialInformationTable{i}.SetSize == 6));
    In_SS8 = intersect(IN,find(TrialInformationTable{i}.SetSize == 8));
    
    Out_SS_4 = intersect(OUT,find(TrialInformationTable{i}.SetSize == 4));
    Out_SS_6 = intersect(OUT,find(TrialInformationTable{i}.SetSize == 6));
    Out_SS_8 = intersect(OUT,find(TrialInformationTable{i}.SetSize == 8));

    CorrectIN_ss4{i} = length(intersect(CorrectIN, nSetSize_4{i}))/numel(In_SS4);
    CorrectIN_ss6{i} = length(intersect(CorrectIN, nSetSize_6{i}))/numel(In_SS6);
    CorrectIN_ss8{i} = length(intersect(CorrectIN, nSetSize_8{i}))/numel(In_SS8);
    
    CorrectOut_ss4{i} = length(intersect(CorrectOut, nSetSize_4{i}))/numel(Out_SS_4);
    CorrectOut_ss6{i} = length(intersect(CorrectOut, nSetSize_6{i}))/numel(Out_SS_6);
    CorrectOut_ss8{i} = length(intersect(CorrectOut, nSetSize_8{i}))/numel(Out_SS_8);
    
   Cowan_K_SS{i}(1) = (CorrectIN_ss4{i} + CorrectOut_ss4{i} - 1)*4; % 
   Cowan_K_SS{i}(2) = (CorrectIN_ss6{i} + CorrectOut_ss6{i} - 1)*6; %
   %Memory Capacity:
   Cowan_K_SS{i}(3) = (CorrectIN_ss8{i} + CorrectOut_ss8{i} - 1)*8; %
   Memory_Capacity{i} = Cowan_K_SS{i}(3);
end

%% Response Time for each set size
for i = 1:length(TrialInformationTable)
    
    Mean_Correct_Response_Time{i} = median(TrialInformationTable{i}.ResponseTime(Correct_Response_Trials{i}));
    STD_Correct_Response_Time{i} = std(TrialInformationTable{i}.ResponseTime(Correct_Response_Trials{i}));
    
    SS4_corr_time{i} = median(TrialInformationTable{i}.ResponseTime(nSetSize_4{i}));
    SS4_corr_time_std{i} = std(TrialInformationTable{i}.ResponseTime(nSetSize_4{i}));
    
    SS6_corr_time{i} = median(TrialInformationTable{i}.ResponseTime(nSetSize_6{i}));
    SS6_corr_time_std{i} = std(TrialInformationTable{i}.ResponseTime(nSetSize_6{i}));
    
    SS8_corr_time{i} = median(TrialInformationTable{i}.ResponseTime(nSetSize_8{i}));
    SS8_corr_time_std{i} = std(TrialInformationTable{i}.ResponseTime(nSetSize_8{i}));
end


%% Correct IN/OUT Decision Times vs Incorrect
for i = 1:length(TrialInformationTable)
    IncorrectTrials{i} = setdiff([1:size(TrialInformationTable{i},1)],Correct_Response_Trials{i});
    Mean_Incorrect_Response_Time{i} =  median(TrialInformationTable{i}.ResponseTime(IncorrectTrials{i}));
    STD_Incorrect_Response_Time{i}  = std(TrialInformationTable{i}.ResponseTime(IncorrectTrials{i}));
    
end

%% Store the results in a table
Patient_Behaviour_Results_Table = array2table(zeros(length(PatientList),1),'VariableNames',{'Patient_Number'});
Patient_Behaviour_Results_Table.Patient_Number = PatientList(1:end)';
Patient_Behaviour_Results_Table.Patient_Initials = {'SO','BP','UC','SN','MC','PN','NP','DG','DS','MJ','SS' }'
Patient_Behaviour_Results_Table.Avg_Correct_Response_Rate =  Correct_Response_Rate(1:end)'
Patient_Behaviour_Results_Table.In_Rate = In_Rate(1:end)'
Patient_Behaviour_Results_Table.Out_Rate = Out_Rate(1:end)'
Patient_Behaviour_Results_Table.SS4_Correct_Rate =  SS4_rate(1:end)'
Patient_Behaviour_Results_Table.SS6_Correct_Rate =  SS6_rate(1:end)'
Patient_Behaviour_Results_Table.SS8_Correct_Rate =  SS8_rate(1:end)'
Patient_Behaviour_Results_Table.Memory_Capacity = Memory_Capacity(1:end)'
Patient_Behaviour_Results_Table.Mean_Correct_Response_Time = Mean_Correct_Response_Time(1:end)'
Patient_Behaviour_Results_Table.Mean_Incorrect_Response_Time = Mean_Incorrect_Response_Time(1:end)'
Patient_Behaviour_Results_Table.SS4_response_time = SS4_corr_time(1:end)'
Patient_Behaviour_Results_Table.SS6_response_time = SS6_corr_time(1:end)'
Patient_Behaviour_Results_Table.SS8_response_time = SS8_corr_time(1:end)'


filename = 'PatientBehaviourResultsTable.xlsx';
strPaths.PathToSave = ['F\',filename]
writetable(Patient_Behaviour_Results_Table,filename,'Sheet','Behaviour Results','WriteVariableNames',true);


%% Plot Set Sizes accuracy
strPaths.BehaviourResults = [strPaths.Results,'Behaviour Results\Figs\'];
mkdir(strPaths.BehaviourResults);
nList = [1:numel(nPid)]  %[4 6 7 9 11] ;% 
SetSizes = {ones(numel(nList),1)*4 ones(numel(nList),1)*6 ones(numel(nList),1)*8}; %replace nList with nPid for all patients
SetSizesAccuracy = [cell2mat(SS4_rate(nList))' cell2mat(SS6_rate(nList))' cell2mat(SS8_rate(nList))']
SetSizeReactionTime = [cell2mat(SS4_corr_time(nList))' cell2mat(SS6_corr_time(nList))' cell2mat(SS8_corr_time(nList))'];

figure;
xoff =0.15;
nLineWidth_Behavior = 2;
Color = {'bo','go','ro'};
MarkerColor = {'b','g','r'};
% c = 0.05;
% indSort = [1:length(nPid)]*2;
fig = figure;
for iSS = 1:length(SetSizes)
    plot([SetSizes{iSS},SetSizes{iSS}],...
        [SetSizesAccuracy(:,iSS),SetSizesAccuracy(:,iSS)],...
        Color{iSS},'MarkerFaceColor',MarkerColor{iSS},'MarkerSize',3);
    hold on;
    xlim([3 9])
    ylim([60 100])
    plot([SetSizes{iSS}(1)-xoff*1.6,SetSizes{iSS}(1)+xoff*1.6],...
        [median(SetSizesAccuracy(:,iSS)),...
        median(SetSizesAccuracy(:,iSS))],...
        'b','LineWidth',nLineWidth_Behavior+1,'Color',[0.4,0.4,1]);
    plot([SetSizes{iSS}(1)-xoff*2,SetSizes{iSS}(1)+xoff*2],...
        [median(SetSizesAccuracy(:,iSS))+std(SetSizesAccuracy(:,iSS))/size(nPid,2),...
        median(SetSizesAccuracy(:,iSS))+std(SetSizesAccuracy(:,iSS))/size(nPid,2)],...
        'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    plot([SetSizes{iSS}(1)-xoff*2,SetSizes{iSS}(1)+xoff*2],...
        [median(SetSizesAccuracy(:,iSS))-std(SetSizesAccuracy(:,iSS))/size(nPid,2),...
        median(SetSizesAccuracy(:,iSS))-std(SetSizesAccuracy(:,iSS))/size(nPid,2)],...
        'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    
end
hold off;
set(gcf,'color','white')
set(gca,'box','off','Ytick',[60 70 80 90 100],...
    'XTick',[4 6 8],'FontSize',12, 'TickDir','out')
xlabel('Set Size');
ylabel('Accuracy (%)');
saveas(fig, [strPaths.BehaviourResults,'SetSize_Accuracy']);
print(fig,[strPaths.BehaviourResults,'SetSize_Accuracy'],'-dpng','-r1000');
print(fig,[strPaths.BehaviourResults,'SetSize_Accuracy'],'-dpdf','-r1000');
print(fig,[strPaths.BehaviourResults,'SetSize_Accuracy'],'-deps','-r1000');

%Set Size - Reaction time 
fig = figure;
for iSS = 1:length(SetSizes)
    plot([SetSizes{iSS},SetSizes{iSS}],...
        [SetSizeReactionTime(:,iSS),SetSizeReactionTime(:,iSS)],...
        Color{iSS},'MarkerFaceColor',MarkerColor{iSS},'MarkerSize',3);
    hold on;
    xlim([3 9])
    ylim([0.4 3])
    plot([SetSizes{iSS}(1)-xoff*1.6,SetSizes{iSS}(1)+xoff*1.6],...
        [median(SetSizeReactionTime(:,iSS)),...
        median(SetSizeReactionTime(:,iSS))],...
        'b','LineWidth',nLineWidth_Behavior+1,'Color',[0.4,0.4,1]);
    plot([SetSizes{iSS}(1)-xoff*2,SetSizes{iSS}(1)+xoff*2],...
        [median(SetSizeReactionTime(:,iSS))+std(SetSizeReactionTime(:,iSS))/size(nPid,2),...
        median(SetSizeReactionTime(:,iSS))+std(SetSizeReactionTime(:,iSS))/size(nPid,2)],...
        'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    plot([SetSizes{iSS}(1)-xoff*2,SetSizes{iSS}(1)+xoff*2],...
        [median(SetSizeReactionTime(:,iSS))-std(SetSizeReactionTime(:,iSS))/size(nPid,2),...
        median(SetSizeReactionTime(:,iSS))-std(SetSizeReactionTime(:,iSS))/size(nPid,2)],...
        'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    
end
hold off;
set(gcf,'color','white')
set(gca,'box','off','Ytick',[0.8 1.2 1.6 2 2.4 2.8 ],...
    'XTick',[4 6 8],'FontSize',12, 'TickDir','out')
xlabel('Set Size');
ylabel('Reaction time (s)');
 
saveas(fig, [strPaths.BehaviourResults,'SetSize_ReactionTime']);
print(fig,[strPaths.BehaviourResults,'SetSize_ReactionTime'],'-dpng','-r1000');
print(fig,[strPaths.BehaviourResults,'SetSize_ReactionTime'],'-dpdf','-r1000');
print(fig,[strPaths.BehaviourResults,'SetSize_ReactionTime'],'-deps','-r1000');
