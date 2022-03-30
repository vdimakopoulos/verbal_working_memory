%% Close all figures, clear variables and command window
close all
clear
clc

%% Paths
strPaths.Main = 'F:\Vasileios\';
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\';
% Change main directory
cd(strPaths.Main)
% Add all subfolders to path
addpath(strPaths.Main)
addpath(genpath(strPaths.Project))
addpath(strPaths.Toolboxes.FieldTrip)
% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Default parameters
stat_alpha = 0.05; %Significance Level
Individ_patientFig = 0; %plot figures for each patient individually
%% Load data
strPaths.BehaviouralData = 'F:\Vasileios\Task Analysis\Data\Monsell Task\presentation Data\';
cd(strPaths.BehaviouralData);
nFiles = dir('*.txt');

for i = 1:length(nFiles)
    clear Condition_Trials;
    TrialInformationTable{i} = readtable(nFiles(i).name);
    TrialInformationTable{i}.Properties.VariableNames = {'nTrial','SetSize','Letters','Condition','In_Out','Probe_letter','Correct_response','Correct','Reaction_time'};
    %% compare reaction time and accuracy between the two conditions (Recent  (r) vs non-recent (nr)
    recent_trials = find(ismember(TrialInformationTable{i}.Condition,'rp'));
    non_recent_trials =  find(ismember(TrialInformationTable{i}.Condition,'nrp'));
    
    ConditionReactionTime{i,1}= TrialInformationTable{i}.Reaction_time(recent_trials);
    ConditionReactionTime{i,2}= TrialInformationTable{i}.Reaction_time(non_recent_trials);
    ConditionReactionTimeVec = [ConditionReactionTime{i,1};ConditionReactionTime{i,2}];
    group = [zeros(1,length(ConditionReactionTime{i,1}))';ones(1,length(ConditionReactionTime{i,2}))'];
    [p,h] = ranksum(ConditionReactionTime{i,1},ConditionReactionTime{i,2});
    if Individ_patientFig
        
        figure;
        boxplot(ConditionReactionTimeVec./1000,group);
        if p < stat_alpha
            annotation(gcf,'textbox',...
                [0.460714285714285 0.894238097153961 0.0669642842507788 0.0880952361793745],...
                'String',{'*'},...
                'FontSize',16,...
                'EdgeColor','none');
        else
            annotation(gcf,'textbox',...
                [0.460714285714285 0.894238097153961 0.0669642842507788 0.0880952361793745],...
                'String',{'n.s.'},...
                'FontSize',16,...
                'EdgeColor','none');
        end
        xlabel('Condition');
        ylabel('Reaction time (s)');
        set(gcf,'Color','white');
        set(gca,'FontSize', 16,'XTick',[1 2],'XTickLabel',{'Recent', 'Non-Recent'},'box','off');
    end
    ConditionAccuracy{i}(1)= (length(find(TrialInformationTable{i}.Correct(recent_trials)))/length(recent_trials))*100;
    ConditionAccuracy{i}(2)= (length(find(TrialInformationTable{i}.Correct(non_recent_trials)))/length(non_recent_trials))*100;
    sprintf(['In patient %d accuracy for recent trials is %.1f%% and for non recent trials is %.1f%%'],i, ConditionAccuracy{i}(1), ConditionAccuracy{i}(2))
    
    %% Compare between rp and nrp / rn and nrn
    for iTrl = 1:size(TrialInformationTable{i},1)
        %% NRN = non recent + OUT
        if strcmp(TrialInformationTable{i}.Condition{iTrl},'nrp') && strcmp(TrialInformationTable{i}.In_Out{iTrl},'OUT')
            Condition_Trials{iTrl} = 'NRN';
            %% NRP = non recent + IN
        elseif strcmp(TrialInformationTable{i}.Condition{iTrl},'nrp') && strcmp(TrialInformationTable{i}.In_Out{iTrl},'IN')
            Condition_Trials{iTrl} = 'NRP';
            %% RP = non recent + IN
        elseif strcmp(TrialInformationTable{i}.Condition{iTrl},'rp') && strcmp(TrialInformationTable{i}.In_Out{iTrl},'IN')
            Condition_Trials{iTrl} = 'RP';
            
            %% RN = non recent + OUT
        elseif strcmp(TrialInformationTable{i}.Condition{iTrl},'rp') && strcmp(TrialInformationTable{i}.In_Out{iTrl},'OUT')
            Condition_Trials{iTrl} = 'RN';
            
        end
    end
    TrialInformationTable{i}.ConditionTrial = Condition_Trials';
    
    rn_trials = find(ismember(TrialInformationTable{i}.ConditionTrial,'RN'));
    nrn_trials =  find(ismember(TrialInformationTable{i}.ConditionTrial,'NRN'));
    rp_trials = find(ismember(TrialInformationTable{i}.ConditionTrial,'RP'));
    nrp_trials =  find(ismember(TrialInformationTable{i}.ConditionTrial,'NRP'));
    %RN vs non RN
    ConditionTrialReactionTime{i,1}= TrialInformationTable{i}.Reaction_time(rn_trials);
    ConditionTrialReactionTime{i,2}= TrialInformationTable{i}.Reaction_time(nrn_trials);
    ConditionTrialReactionTimeVec1 = [ConditionTrialReactionTime{i,1};ConditionTrialReactionTime{i,2}];
    group1 = [zeros(1,length(ConditionTrialReactionTime{i,1}))';ones(1,length(ConditionTrialReactionTime{i,2}))'];
    %RP vs NRP
    ConditionTrialReactionTime{i,3}= TrialInformationTable{i}.Reaction_time(rp_trials);
    ConditionTrialReactionTime{i,4}= TrialInformationTable{i}.Reaction_time(nrp_trials);
    ConditionTrialReactionTimeVec2 = [ConditionTrialReactionTime{i,3};ConditionTrialReactionTime{i,4}];
    group2 = [zeros(1,length(ConditionTrialReactionTime{i,3}))';ones(1,length(ConditionTrialReactionTime{i,4}))'];
    [p1,h1] = ranksum(ConditionTrialReactionTime{i,1},ConditionTrialReactionTime{i,2});
    [p2,h2] = ranksum(ConditionTrialReactionTime{i,3},ConditionTrialReactionTime{i,4});
    if Individ_patientFig
        figure;
        boxplot([ConditionTrialReactionTimeVec1; ConditionTrialReactionTimeVec2]./1000,[group1; group2+2]);
        if p1 < stat_alpha
            annotation(gcf,'textbox',...
                [0.280357142857142,0.770428573344439,0.066964284250779,0.088095236179375],...
                'String',{'*'},...
                'FontSize',16,...
                'EdgeColor','none');
        else
            annotation(gcf,'textbox',...
                [0.280357142857142,0.770428573344439,0.066964284250779,0.088095236179375],...
                'String',{'n.s.'},...
                'FontSize',16,...
                'EdgeColor','none');
        end
        
        if p2 < stat_alpha
            annotation(gcf,'textbox',...
                [0.699999999999998,0.768047620963488,0.066964284250779,0.088095236179375],...
                'String',{'*'},...
                'FontSize',16,...
                'EdgeColor','none');
        else
            annotation(gcf,'textbox',...
                [0.699999999999998,0.768047620963488,0.066964284250779,0.088095236179375],...
                'String',{'n.s.'},...
                'FontSize',16,...
                'EdgeColor','none');
        end
        xlabel('Condition');
        ylabel('Reaction time (s)');
        set(gcf,'Color','white');
        set(gca,'FontSize', 14,'XTick',[1:4],'XTickLabel',...
            {'Recent negative', 'Non-recent negative','Recent positive',...
            'Non-recent positive'},'box','off','XTickLabelRotation',45);
    end
    
    ConditionTrialAccuracy{i}(1)= (length(find(TrialInformationTable{i}.Correct(rn_trials)))/length(rn_trials))*100;
    ConditionTrialAccuracy{i}(2)= (length(find(TrialInformationTable{i}.Correct(nrn_trials)))/length(nrn_trials))*100;
    
    ConditionTrialAccuracy{i}(3)= (length(find(TrialInformationTable{i}.Correct(rp_trials)))/length(rp_trials))*100;
    ConditionTrialAccuracy{i}(4)= (length(find(TrialInformationTable{i}.Correct(nrp_trials)))/length(nrp_trials))*100;
    sprintf(['In patient %d accuracy for RN trials is %.1f%% and for NRN trials is %.1f%%'],i, ConditionTrialAccuracy{i}(1), ConditionTrialAccuracy{i}(2))
    sprintf(['In patient %d accuracy for RP trials is %.1f%% and for NRP trials is %.1f%%'],i, ConditionTrialAccuracy{i}(3), ConditionTrialAccuracy{i}(4))
    
    
    %% compare between in and out trials reaction time and accuracy
    IN_trials = find(ismember(TrialInformationTable{i}.In_Out,'IN'));
    OUT_trials =  find(ismember(TrialInformationTable{i}.In_Out,'OUT'));
    
    In_OutReactionTime{i,1}= TrialInformationTable{i}.Reaction_time(IN_trials);
    In_OutReactionTime{i,2}= TrialInformationTable{i}.Reaction_time(OUT_trials);
    In_OutReactionTimeVec = [In_OutReactionTime{i,1};In_OutReactionTime{i,2}];
    group = [zeros(1,length(In_OutReactionTime{i,1}))';ones(1,length(In_OutReactionTime{i,2}))'];
    [p,h] = ranksum(In_OutReactionTime{i,1},In_OutReactionTime{i,2});
    if Individ_patientFig
        
        figure;
        boxplot(In_OutReactionTimeVec./1000,group);
        if p < stat_alpha
            annotation(gcf,'textbox',...
                [0.460714285714285 0.894238097153961 0.0669642842507788 0.0880952361793745],...
                'String',{'*'},...
                'FontSize',16,...
                'EdgeColor','none');
        else
            annotation(gcf,'textbox',...
                [0.460714285714285 0.894238097153961 0.0669642842507788 0.0880952361793745],...
                'String',{'n.s.'},...
                'FontSize',16,...
                'EdgeColor','none');
        end
        xlabel('IN/OUT');
        ylabel('Reaction time (s)');
        set(gcf,'Color','white');
        set(gca,'FontSize', 16,'XTick',[1 2],'XTickLabel',{'IN', 'OUT'},'box','off');
    end
    InOutAccuracy{i}(1)= (length(find(TrialInformationTable{i}.Correct(IN_trials)))/length(IN_trials))*100;
    InOutAccuracy{i}(2)= (length(find(TrialInformationTable{i}.Correct(OUT_trials)))/length(OUT_trials))*100;
    sprintf(['In patient %d accuracy for IN trials is %.1f%% and for OUT trials is %.1f%%'],i, InOutAccuracy{i}(1), InOutAccuracy{i}(2))
    
    %%
    %% compare  reaction time and accuracy among different set size
    SetSize_trials{1} = find(ismember(TrialInformationTable{i}.SetSize,4));
    SetSize_trials{2} =  find(ismember(TrialInformationTable{i}.SetSize,6));
    SetSize_trials{3} =  find(ismember(TrialInformationTable{i}.SetSize,8));
    SetSize_trials{4} =  find(ismember(TrialInformationTable{i}.SetSize,[6 8]));
    SetSize_trials{5} =  find(ismember(TrialInformationTable{i}.SetSize,[4 6 8]));
    SetSizeReactionTimeVec = [];
    for ss = 1:length(SetSize_trials)-2 % Only trials with set size 4/6/8 letters
        SetSizeReactionTime{ss,i}= TrialInformationTable{i}.Reaction_time(SetSize_trials{ss});
        SetSizeReactionTimeVec = [SetSizeReactionTimeVec;SetSizeReactionTime{ss,i}];
    end
    group = [zeros(1,length(SetSizeReactionTime{1,i}))';ones(1,length(SetSizeReactionTime{2,i}))';2*ones(1,length(SetSizeReactionTime{3,i}))'];
    if Individ_patientFig
        figure;
        boxplot(SetSizeReactionTimeVec'./1000,group);
        xlabel('Set Size');
        ylabel('Reaction time (s)');
        set(gcf,'Color','white');
        set(gca,'FontSize', 16,'XTick',[1:3],'XTickLabel',{'4', '6','8'},'box','off');
    end
    SetSizeAccuracy{i}(1)= (length(find(TrialInformationTable{i}.Correct(SetSize_trials{1})))/length(SetSize_trials{1}))*100;
    SetSizeAccuracy{i}(2)= (length(find(TrialInformationTable{i}.Correct(SetSize_trials{2})))/length(SetSize_trials{2}))*100;
    SetSizeAccuracy{i}(3)= (length(find(TrialInformationTable{i}.Correct(SetSize_trials{3})))/length(SetSize_trials{3}))*100;
    
    sprintf(['In patient %d accuracy for 4 letters is %.1f%%,  for 6 letters is %.1f%% and for 8 letters is %.1f%%'],i, SetSizeAccuracy{i}(1), SetSizeAccuracy{i}(2),SetSizeAccuracy{i}(3))
    
    
    
    
    
    
end

%% Accuracy/reaction time for each condition in the group level
%reaction time
%RN
ConditionTrialRT{1} = cell2mat({cellfun(@median,{ConditionTrialReactionTime{:,1}}','UniformOutput',1)});
%NRN
ConditionTrialRT{2} = cell2mat({cellfun(@median,{ConditionTrialReactionTime{:,2}}','UniformOutput',1)});
%RP
ConditionTrialRT{3} = cell2mat({cellfun(@median,{ConditionTrialReactionTime{:,3}}','UniformOutput',1)});
%NRP
ConditionTrialRT{4} = cell2mat({cellfun(@median,{ConditionTrialReactionTime{:,4}}','UniformOutput',1)});

ConditionTrialRTVec =cell2mat(ConditionTrialRT');
group = [zeros(1,length(ConditionTrialRT{1}))';ones(1,length(ConditionTrialRT{2}))';ones(1,length(ConditionTrialRT{3}))'*2;ones(1,length(ConditionTrialRT{4}))'*3];
figure;
boxplot(ConditionTrialRTVec./1000,group);
hold on;
scatter(ones(length(ConditionTrialRT{1}),1),ConditionTrialRT{1}./1000,'filled','r');
scatter(ones(length(ConditionTrialRT{2}),1)*2,ConditionTrialRT{2}./1000,'filled','g');
scatter(ones(length(ConditionTrialRT{3}),1)*3,ConditionTrialRT{3}./1000,'filled','b');
scatter(ones(length(ConditionTrialRT{4}),1)*4,ConditionTrialRT{4}./1000,'filled','k');

xlabel('Condition');
ylabel('Reaction Time (s)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:4],'XTickLabel',{'RN', 'NRN','RP','NRP'},'box','off');


%accuracy
temp_CondAcc = cell2mat(ConditionTrialAccuracy);
ConditionAccParticipants(:,1) = temp_CondAcc(1:4:end); %rn
ConditionAccParticipants(:,2) = temp_CondAcc(2:4:end); %nrn
ConditionAccParticipants(:,3) = temp_CondAcc(3:4:end);%rp
ConditionAccParticipants(:,4) = temp_CondAcc(4:4:end);%nrp


figure;
boxplot(ConditionAccParticipants);
hold on;
scatter(ones(length(ConditionAccParticipants(:,1)),1),ConditionAccParticipants(:,1),'filled','r');
scatter(ones(length(ConditionAccParticipants(:,2)),1)*2,ConditionAccParticipants(:,2),'filled','g');
scatter(ones(length(ConditionAccParticipants(:,3)),1)*3,ConditionAccParticipants(:,3),'filled','b');
scatter(ones(length(ConditionAccParticipants(:,4)),1)*4,ConditionAccParticipants(:,4),'filled','k');


xlabel('Condition');
ylabel('Accuracy (%)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:4],'XTickLabel',{'RN', 'NRN','RP', 'NRP'},'box','off');
%% IN/OUT accuracy in the group of patients
%reaction time
%IN
InOutRTAcc{1} = cell2mat({cellfun(@median,{In_OutReactionTime{:,1}}','UniformOutput',1)});
%OUT
InOutRTAcc{2} = cell2mat({cellfun(@median,{In_OutReactionTime{:,2}}','UniformOutput',1)});
InOutRTAccVec =cell2mat(InOutRTAcc');
group = [zeros(1,length(InOutRTAcc{1}))';ones(1,length(InOutRTAcc{2}))'];
figure;
boxplot(InOutRTAccVec./1000,group);
hold on;
scatter(ones(length(InOutRTAcc{1}),1),InOutRTAcc{1}./1000,'filled','r');
scatter(ones(length(InOutRTAcc{2}),1)*2,InOutRTAcc{2}./1000,'filled','g');

xlabel('IN OUT');
ylabel('Reaction Time (s)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:2],'XTickLabel',{'IN', 'OUT'},'box','off');



%accuracy
temp_IO = cell2mat(InOutAccuracy);
In_OutAccParticipants(:,1) = temp_IO(1:2:end);
In_OutAccParticipants(:,2) = temp_IO(2:2:end);

figure;
boxplot(In_OutAccParticipants);
hold on;
scatter(ones(length(In_OutAccParticipants(:,1)),1),In_OutAccParticipants(:,1),'filled','r');
scatter(ones(length(In_OutAccParticipants(:,2)),1)*2,In_OutAccParticipants(:,2),'filled','g');


xlabel('In-Out');
ylabel('Accuracy (%)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:2],'XTickLabel',{'IN', 'OUT'},'box','off');

%% Set Size Reaction time and accuracy in the group of patients
%Reaction time
%ss 4
GroupSSAcc{1} = cell2mat({cellfun(@median,{SetSizeReactionTime{1,:}}','UniformOutput',1)});
%cell2mat({SetSizeReactionTime{1,:}}');
%ss 6
GroupSSAcc{2} = cell2mat({cellfun(@median,{SetSizeReactionTime{2,:}}','UniformOutput',1)});
%ss8
GroupSSAcc{3} = cell2mat({cellfun(@median,{SetSizeReactionTime{3,:}}','UniformOutput',1)});
GroupSSAccVec =cell2mat(GroupSSAcc');
group = [zeros(1,length(GroupSSAcc{1}))';ones(1,length(GroupSSAcc{2}))';2*ones(1,length(GroupSSAcc{3}))'];
figure;
boxplot(GroupSSAccVec./1000,group);
hold on;
scatter(ones(length(GroupSSAcc{1}),1),GroupSSAcc{1}./1000,'filled','r');
scatter(ones(length(GroupSSAcc{2}),1)*2,GroupSSAcc{2}./1000,'filled','g');
scatter(ones(length(GroupSSAcc{3}),1)*3,GroupSSAcc{3}./1000,'filled','b');

xlabel('Set Size');
ylabel('Reaction Time (s)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:3],'XTickLabel',{'4', '6','8'},'box','off');


% Accuracy
temp_SS = cell2mat(SetSizeAccuracy);
SetSizeAccParticipants(:,1) = temp_SS(1:3:end);
SetSizeAccParticipants(:,2) = temp_SS(2:3:end);
SetSizeAccParticipants(:,3) = temp_SS(3:3:end);

figure;
boxplot(SetSizeAccParticipants);
hold on;
scatter(ones(length(SetSizeAccParticipants(:,1)),1),SetSizeAccParticipants(:,1),'filled','r');
scatter(ones(length(SetSizeAccParticipants(:,2)),1)*2,SetSizeAccParticipants(:,2),'filled','g');
scatter(ones(length(SetSizeAccParticipants(:,3)),1)*3,SetSizeAccParticipants(:,3),'filled','b');


xlabel('Set Size');
ylabel('Accuracy (%)');
set(gcf,'Color','white');
set(gca,'FontSize', 16,'XTick',[1:3],'XTickLabel',{'4', '6','8'},'box','off');
