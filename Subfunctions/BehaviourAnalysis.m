addpath(genpath('F:\Vasileios\Toolboxes\eeglab14_1_1b\functions'))

Path = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Performance Data\42 DS\';

strFilename_Session1 = [Path,'gelflog 1.txt']
strFilename_Session2 = [Path,'gelflog 2.txt']

Behav_res_s1 = readtable(strFilename_Session1,'Format','auto');
Behav_res_s2 = readtable(strFilename_Session2,'Format','auto');
%Remove the column with the boolean flag of response
Behav_res_s2 = removevars(Behav_res_s2,8);
Behav_res_s2.Properties.VariableNames{8} = 'Var8';
%Fix trial numbers
Behav_res_s2.Var1 = Behav_res_s2.Var1+50; 
Behav_res_Merged = [Behav_res_s1; Behav_res_s2]

%Assign the correct names on the vars
Behav_res_Merged.Properties.VariableNames = {'Trial Number', 'Set Size', 'Sequence', 'In', 'Out',...
    'Probed Letter','Response','Time'}

OutputFile = [Path, 'gelflog_merged.mat']
save(OutputFile,'Behav_res_Merged')

%% Performance for correct response on the IN/OUT case
indIn = find(strcmp(Behav_res_Merged.Properties.VariableNames,'In'));
indOut = find(strcmp(Behav_res_Merged.Properties.VariableNames,'Out'));
response = find(~cellfun(@isempty,strfind(Behav_res_Merged.Response,'right')));
wrong_response = setdiff([1:length(Behav_res_Merged.("Trial Number"))],response);
IN = find(table2array(Behav_res_Merged(:,indIn))==1);

Correct_IN_response = intersect(IN,response);
IncorrectIN = intersect(wrong_response,IN);
Performance_IN = (length(Correct_IN_response)/length(IN))*100;

OUT = find(table2array(Behav_res_Merged(:,indOut))==1);
IncorrectOut = intersect(wrong_response,OUT);
Correct_OUT_response = intersect(OUT,response);
Performance_Out = (length(Correct_OUT_response)/length(OUT))*100;

sprintf('The correct response rate was %.f%% for IN  and %.f%% for OUT',...
    Performance_IN,Performance_Out)
Cowan_K = ((Performance_IN/100 + Performance_Out/100) -1)*[4 6 8];
sprintf('Across all sessions the capacity averaged %.2f which indicates that \nthe subject was able to maintain at least %d letters in memory.',Cowan_K(3), floor(Cowan_K(3)))
%% Set Size Performance 

SS_4 = find(Behav_res_Merged.("Set Size")==4);
SS_4_correct = intersect(SS_4,response);
SS_6 = find(Behav_res_Merged.("Set Size")==6);
SS_6_correct = intersect(SS_6,response);
SS_8 = find(Behav_res_Merged.("Set Size")==8);
SS_8_correct = intersect(SS_8,response);

Performance_SS4 = (length(SS_4_correct)/length(SS_4))*100;
Performance_SS6 = (length(SS_6_correct)/length(SS_6))*100;
Performance_SS8 = (length(SS_8_correct)/length(SS_8))*100;


sprintf('The  rate of correct responses decreased with set size\nfrom a set size of 4 (%.f%% correct responses),\n to set sizes of 6 (%.f%%) and 8 (%.f%%)',...
    Performance_SS4,Performance_SS6,Performance_SS8)
%% Cowan K capacity per set size
Correct_SS_IN_rate = {};
Correct_SS_IN_rate{1} = length(intersect(SS_4,Correct_IN_response))/length(intersect(SS_4,IN));
Correct_SS_IN_rate{2} = length(intersect(SS_6,Correct_IN_response))/length(intersect(SS_6,IN));
Correct_SS_IN_rate{3} = length(intersect(SS_8,Correct_IN_response))/length(intersect(SS_8,IN));

Correct_SS_Out_rate = {};
Correct_SS_Out_rate{1} = length(intersect(SS_4,Correct_OUT_response))/length(intersect(SS_4,OUT));
Correct_SS_Out_rate{2} = length(intersect(SS_6,Correct_OUT_response))/length(intersect(SS_6,OUT));
Correct_SS_Out_rate{3} = length(intersect(SS_8,Correct_OUT_response))/length(intersect(SS_8,OUT));

Cowan_K_SS = {}
Cowan_K_SS{1} = (Correct_SS_IN_rate{1}+Correct_SS_Out_rate{1} - 1)*4 % set size 4
Cowan_K_SS{2} = (Correct_SS_IN_rate{2}+Correct_SS_Out_rate{2} - 1)*6 % set size 6
Cowan_K_SS{3} = (Correct_SS_IN_rate{3}+Correct_SS_Out_rate{3} - 1)*8 % set size 8



%% Set Size Response Time
Corr_TrialNum  = sort([Correct_IN_response',Correct_OUT_response'])
CorrectTrialsTime = Behav_res_Merged.Time(Corr_TrialNum);

SS4_corr_time = Behav_res_Merged.Time(SS_4_correct)/1000;
SS6_corr_time = Behav_res_Merged.Time(SS_6_correct)/1000;
SS8_corr_time = Behav_res_Merged.Time(SS_8_correct)/1000;
sprintf('The mean response time (RT) for correct trials (75 out of 100) is %.2f +- %.2f seconds and is increased with the workload from set size 4 (%.2f +-%.2f s) to 6 ((%.2f +-%.2f s)) and 8((%.2f +-%.2f s) respectively.)', ...
    median(CorrectTrialsTime)/1000,std(CorrectTrialsTime)/1000, median(SS4_corr_time),std(SS4_corr_time), median(SS6_corr_time), ...
    std(SS6_corr_time),median(SS8_corr_time),std(SS8_corr_time))

%% Slope for reaction time per item
SS_correct = sort([SS_4_correct' SS_6_correct' SS_8_correct']);
varSS = Behav_res_Merged.("Set Size")(SS_correct);
varRT = Behav_res_Merged.Time(SS_correct)./1000; %in seconds

LM = fitlm(varSS,varRT,'linear');
figure;
yaxis_t = 0:5;
boxplot(varRT',varSS')
hold on;
% plot(yaxis_t,LM.Coefficients.Estimate(2)*yaxis_t+LM.Coefficients.Estimate(1),'r')
slope = LM.Coefficients.Estimate(2)*1000

%% Decision Times

CorrectIN_time = Behav_res_Merged.Time(Correct_IN_response);
IncorrectIN_time = Behav_res_Merged.Time(IncorrectIN);
MeanCorrectIN_time = mean(CorrectIN_time);
StdCorrectIN_time = std(CorrectIN_time);
MeanIncorrectIN_time = mean(IncorrectIN_time);
StdIncorrectIN_time = std(IncorrectIN_time);

CorrectOut_time = Behav_res_Merged.Time(Correct_OUT_response);
IncorrectOut_time = Behav_res_Merged.Time(IncorrectOut);
MeanCorrectOut_time = mean(CorrectOut_time);
StdCorrectOut_time = std(CorrectOut_time);
MeanIncorrectOut_time = mean(IncorrectOut_time);
StdIncorrectOut_time = std(IncorrectOut_time);

CorrectResponse_Time = [CorrectIN_time; CorrectOut_time];
IncorrectResponse_Time = [IncorrectIN_time;IncorrectOut_time];



sprintf('Correct IN/OUT decisions were made more rapidly than incorrect decisions (%.2f +- %.2f versus %.2f +- %.2f s) respectively.',...
    median(CorrectResponse_Time)/1000,std(CorrectResponse_Time)/1000, median(IncorrectResponse_Time)/1000,...
    std(IncorrectResponse_Time)/1000)


sprintf('Correct IN responses were made more rapidly in average  %.2f  +- %.1f seconds \n versus incorrect responses %.2f +- %.1f s',...
    MeanCorrectIN_time/1000,StdCorrectIN_time/1000 ,MeanIncorrectIN_time/1000,StdIncorrectIN_time/1000)


sprintf('Correct OUT responses were made slower in average  %.2f  +- %.1f seconds \n versus incorrect responses %.2f +- %.1f s',...
    MeanCorrectOut_time/1000,StdCorrectOut_time/1000 ,MeanIncorrectOut_time/1000,StdIncorrectIN_time/1000)




