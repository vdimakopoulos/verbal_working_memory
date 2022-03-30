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

%% Load data for 6 sessions - Macro Data
patient = '37 PN';

strMacroDataDir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\37 PN\';
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end
% EDF data
strEDFdir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\EDF Data\';
cd (strEDFdir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_sesEDF(i) = load(files(i).name); %load the files
    data_sesEDF(i).dataEDF.label = strrep(data_sesEDF(i).dataEDF.label,' ','')
end


%% Load data for 6 sessions - Scalp Data
strScalpDataDir = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\37 PN\';
cd (strScalpDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
 data_Scalp_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable_Scalp = [];
for i=1:length(data_Scalp_ses)
    TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
end


%% Merge Sessions
dataBipolar = data_sesEDF(1).dataEDF;
for nSes = 2:length(data_sesEDF)
    cfg = [];
    dataBipolar = ft_appenddata(cfg,dataBipolar,data_sesEDF(nSes).dataEDF);
    
end



dataBipolarScalp = data_Scalp_ses(1).data;
for nSes = 2:length(data_Scalp_ses)
    if nSes == 6 %% remove the two redundant channels which are present only on this session
        cfg =[];
        cfg.channel = {'all','-Subm1','-Subm2'};
        data_Scalp_ses(nSes).data = ft_selectdata(cfg, data_Scalp_ses(nSes).data);
    end
    cfg = [];
    dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
    
    
end

%% Split to scalp and macro data
cfg = [];
cfg.channel = [1:19];
dataBipolarScalp = ft_selectdata(cfg,dataBipolar);

cfg = [];
cfg.channel = [20:65];
dataBipolar = ft_selectdata(cfg,dataBipolar);


%% Reference
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =  [4] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

A1_chan = find(contains(dataBipolarScalp.label,'A1'));
A2_chan = find(contains(dataBipolarScalp.label,'A2'));
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =  [A1_chan A2_chan] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolarScalp  = ft_preprocessing(cfg,dataBipolarScalp);


%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 15;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 33 34 33 35 34 35] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end
%Append the bipolar channel labels to the reference channels labels
for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
montage.tra             = montage_matrix;
% dataBipolar = ft_apply_montage(dataBipolar,montage);

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod  = 'bipolar'
cfg.refchannel = [47 61]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;


%% Downsample to 200 Hz
cfg = [];
cfg.resamplefs = 200;
dataBipolarResampled = ft_resampledata(cfg,dataBipolar);
dataBipolarScalpResampled = ft_resampledata(cfg,dataBipolarScalp);


%% Append again the data from scalp and macro
cfg = [];
data_bipolar_appended = ft_appenddata(cfg,dataBipolar,dataBipolarScalp);

%% Select correct Trials
[dataBipolar_SS,TrialInformationTableCorrect] = Get_Only_Correct_Trials_FieldTrip(data_bipolar_appended,TrialInformationTable);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTableCorrect);



%% Select the latencies for every task period for each task period
[dataBipolar_Ret_SS{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
[dataBipolar_Ret_SS{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
[dataBipolar_Ret_SS{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance



%% Channel pairs to run
%To be configured for every combination of electrode pairs you have
clear nChannelPairs;

nScalpChans = length(dataBipolarScalp.label);

nChannelPairs = [];

for i=1:length(Bip_chans)
    nChannelPairs = [nChannelPairs;[i+size(montage.labelold,1)+zeros(nScalpChans,1),(1:nScalpChans)'+Bip_chans(end) ]];
end
strChannelNameList = dataBipolar_Ret_SS{1}{1}.label;

hipp_chans= find(contains(strChannelNameList(nChannelPairs(:,1)),'H'));
hipp_chanPairs = nChannelPairs(hipp_chans,:);

%Validate that the channel pairs are correct
disp(strChannelNameList(hipp_chanPairs))


%% Granger
clear gdata PowCorr;
tStart = tic;
nSet_Size = size(Set_Sizes,1);
Fix_Data = dataBipolar_Ret_SS{1};
Enc_Data = dataBipolar_Ret_SS{2};
Maint_Data = dataBipolar_Ret_SS{3};
for iPair = 1:size(hipp_chanPairs,1)
    fprintf('Calculating Granger for pair %d for all set sizes\n',iPair)
    for iSS = 4:nSet_Size%1:nSet_Size
        freq                   = [];
        freq.freqcfg           = [];
        freq.freqcfg.method    = 'mtmfft';
        freq.freqcfg.foi       = [1:1:100];%[4:1:100];
        freq.freqcfg.output    = 'fourier';
        freq.freqcfg.tapsmofrq = 2;
        Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS});
        Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS});
        Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS});
        
        grangercfg = [];
        grangercfg.method  = 'granger';
        grangercfg.granger.conditional = 'no';
        grangercfg.granger.sfmethod = 'bivariate';
        grngChan1 = strChannelNameList{hipp_chanPairs(iPair,1)};
        grngChan2 = strChannelNameList{hipp_chanPairs(iPair,2)};
        grangercfg.channelcmb = {grngChan1,grngChan2};
        gdata.Maint{iPair,iSS}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
        gdata.Enc{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Encoding_freq);
        gdata.Fix{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Fixation_freq);
        
        
    end
    tStop = toc(tStart);
    fprintf('\n\n\n\n\n\n\n\n\n\n Pair %d %f seconds elapsed \n',iPair,tStop)
end
%% Save results
strSavePath = 'F:/Vasileios/Task Analysis/Data/Analysis Data/Granger for EDF dataset/';
strPatSavePath = [strSavePath,'Patient ', patient,'\'];
mkdir(strPatSavePath)
save([strPatSavePath,'Granger_depth_scalp_all_channel_pairs_all_set_sizes'],'gdata','-v7.3');

%% Visualize the Granger for one channel Pair

hipp_chan_to_plot = Bip_chans(12);% find(contains(strChannelNameList(Bip_chans),'H'))
freq_ax    =  gdata.Fix{1,4}.freq%gdata.Fix{1,1}.freq;
iSS = 4;
light_blue =  [0.30,0.75,0.93];
light_red  =  [1 0.45 0.45];

Colors = {'b',light_blue,'r',light_red};
chanPair = find(hipp_chanPairs==hipp_chan_to_plot);
figure;
semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
xlim([4 30])


%% Visualize the results (granger spectra in scalp eeg arrangement)
%Select Bipolar channel to plot
hipp_chan_to_plot = find(contains(strChannelNameList(Bip_chans),'H'));
Bip_chans_to_plot = Bip_chans(12)%hipp_chan_to_plot);
%Select set size to plot
iSS = 4;
% colors
light_blue =  [0.30,0.75,0.93];
light_red  =  [1 0.45 0.45];

Colors = {'b',light_blue,'r',light_red};

freq_ax    =  gdata.Fix{1,4}.freq%gdata.Fix{1,1}.freq;

for i = 1:numel(Bip_chans_to_plot)
    indBipChan = find(hipp_chanPairs==Bip_chans_to_plot(i))
    figure('units','normalized','outerposition',[0 0 1 1])
    ha = tight_subplot(5,9,[.07 .03],[.1 .01],[.01 .01])
    set(gcf,'color','white')
    figLayout = Get_Figure_Scalp_Layout_Information();
    for j = 1:numel(indBipChan)
        chanPair = indBipChan(j);
        Scalp_channel = strChannelNameList{hipp_chanPairs(chanPair,2)};
        indSubplot = find(strcmp(figLayout(:,1),Scalp_channel));
        if ~isempty(indSubplot)
            nSubplot(j) = figLayout{indSubplot,2};
            axes(ha(nSubplot(j)));
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',3);
            hold on;
            semilogx(freq_ax,gdata.Enc{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',3);
            semilogx(freq_ax,gdata.Maint{chanPair,iSS}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',3);
            xlim([4 30])
            ylabel(Scalp_channel)
            set(gca,'box','off')
            %             ylim([0 0.2])
        end
    end
    for k = 1:size(ha,1)
        
        if ~ismember(k,nSubplot)
            set(ha(k),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    suptitle(sprintf('Depth Electrode %s - Scalp Granger',dataBipolar_Ret_SS{1}{iSS}.label{Bip_chans_to_plot(i)}))
end