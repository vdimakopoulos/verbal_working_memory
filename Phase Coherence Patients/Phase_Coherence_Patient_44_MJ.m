%% Close all figures, clear variables and command window
close all
clear
clc

%% Paths
Drive_Letter = 'F:\';
strPaths.Main = [Drive_Letter,'Vasileios\'];
strPaths.Project = [strPaths.Main, 'Task Analysis\'];
strPaths.GeneralFunctions = [Drive_Letter,'Vasileios\Task Analysis\Code\'];%'F:\Vasileios\Task Analysis\Code\';
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
strPaths.Toolboxes.FieldTrip            = [Drive_Letter,'Vasileios\Toolboxes\fieldtrip-20200315\'];

% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = [Drive_Letter,'Vasileios\Toolboxes\eeglab14_1_1b\'];

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

%% Load data for 4 sessions - Macro Data
strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\44 MJ\'];
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end


%% Load data for 4 sessions - Scalp Data
strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\44 MJ\'];
cd (strScalpDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_Scalp_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable_Scalp = [];
for i=1:length(data_Scalp_ses)
    TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
end


%% Merge sessions
dataBipolar = data_ses(1).data;
for nSes = 2:length(data_ses)
    cfg = [];
    dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
    
end

dataBipolarScalp = data_Scalp_ses(1).data;
for nSes = 2:length(data_Scalp_ses)
    cfg = [];
    dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
    
end
cd(strPaths.Project)
%% Flags and Parameters
TaskPeriod = 'maint'; %1 'encod' for encoding, 'maint' for maintenance, 'retr' for retrieval,'fix' for fixation
% TaskPeriod = 'encod';
% TaskPeriod = 'retr';
Scalp_Hippocampus_coupling = 0;
Grid_Hippocampus_coupling = 1;
Scalp_Grid_coupling = 0;
GammaFreqPLV = 0;
EEG_Preproc = 1;
Visual_Evoked_Potential = 0;
PSD_calc = 0;


if Scalp_Hippocampus_coupling
    coupling_type = 'Scalp-Hippocampus';
    strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
elseif Grid_Hippocampus_coupling
    coupling_type = 'Grid-Hippocampus';
    strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\'];
elseif Scalp_Grid_coupling
    coupling_type = 'Scalp-Grid';
    strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];
end

%% Rereferencing
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =  [3] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);


%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 24;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 33 34 33 35 34 35 41 42 41 43 42 43 ...
    49 50 49 51 50 51 57 58 57 59 58 59] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
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
cfg.refchannel = [66:1:88]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);

%% macro data
% Select only correct trials
[dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);


%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable);
%%
if EEG_Preproc
    
    
    %
    %% Resample
    cfg = [];
    cfg.resamplefs = 500;
    dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);
    
    %% Visualize scalp raw data
    cfg =[];
    cfg.viewmode = 'vertical';
    cfg.preproc.lpfilter = 'yes'
    cfg.preproc.lpfreq = 45;
    cfg.trial = 5;
    ft_databrowser(cfg,dataBipolarScalp);
    %% ICA
    cfg = [];
    cfg.method       = 'runica'
    cfg.channel      = {'all' '-Submm','-Submp','-EKGm','-EKGp'};
    cfg.numcomponent = 'all';
    cfg.demean       = 'yes';
    cfg.feedback     =  'text'
    IC_components = ft_componentanalysis(cfg,dataBipolarScalp)
    
    %% Plot the IC components of ICA using topoplot
    figure;
    cfg = [];
    cfg.component  = [1:size(IC_components.topo,1)];
    cfg.layout      = 'elec1020.lay';
    cfg.colormap    = 'jet';
    cfg.colorbar    = 'no'
    cfg.comment = 'no';
    ft_topoplotIC(cfg,IC_components);
    
    %%
    
    cfg = [];
    cfg.layout = 'elec1020.lay';
    cfg.viewmode = 'component';
    
    ft_databrowser(cfg, IC_components);
    
    %% Reject IC componentss
    cfg = [];
    cfg.component = [1 2 5 6 9 12 13 14 15 16 18 20 21];
    cfg.demean = 'yes';
    dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);
    
    cfg =[];
    cfg.channel = {'all' '-Submm','-Submp','-EKGm','-EKGp'};
    dataBipolarScalp_Reconstructed = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    
    
    %% re-ref EEG
    cfg               = []
    cfg.channel       = {'all'}
    cfg.reref         =  'yes'
    cfg.refchannel    =   [1, 2]
    cfg.refmethod     =  'avg'
    dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp_Reconstructed);
    
    %%
    cfg =[];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg,dataBipolarScalp);
    
    
    %% Select only correct trials
    [dataBipolar_Scalp_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable);
    
    %% Divide data into set sizes
    %   [4] -> Set Size 1
    %       [6] -> Set Size 2
    %           [8] -> Set Size 3
    %               [6 8] -> Set Size 4
    %                   [4 6 8] -> Set Size 5
    Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
    [dataBipolar_Scalp_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_Scalp_SS,TrialInformationTable);
    
    
end


%% Select the latencies for every task period and calculate the PSD for each task period
if strcmp(TaskPeriod,'fix') %Fixation
    %% Extract 2 seconds of fixation latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
        cfg =[];
        cfg.resamplefs = 1000;
        dataBipolar_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Ret_SS{iSS})
        dataBipolar_Scalp_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Scalp_Ret_SS{iSS})
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_fix =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    fr_fix_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'encod') % Encoding
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_enc =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    fr_enc_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'maint') %Maintenance
    %% Extract 2 seconds of maintenance latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    fr_maint_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
else %Retrieval
    %% Extract 2 seconds of retrieval latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [0,2-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        %         dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
    end
end

%% Power Spectra of Depth electrodes
% Select 1 or more bipolar chans
Bip_Chan_to_Plot = Bip_chans(3)%1:end);

for i = 1:length(Bip_Chan_to_Plot)
    fig = figure(i);
    semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(Bip_Chan_to_Plot(i),:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(Bip_Chan_to_Plot(i),:)),'g','LineWidth',3);
    semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(Bip_Chan_to_Plot(i),:)),'r','LineWidth',3);
    xlim([4 100]);
    ylim([-20 40])
    strTitle = ['PSD of Bipolar Channel ',dataBipolar.label{Bip_Chan_to_Plot(i)}];
    title(strTitle)
    hold off;
end



%% Scalp PLV with Depth Electrodes during fixation
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_fix_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P44(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_fix_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else%In case you forgot to declare the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_fix_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end


%% Scalp PLV with Depth Electrodes during encoding
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_enc_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P44(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end


%% Scalp PLV with Depth Electrodes during maintenance
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_maint_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P44(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_maint_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_maint_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end


%% Load the PLV vars
strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
load([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat']);
load([strVariableFolder,'Patient_',num2str(44,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_f_1_1_30_5_100.mat']);


%% Plot Scalp Depth PLV for one or more bipolar channels
Bipolar_Channels = unique(nChannelPairs_Sclp_depth(:,1));
Bip_chans_to_plot = Bipolar_Channels(23)%1:end);
iSS = 4; % Set Size [6 8]
% Scalp topoplot related
FreqBand_Sclp = [6,7];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

for i = 1:numel(Bip_chans_to_plot)
    nFig = figure
    indBipChan = find(nChannelPairs_Sclp_depth==Bip_chans_to_plot(i))
    ha =  tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
    set(gcf,'color','white')
    figLayout = Get_Figure_Scalp_Layout_Information();
    for j = 1:numel(indBipChan)
        Scalp_channel = strcat({' '},dataBipolar_appended{iSS}.label{j});
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot(j) = figLayout{indSubplot,2};
        axes(ha(nSubplot(j)));
        datavector{i}(j) = median(abs(PLV_Pairs_Scalp_maint_depth{indBipChan(j),iSS}(indFreq1:indFreq2)));
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_depth{indBipChan(j),iSS}),'k','LineWidth',3);
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_depth{indBipChan(j),iSS}),'r','LineWidth',3);
        ylim([0 1]);
        xlim([4 100]);
        ylabel(Scalp_channel);
        hold off;
    end
    suptitle(sprintf('Depth Electrode %s - Scalp PLV',dataBipolar_appended{iSS}.label{Bip_chans_to_plot(i)}))
    for k = 1:size(ha,1)
        
        if ~ismember(k,nSubplot)
            set(ha(k),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    nFig = figure
    Scalp_topoplot(strPaths.Toolboxes,FreqBand_Sclp,freqAxis_sclp,datavector{i},dataBipolarScalp,4,nChannelPairs_Sclp_depth,1)
end
