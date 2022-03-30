%%
clc;
clear;
close all;

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


%% Load data for 7 sessions - Macro Data
strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\49 TK\'];
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(files)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end


%% Load data for 7 sessions - Scalp Data
strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\49 TK\'];
cd (strScalpDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_Scalp_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable_Scalp = [];
for i=1:length(files)
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
cfg.refchannel    =  [3]% 12 22 29 37] %add white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);



%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 9;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [1 3 2 3 9 13 10 13 11 13 12 13 17 18 49 51 50 51]%Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc. need wm contacts
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
cfg.refchannel = [65 73] % the number of the bipolar channels
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

%% Select the latencies for every task period and calculate the PSD for each task period
if strcmp(TaskPeriod,'fix') %Fixation
    %% Extract 2 seconds of fixation latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        %         dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
        cfg =[];
        cfg.resamplefs = 1000;
        dataBipolar_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Ret_SS{iSS})
        %         dataBipolar_Scalp_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Scalp_Ret_SS{iSS})
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_fix =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    %     fr_fix_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'encod') % Encoding
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        %         dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_enc =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    %     fr_enc_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'maint') %Maintenance
    %% Extract 2 seconds of maintenance latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        %         dataBipolar_Scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_Scalp_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    %     fr_maint_scalp =ft_freqanalysis(cfg,dataBipolar_Scalp_Ret_SS{4});
    
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


%% Channel pairs to run
% Hippocampus depth
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;
MTL_Bip_chans = chans_to_be_bipolar(1:end);

nDepth = size(montage.labelold,1);
nChannelPairs = [];
cortical_chans = [25:48 57:64];

for i=1:2:length(MTL_Bip_chans)
%     nChannelPairs = [nChannelPairs;[round(i/2)+nDepth+zeros(length(Bip_chans(end-2):Bip_chans(end)),1),([Bip_chans(end-2):Bip_chans(end)])']];
    nChannelPairs = [nChannelPairs;[round(i/2)+nDepth+zeros(length(cortical_chans),1),cortical_chans']];

end
strChannelNameList = dataBipolar.label;


%% PLV Analysis
tStart = tic;
for iPair = 1:size(nChannelPairs,1)
    %% Channel numbers
    nChannel_1 = nChannelPairs(iPair,1);
    nChannel_2 = nChannelPairs(iPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,5);
    for iSS = 4:nSet_Size
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});%ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
    end
    
    %% Phase coherence
    
    for iSS = 4:nSet_Size
        %% Frequency analysis
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
        %         cfg.foi         = 0.5:0.5:30;
        if GammaFreqPLV
            cfg.foi         = 0.5:10:120;
        else
            cfg.foi         = [1:1:29, 30:5:100];
            
        end
        cfg.tapsmofrq   = 2;
        cfg.channel     = 1:2;
        Freq            = ft_freqanalysis(cfg,data_SinglePair_SS{iSS});
        
        cfgFreq = cfg;
        
        %% PLV analysis
        cfg             = [];
        cfg.method      = 'plv';
        cfg.complex     = 'complex';
        cfg.feedback    = 'none';
        PLV_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
        
        %% Coherence
        cfg             = [];
        cfg.method      = 'coh';
        cfg.complex     = 'complex';
        cfg.feedback    = 'none';
        Coh_SS{iSS} = ft_connectivityanalysis(cfg,Freq);
        
    end
    %% Store results for the channel pair
    strPairLabels = PLV_SS{iSS}.label';
    freqAxis = PLV_SS{iSS}.freq;
    for iSS = 4:nSet_Size
        PLV_Spectrum_SS{iSS} = squeeze(PLV_SS{iSS}.plvspctrm(1,2,:))';
        Coh_Spectrum_SS{iSS} = squeeze(Coh_SS{iSS}.cohspctrm(1,2,:))';
    end
    TimeInterval = [-2,-1/dataBipolar_SS{iSS}.fsample];
    
    %%
    PLV_Pairs(iPair,:) = PLV_Spectrum_SS;
    Coh_Pairs(iPair,:) = Coh_Spectrum_SS;
    
    tStop = toc(tStart);
    fprintf('\n\n\n\n\n\n\n\n\n\n Pair %d %f seconds elapsed \n',iPair,tStop)
end
% Save PLV Pairs
if 1
    if Grid_Hippocampus_coupling
        strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\']
        mkdir(strVariableFolder);
        
        if ~GammaFreqPLV
            if strcmp(TaskPeriod,'maint')
                save([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_maintenance2_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_retrieval.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_fixation2_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
                
            end
        end
    end
end

%% Plot all the PLV results

PLV_all_periods{1} = load([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_fixation2_f_1_1_30_5_100.mat'])
freqAxis = PLV_all_periods{1}.freqAxis;
nChannelPairs = PLV_all_periods{1}.nChannelPairs;
PLV_all_periods{1} = PLV_all_periods{1}.PLV_Pairs;

% PLV_all_periods{2} = load([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'])
% PLV_all_periods{2} = PLV_all_periods{2}.PLV_Pairs;

PLV_all_periods{3} = load([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_PLV_Pairs_maintenance2_f_1_1_30_5_100.mat'])
PLV_all_periods{3} = PLV_all_periods{3}.PLV_Pairs;

iSS = [4 5];
strSaveFigpath = 'F:/Vasileios/Task Analysis/Analysis Results/49 TK/Figures/PLV/';
mkdir(strSaveFigpath);
for i = 1:length(nChannelPairs)
    figure;
    for j = 1:length(iSS)
        subplot(2,1,j)
        semilogx(freqAxis,abs(PLV_all_periods{1}{i,iSS(j)}),'k','LineWidth',3)
        hold on;
%         semilogx(freqAxis,abs(PLV_all_periods{2}{i,iSS(j)}),'g','LineWidth',3)
        semilogx(freqAxis,abs(PLV_all_periods{3}{i,iSS(j)}),'r','LineWidth',3)
        xlim([4 100])
        set(gca,'box','off','XTick',[4 10 20 30 100],'XTickLabel',[4 10 20 30 100],'FontSize',12)
        xlabel('Frequency (Hz)')
        ylabel('PLV')
    end
    ttl = ['Channel Pair ',sprintf('%s--%s',strChannelNameList{nChannelPairs(i,1)},strChannelNameList{nChannelPairs(i,2)})]
    suptitle(ttl)
    figName = [strSaveFigpath, sprintf('Patient49_PLV_pair %s--%s',strChannelNameList{nChannelPairs(i,1)},strChannelNameList{nChannelPairs(i,2)})];
    saveas(gcf,figName);
    saveas(gcf,figName,'png')
end



%% Granger
nSet_Size = size(dataBipolar_SS,2);
Maint_Data = {};
iSS_current = 4; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]

for iSS =iSS_current% 1:nSet_Size
    cfg = [];
    cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
    Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

nSet_Size = size(dataBipolar_SS,2);
Enc_Data = {};
for iSS = iSS_current%1:nSet_Size
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
    Enc_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

nSet_Size = size(dataBipolar_SS,2);
Fix_Data = {};
for iSS = iSS_current%1:nSet_Size
    cfg = [];
    cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample];
    Fix_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end
iSS = 4;
freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.foi       = [1:1:100];%[4:1:100];
freq.freqcfg.output    = 'fourier';
freq.freqcfg.tapsmofrq = 2;
Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS});
Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS});
Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS});





gdata_Pairs = [];

for i = 1:length(nChannelPairs)
    
    % Channel numbers
    nChannel_1 = nChannelPairs(i,1);
    nChannel_2 = nChannelPairs(i,2);
    
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    grangercfg.channelcmb = {char(strChannelNameList(nChannel_1)), char(strChannelNameList(nChannel_2))};%{'_AHL2-_AHL3','_GL_C2'}; %{'_GL_C2','_AHL2-_AHL3'}
    grangercfg.granger.conditional = 'no';
    % grangercfg.complex = 'yes';
    grangercfg.granger.sfmethod = 'bivariate';
    
    gdata_Pairs.Maint{i}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata_Pairs.Enc{i}      = ft_connectivityanalysis(grangercfg, Encoding_freq);
    gdata_Pairs.Fix{i}      = ft_connectivityanalysis(grangercfg, Fixation_freq);
    
end
strVariableFolder = 'F:/Vasileios/Task Analysis/Analysis Results/49 TK/Analysis Data/Granger/';
mkdir(strVariableFolder)
save([strVariableFolder,'Patient_',num2str(49,'%.2d '),'_','Depth_Cortex','_Granger_f_1_1_30_5_100_','setSize_',num2str(iSS_current,'%.2d'),'.mat'],'gdata_Pairs','-v7.3')



%% Visualize Granger for all pairs

light_blue = [0.30,0.75,0.93];
light_red       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',light_red};
freq_ax    = gdata_Pairs.Fix{1}.freq;

strSaveFigpath = 'F:/Vasileios/Task Analysis/Analysis Results/49 TK/Figures/Granger/Spectra/';
mkdir(strSaveFigpath)
for i = 1:length(nChannelPairs)
    figure
    
    semilogx(freq_ax,gdata_Pairs.Enc{i}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',2);
    hold on;
    semilogx(freq_ax,gdata_Pairs.Enc{i}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',2);
    semilogx(freq_ax,gdata_Pairs.Maint{i}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',2);
    semilogx(freq_ax,gdata_Pairs.Maint{i}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',2);
    
    xlim([4 100]);
    set(gca,'box','off','XTick',[4 10 20 30 100],'XTickLabel',[4 10 20 30 100],'FontSize',12)
    xlabel('Frequency (Hz)')
    ylabel('Granger')
    ttl = ['Channel Pair ',sprintf('%s--%s',strChannelNameList{nChannelPairs(i,1)},strChannelNameList{nChannelPairs(i,2)})]
    suptitle(ttl)
    figName = [strSaveFigpath, sprintf('Patient49_PLV_pair %s--%s',strChannelNameList{nChannelPairs(i,1)},strChannelNameList{nChannelPairs(i,2)})];
    saveas(gcf,figName);
    saveas(gcf,figName,'png')
end


