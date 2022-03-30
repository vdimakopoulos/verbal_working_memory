%%
clc;
clear;
close all;

%% Paths
Drive_Letter = 'E:\';
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


%% Load data for 2 sessions - Macro Data
strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\14 UC\'];
cd (strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
    data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end

%% Merge sessions
dataBipolar = data_ses(1).data;
for nSes = 2:length(data_ses)
    cfg = [];
    dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
end


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
cfg.refchannel    =  [41 42] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);


%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 9;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [33 41 35 36 36 41 37 38 38 39 39 40 40 41 41 42 42 34 ] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
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
cfg.refchannel = [43 51]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;


%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);

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


%% Run once

strChannelNameList = dataBipolar.label;
indGL = find(~cellfun(@isempty,strfind(strChannelNameList()','GL')));
indTL = find(~cellfun(@isempty,strfind(strChannelNameList()','TL')));
indTL = indTL([1:size(indTL,1)-num_bipolar_chans]); %Exclude Bipolar Channels

for i =1:length(indGL)
    
    nGridChanNum =i;
    strChanName = strcat('GL',int2str(nGridChanNum));
    gridGL_plot_order(i) = min(find(~cellfun(@isempty,strfind(strChannelNameList()', strChanName))));
    
end




for i =1:length(indTL)
    
    nGridChanNum =i;
    strChanName = strcat('TL',int2str(nGridChanNum));
    depthTL_plot_order(i) = min(find(~cellfun(@isempty,strfind(strChannelNameList()', strChanName))));
    
end



nGridChans = length(find(contains(dataBipolar.label,'GL')));
% k = 0;
strGridRows = {'GL A','GL B','GL C','GL D'};
for i = 1:nGridChans
    
    nGridNum{i} = str2double(strsplit(dataBipolar.label{i},'GL'))
    
    nChanNum(i) = nGridNum{i}(2);
    nRowLetter(i) = ceil(nChanNum(i)/8);
    dataBipolar.label{i}  = strcat(strGridRows{nRowLetter(i)},int2str((nChanNum(i)-(nRowLetter(i)-1)*8)));
    
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
        
        cfg =[];
        cfg.resamplefs = 2000;
        dataBipolar_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Ret_SS{iSS})
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:200;
    %     cfg.pad='nextpow2'
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_fix =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'encod') % Encoding
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_enc =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    
    
elseif strcmp(TaskPeriod,'maint') %Maintenance
    %% Extract 2 seconds of maintenance latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        
    end
    
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
    
else %Retrieval
    %% Extract 2 seconds of retrieval latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [0,2-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        
    end
end


%% Power Spectra of grid channels

%Grid GL
strPaths.Figures = [strPaths.Results,'14 UC Results\PSD\Grid Chans\'];
mkdir(strPaths.Figures)


fig = figure;
ha = tight_subplot(8,4,[.03 .03],[.03 .03],[.03 .03])
order = [29:-4:1 30:-4:2 31:-4:3 32:-4:4]
for i = 1:length(indGL)
    k = order(i);
    
    %    axes(ha(k));
    semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(gridGL_plot_order(i),:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(gridGL_plot_order(i),:)),'g','LineWidth',3);
    semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(gridGL_plot_order(i),:)),'r','LineWidth',3);
    xlim([4 40]);
    ylim([-20 30])
    ylabel(dataBipolar.label{gridGL_plot_order(i)});
end
suptitle('PSD of grid GL channels');
saveas(fig,[strPaths.Figures,'PSD of grid GL Chans'],'fig')
saveas(fig,[strPaths.Figures,'PSD of grid GL Chans'],'png')



%Depth TL
strPaths.Figures = [strPaths.Results,'14 UC Results\PSD\Depth Chans\'];
mkdir(strPaths.Figures)


fig = figure;
ha = tight_subplot(2,5,[.03 .03],[.03 .03],[.03 .03])
% order = [29:-4:1 30:-4:2 31:-4:3 32:-4:4]
for i = 1:num_bipolar_chans
    k = i;
    
    axes(ha(k));
    semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(i+size(montage.labelold,1),:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(i+size(montage.labelold,1),:)),'g','LineWidth',3);
    semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(i+size(montage.labelold,1),:)),'r','LineWidth',3);
    xlim([4 40]);
    ylim([-30 30])
    ylabel(dataBipolar.label{i+size(montage.labelold,1)});
end
suptitle('PSD of depth TL channels');
saveas(fig,[strPaths.Figures,'PSD of depth TL Chans'],'fig')
saveas(fig,[strPaths.Figures,'PSD of depth TL Chans'],'png')



%% TFR of Grid

%Grid GL
strPaths.Figures = [strPaths.Results,'14 UC Results\PSD\Grid Chans\'];
mkdir(strPaths.Figures)


Chan_to_plot = 1:nGridChans;
clim = [-1 1];
subplot_flag = 0;

cfg = [];
cfg.baseline = [TFR2.time(1) TFR2.time(5)]%[TFR2.time(8) TFR2.time(9)];
cfg.baselinetype = 'absolute';
TFR_baselined = ft_freqbaseline(cfg,TFR2)


fig = figure;
ha = tight_subplot(8,4,[.03 .03],[.03 .03],[.03 .03])
set(gcf,'color','white')
order = [29:-4:1 30:-4:2 31:-4:3 32:-4:4]

for i = Chan_to_plot(1):Chan_to_plot(end)
    k = order(i); %Plot order
    indGrid = gridGL_plot_order(i);
    
    axes(ha(k));
    
    TFR_psd_base = squeeze(TFR_baselined.powspctrm(indGrid,:,:))/3%
    contourf(TFR_psd.time,TFR_psd.freq,TFR_psd_base,100,'LineColor','none')
    
    
    %Axes properties
    set(gca,'clim',clim,'yscale','log')
    set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
    set(gca,'FontSize',6)
    set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
    colormap jet%(flipud(jet))
    colorbar
    ylabel(dataBipolar.label(indGrid));
    
end
suptitle('TFR PSD of grid GL channels');
saveas(fig,[strPaths.Figures,'TFR PSD of grid GL Chans'],'fig')
saveas(fig,[strPaths.Figures,'TFR PSD of grid GL Chans'],'png')


%% Channel pairs to run
% Hippocampus strip
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;

nGridChansGL = length(find(contains(dataBipolar_SS{4}.label,'GL')));
indTL_start = size(montage.labelold,1)+1;
nChannelPairs = [];%[[ones(nStripChans,1),(1:nStripChans)'+AHL_chans];[nStripChans+AHL_chans+ones(nStripChans,1),(1:nStripChans)'+AHL_chans]];

for i=1:2:length(chans_to_be_bipolar)
    nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChansGL,1),([1:nGridChansGL])']];
end


%% PLV Analysis
tStart = tic;
for iPair = 1:size(nChannelPairs,1)
    %% Channel numbers
    nChannel_1 = nChannelPairs(iPair,1);
    nChannel_2 = nChannelPairs(iPair,2);
    
    %% Create data structure with the selected channels for each set size
    data_SinglePair_SS = cell(1,5);
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.channel = [nChannel_1,nChannel_2];
        data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});%ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});
    end
    
    %% Phase coherence
    
    for iSS = 1:nSet_Size
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
    for iSS = 1:nSet_Size
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
                save([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_retrieval.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
                
            end
        end
    end
end


%% Plot the PLV results

iSS_to_Plot = 4;
k = 0;
strPaths.Figures = [strPaths.Results,'14 UC Results\PLV\Grid Depth\'];
mkdir(strPaths.Figures);
% Load the PLV results
if ~GammaFreqPLV
    
    PLV_Vars = load([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'])
    PLV_Pairs_maint = PLV_Vars.PLV_Pairs;
    freqAxis = PLV_Vars.freqAxis;
    
    PLV_Vars = load([strVariableFolder,'Patient_',num2str(14,'%.2d '),'_','Grid_Depth','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'])
    PLV_Pairs_fix = PLV_Vars.PLV_Pairs;
    nChannelPairs = PLV_Vars.nChannelPairs;
    freqAxis_fix = PLV_Vars.freqAxis;
end

%Grid PR
for i = 1:size(nChannelPairs,1)
    if nChannelPairs(i,2)<= 32
        if i == 1 || mod(i-1,length(indGL)) == 0
            fig = figure;
            ha = tight_subplot(8,4,[.03 .03],[.1 .2],[.03 .03])
            strTitle =['PLV of bipolar channel ',dataBipolar.label(nChannelPairs(i,1)),' and grid GL, ss ',num2str(Set_Sizes{iSS_to_Plot})]
            suptitle(strTitle);
            if i~=1
                k = floor(i/length(indGL));
            end
        end
        ax = order(nChannelPairs(i,2)); %Plot order
        axes(ha(ax));
        ind = (gridGL_plot_order(nChannelPairs(i,2))+k*nGridChans);
        semilogx(freqAxis_fix,abs(PLV_Pairs_fix{ind,iSS_to_Plot}),'k','LineWidth',3);
        hold on;
        semilogx(freqAxis,abs(PLV_Pairs_maint{ind,iSS_to_Plot}),'r','LineWidth',3);
        xlim([4 100]);
        ylim([0 1]);
        ylabel(dataBipolar.label(nChannelPairs(ind,2)));
        
    end
    
end

for i = 1:fig.Number
    figr = figure(i);
    saveas(gcf,char(strcat(strPaths.Figures,'PLV of bipolar channel_',dataBipolar.label(Bip_chans(i)),'_and grid GL, ss ',num2str(Set_Sizes{iSS_to_Plot}))),'fig')
    saveas(gcf,char(strcat(strPaths.Figures,'PLV of bipolar channel_',dataBipolar.label(Bip_chans(i)),'_and grid GL, ss ',num2str(Set_Sizes{iSS_to_Plot}))),'png')
    
end

