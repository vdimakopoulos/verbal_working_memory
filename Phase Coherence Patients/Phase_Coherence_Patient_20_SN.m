%% Close all figures, clear variables and command window
close all
clear
clc

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
strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\20 SN\'];
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
cfg.refchannel    =  [48] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);


%% Downsample to 500 Hz
cfg =[];
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

%% Run Once

strChannelNameList = dataBipolar.label;
indPL = find(~cellfun(@isempty,strfind(strChannelNameList()','PL')));
indTL = find(~cellfun(@isempty,strfind(strChannelNameList()','TL')));

for i =1:length(indPL)
    
   nGridChanNum =i;
   strChanName = strcat('PL',int2str(nGridChanNum));
   gridPL_plot_order(i) = min(find(~cellfun(@isempty,strfind(strChannelNameList()', strChanName))));
    
end


for i =1:length(indTL)
    
   nGridChanNum =i;
   strChanName = strcat('TL',int2str(nGridChanNum));
   gridTL_plot_order(i) = min(find(~cellfun(@isempty,strfind(strChannelNameList()', strChanName))));
    
end


nGridChans = length(find(contains(dataBipolar.label,'PL')));
% k = 0;
strGridRows = {'A','B','C','D'};
for i = 1:nGridChans
    
    nGridNum{i} = str2double(strsplit(dataBipolar.label{i},'PL'))
    
    nChanNum(i) = nGridNum{i}(2);
    nRowLetter(i) = ceil(nChanNum(i)/8);
    dataBipolar.label{i}  = strcat(strGridRows{nRowLetter(i)},int2str((nChanNum(i)-(nRowLetter(i)-1)*8)));
    
end

nGridChansTL = length(find(contains(dataBipolar.label,'TL')));
% k = 0;
strGridRows = {'A','B'};
for i = nGridChans+1:nGridChans+nGridChansTL
    
    nGridNum{i} = str2double(strsplit(dataBipolar.label{i},'TL'))
    
    nChanNum2(i) = nGridNum{i}(2);
    nRowLetter2(i) = ceil(nChanNum2(i)/8);
    dataBipolar.label{i}  = strcat(strGridRows{nRowLetter2(i)},int2str((nChanNum2(i)-(nRowLetter2(i)-1)*8)));
    
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

%Grid PL
strPaths.Figures = [strPaths.Results,'20 SN Results\PSD\Grid Chans\'];
mkdir(strPaths.Figures)

fig = figure;
ha = tight_subplot(8,4,[.03 .03],[.1 .1],[.1 .1])
plot_order =[29:-4:1 30:-4:2 31:-4:3 32:-4:4]
for i = 1:length(indPL)
    axes(ha(plot_order(i)));
    semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(gridPL_plot_order(i),:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(gridPL_plot_order(i),:)),'g','LineWidth',3);
    semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(gridPL_plot_order(i),:)),'r','LineWidth',3);
    xlim([4 100]);
    ylim([-20 30])
    ylabel(dataBipolar.label{gridPL_plot_order(i)});
end
suptitle('PSD of grid PL channels');
saveas(fig,[strPaths.Figures,'PSD of grid PL Chans'],'fig')
saveas(fig,[strPaths.Figures,'PSD of grid PL Chans'],'png')


%Grid TL
fig = figure;
ha = tight_subplot(2,8,[.03 .01],[.1 .2],[.02 .02])
for i = 1:length(indTL)
    axes(ha(i));
    semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(gridTL_plot_order(i),:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(gridTL_plot_order(i),:)),'g','LineWidth',3);
    semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(gridTL_plot_order(i),:)),'r','LineWidth',3);
    xlim([4 40]);
    ylim([-10 30])
    ylabel(dataBipolar.label{gridTL_plot_order(i)});
end
suptitle('PSD of grid TL channels');
saveas(fig,[strPaths.Figures,'PSD of grid TL Chans'],'fig')
saveas(fig,[strPaths.Figures,'PSD of grid TL Chans'],'png')

%% TFR GRID
%Grid PL
strPaths.Figures = [strPaths.Results,'20 SN Results\TFR PSD\Grid Chans\'];
mkdir(strPaths.Figures)

Chan_to_plot = 1:32
clim = [-1 0.5];
subplot_flag = 0;

fig = figure;
ha = tight_subplot(8,4,[.03 .03],[.1 .1],[.1 .1])

cfg = [];
cfg.baseline = [TFR2.time(1) TFR2.time(5)]%[TFR2.time(8) TFR2.time(9)];
cfg.baselinetype = 'absolute';
TFR_baselined = ft_freqbaseline(cfg,TFR2)

set(gcf,'color','white')

for i = Chan_to_plot(1):Chan_to_plot(end)
    indGrid = find(gridPL_plot_order == i)
    
    axes(ha(plot_order(indGrid)))
    TFR_psd_base = squeeze(TFR_baselined.powspctrm(i,:,:))/5-1%/3-1%
    contourf(TFR_psd.time,TFR_psd.freq,TFR_psd_base,100,'LineColor','none')
    
     
    %Axes properties
    set(gca,'clim',clim,'yscale','log')
    set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
    set(gca,'FontSize',6)
    set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
    colormap jet%(flipud(jet))
    colorbar
    ylabel(dataBipolar.label(i));
    
end
suptitle('TFR PSD of grid PL channels');
saveas(fig,[strPaths.Figures,'TFR PSD of grid PL Chans'],'fig')
saveas(fig,[strPaths.Figures,'TFR PSD of grid PL Chans'],'png')


%GridTL
Chan_to_plot = 33:length(dataBipolar.label)
clim = [-1 0];
subplot_flag = 0;

fig = figure;
ha = tight_subplot(2,8,[.03 .01],[.1 .2],[.02 .02])

cfg = [];
cfg.baseline = [TFR2.time(1) TFR2.time(5)]%[TFR2.time(8) TFR2.time(9)];
cfg.baselinetype = 'absolute';
TFR_baselined = ft_freqbaseline(cfg,TFR2)

set(gcf,'color','white')

for i = Chan_to_plot(1):Chan_to_plot(end)
    indGrid = find(gridTL_plot_order == i)
    
    axes(ha(indGrid))
    TFR_psd_base = squeeze(TFR_baselined.powspctrm(i,:,:))/5-1%/3-1%
    contourf(TFR_psd.time,TFR_psd.freq,TFR_psd_base,100,'LineColor','none')
    
     
    %Axes properties
    set(gca,'clim',clim,'yscale','log')
    set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
    set(gca,'FontSize',6)
    set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
    colormap jet%(flipud(jet))
    colorbar
    ylabel(dataBipolar.label(i));
    
end
suptitle('TFR PSD of grid TL channels');
saveas(fig,[strPaths.Figures,'TFR PSD of grid TL Chans'],'fig')
saveas(fig,[strPaths.Figures,'TFR PSD of grid TL Chans'],'png')


%% Channel pairs to run
%Grid PL - Grid TL
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;

nGridChans = length(indPL);
nChannelPairs = [];

for i=nGridChans+1:nGridChans+length(indTL)
    nChannelPairs = [nChannelPairs;[i+zeros(nGridChans,1),(1:nGridChans)']];
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
                save([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_retrieval.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
                
            end
        end
    end
end


%% Plot the PLV results

iSS_to_Plot = 4;
k = 0;

% Load the PLV results
if ~GammaFreqPLV
    
    PLV_Vars = load([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'])
    PLV_Pairs_maint = PLV_Vars.PLV_Pairs;
    freqAxis = PLV_Vars.freqAxis;

    
%     PLV_Vars = load([strVariableFolder,'Patient_',num2str(12,'%.2d '),'_','Grid_Depth','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'])
%     PLV_Pairs_encod = PLV_Vars.PLV_Pairs;
%     freqAxis = PLV_Vars.freqAxis;
    
    PLV_Vars = load([strVariableFolder,'Patient_',num2str(20,'%.2d '),'_','Grid_Depth','_PLV_Pairs_fixation_f_1_1_30_5_100.mat'])
    PLV_Pairs_fix = PLV_Vars.PLV_Pairs;
    nChannelPairs = PLV_Vars.nChannelPairs;
    freqAxis_fix = PLV_Vars.freqAxis;
end



for i = 129:160%1:size(nChannelPairs,1)
    if i == 1 || mod(i-1,length(indPL)) == 0
        fig = figure;
        ha = tight_subplot(8,4,[.03 .03],[.1 .1],[.1 .1])
        strTitle =['PLV of grid TL ',dataBipolar.label(nChannelPairs(i,1)),' and grid PL, ss ',num2str(Set_Sizes{iSS_to_Plot})]
        suptitle(strTitle);
        
        if i~=1
            k = floor(i/length(indPL));
        end
    end
    ind  = find(gridPL_plot_order == nChannelPairs(i,2))
    axes(ha(plot_order(ind)));
    semilogx(freqAxis_fix,abs(PLV_Pairs_fix{i,iSS_to_Plot}),'k','LineWidth',3);
    hold on;
    %     semilogx(freqAxis,abs(PLV_Pairs_encod{i,iSS_to_Plot}),'g','LineWidth',3);
    semilogx(freqAxis,abs(PLV_Pairs_maint{i,iSS_to_Plot}),'r','LineWidth',3);
    xlim([4 100]);
    ylim([0 1]);
    ylabel(dataBipolar.label(nChannelPairs(i,2)));
    
    
end