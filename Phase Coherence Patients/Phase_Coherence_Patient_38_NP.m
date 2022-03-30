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

%% Load data for 6 sessions - Macro Data
strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\38 NP\'];
cd(strMacroDataDir);
files = dir('*.mat'); %try to look for all the .mat files under the folder
for i=1:length(files)
 data_ses(i) = load(files(i).name); %load the files
end

TrialInformationTable = []
for i=1:length(data_ses)
    TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
end


%% Load data for 6 sessions - Scalp Data
strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\38 NP\'];
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
for nSes = 1:length(data_ses)
    if nSes<=2
        cfg = [];
        cfg.channel = {'all','-mTLIL*','-mTLSL5','-mTLSL6'}
        data_ses(nSes).data = ft_selectdata(cfg,data_ses(nSes).data)
    end
    if nSes == 1
        dataBipolar = data_ses(nSes).data;
        
    else
        cfg = [];
        dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
    end
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

%% Reref to different sources
strChannelNameList = dataBipolar.label;
AhippChans = find(contains(strChannelNameList,'AH'));
PhippChans = find(contains(strChannelNameList,'PH'));
GridChans = find(contains(strChannelNameList,'T'));
dataRerefAhipp = reref_toSeparate_Chan(dataBipolar, 3, 'all','avg');
dataRerefGrid = reref_toSeparate_Chan(dataBipolar, 19, 'all','avg');
%selection of channels
cfg = [];
cfg.channel = [AhippChans; PhippChans];
dataRerefAhipp = ft_selectdata(cfg,dataRerefAhipp);
cfg = [];
cfg.channel = GridChans;
dataRerefGrid = ft_selectdata(cfg,dataRerefGrid);


% append them again
cfg = [];
dataBipolar = ft_appenddata(cfg,dataRerefAhipp,dataRerefGrid);
%% Rereferencing
cfg               =   [];
cfg.reref         =  'yes' 
cfg.refchannel    =  [19] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

dataBipolar.label = strrep(dataBipolar.label,'m','');

%% Apply montage %%
clear montage;
montage.labelold        = dataBipolar.label;
num_bipolar_chans       = 9;
num_reference_chans     = size(montage.labelold,1);

%prepare the montage matrix
montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
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
cfg.refchannel = [37:45]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
macro_data = dataBipolar;


%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 80;%500;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);


cfg = [];
cfg.channel = {'all' '-Subm1','-Subm2'};
dataBipolarScalp = ft_selectdata(cfg,dataBipolarScalp)
%% macro data
% Select only correct trials
[dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
% 
% cfg = [];
% cfg.trials = setdiff([1:length(dataBipolar_SS.trialinfo)],[177 178])
% dataBipolar_SS = ft_selectdata(cfg,dataBipolar_SS);
%Update trial information table
% TrialInformationTable = TrialInformationTable([1:176,179:size(TrialInformationTable,1)],:);

%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable);


%% Visual inspection
cfg = [];
cfg.viewmode = 'vertical'
ft_databrowser(cfg,dataBipolar_SS{4});

cfg          = [];
cfg.method   = 'trial';
cfg.keeptrials ='nan';
dummy        = ft_rejectvisual(cfg,dataBipolar_SS{4});

%% Reject trial
cfg =[];
cfg.trials = setdiff([1:length(dataBipolar_SS{4}.trialinfo)],[1 2:2:8 9 11:1:22 25 27 28 31 32 34:1:36 39:1:42 47 52 54:1:57 59 62 65 66 69:1:71 73 75 79 85 88 102 107 108 115 116 121 125 131 138 143 144])';
dataBipolar_SS{4} =  ft_selectdata(cfg,dataBipolar_SS{4});
%% 
if EEG_Preproc
%     
    %% Resample 
    cfg = [];
    cfg.resamplefs = 500;
    dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);
    
    %% Visualize scalp raw data
    cfg =[];
    cfg.viewmode = 'vertical';d
    ft_databrowser(cfg,dataBipolarScalp);
    
    %% Select only correct trials
    [dataBipolar_Scalp_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp,TrialInformationTable_Scalp);
    
    cfg = [];
    cfg.trials = setdiff([1:length(dataBipolar_Scalp_SS.trialinfo)],[177 178])
    dataBipolar_Scalp_SS = ft_selectdata(cfg,dataBipolar_Scalp_SS);
    %Update trial information table
    TrialInformationTable_Scalp = TrialInformationTable_Scalp([1:176,179:size(TrialInformationTable_Scalp,1)],:)

  %% ICA  
    cfg = [];
    cfg.method       = 'runica'
    cfg.channel      = {'all' '-_Submm','-_Submp'};
    cfg.numcomponent = 'all';
    cfg.demean       = 'yes';
    cfg.feedback     =  'text'
    IC_components = ft_componentanalysis(cfg,dataBipolar_Scalp_SS)
    
    
    %% IC power
    cfg              = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.tapsmofrq   = 2;
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    freq = ft_freqanalysis(cfg, IC_components);
    
    nsubplots = 25;
    nbyn = sqrt(nsubplots);% sqrt(nsubplots) should not contain decimals,
    
    Nfigs = ceil(size(IC_components.topo,1)/nsubplots);
    tot = Nfigs*nsubplots;
    
    rptvect = 1:size(IC_components.topo,1);
    rptvect = padarray(rptvect, [0 tot-size(IC_components.topo,1)], 0,'post');
    rptvect = reshape(rptvect,nsubplots,Nfigs)';
    
    for r=1:size(rptvect,1);
        figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        k=0;
        for j=1:size(rptvect,2);
            if~(rptvect(r,j)==0);
                k=k+1;
                cfg=[];
                cfg.channel = rptvect(r,j);
                subplot(nbyn,nbyn,k);ft_singleplotER(cfg,freq);
            end
        end
    end
    
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
    cfg.component = [1:7 10 12 13 17 20];% maybe add 18
    cfg.demean = 'yes';
    dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);
        
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
        
    
    %% Divide data into set sizes
    %   [4] -> Set Size 1
    %       [6] -> Set Size 2
    %           [8] -> Set Size 3
    %               [6 8] -> Set Size 4
    %                   [4 6 8] -> Set Size 5
  
    Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
    [dataBipolar_Scalp_SS,TrialInformationTableScalp_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_Scalp_SS,TrialInformationTable_Scalp);
    %% Reject trials
    cfg =[];
    cfg.trials = setdiff([1:length(dataBipolar_Scalp_SS{4}.trialinfo)],[1 2:2:8 9 11:1:22 25 27 28 31 32 34:1:36 39:1:42 47 52 54:1:57 59 62 65 66 69:1:71 73 75 79 85 88 102 107 108 115 116 121 125 131 138 143 144])';
    dataBipolar_Scalp_SS{4} =  ft_selectdata(cfg,dataBipolar_Scalp_SS{4});
    
   %%
    cfg =[];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg, dataBipolar_Scalp_SS{4});
    
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
        dataBipolar_Scalp_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Scalp_SS{iSS});

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
    cfg.foi         = [0 40]%0.5:0.5:100;
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

    end
end



%% Power Spectra of Bipolar Depth electrodes
% Select 1 or more bipolar chans
Bip_Chan_to_Plot = Bip_chans(1:end);

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


%% Channel pairs to run
% Hippocampus strip
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;

nStripChans = length(find(contains(dataBipolar.label,'T')));
PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
AHR_chans = length(find(contains(montage.labelold, 'AHR')==1));
PHR_chans = length(find(contains(montage.labelold, 'PHR')==1));
Depth_electrodes = PHL_chans+AHR_chans+PHR_chans;

nChannelPairs = [];%[[ones(nStripChans,1),(1:nStripChans)'+AHL_chans];[nStripChans+AHL_chans+ones(nStripChans,1),(1:nStripChans)'+AHL_chans]];

for i=1:2:length(chans_to_be_bipolar)
    nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nStripChans,1),(1:nStripChans)'+ Depth_electrodes]];
end
strChannelNameList = dataBipolar.label;


%% Comment this if you have bipolar hipp chans
chans_to_be_bipolar = [1 2 9 10 17 18];
nChannelPairs = []
for i = 1:length(chans_to_be_bipolar)
    nChannelPairs = [nChannelPairs; [[repelem(chans_to_be_bipolar(i),length(GridChans))]' [GridChans]]]
end
strChannelNameList = strrep(dataBipolar.label,'m','');
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
                save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100_new.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_1_30_5_100_new.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_retrieval_new.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
                
            end
        end
    end
end


%% Plot the PLV results for every bipolar channel
iSS_to_Plot = 4;
k =0;
% Load the PLV results
if ~GammaFreqPLV
    strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\']

    PLV_Vars = load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100_new.mat'])
    PLV_Pairs_maint = PLV_Vars.PLV_Pairs;
    
%     PLV_Vars = load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_1_30_5_100.mat'])
    PLV_Pairs_encod = PLV_Vars.PLV_Pairs;
    freqAxis = PLV_Vars.freqAxis;
    
    PLV_Vars = load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat'])
    PLV_Pairs_fix = PLV_Vars.PLV_Pairs;
    nChannelPairs = PLV_Vars.nChannelPairs;
    freqAxis_fix = PLV_Vars.freqAxis;
end


for i=1:size(nChannelPairs,1)
    
    if i==1 | ~mod(i-1,nStripChans/3)
        fig = figure;
        stripPL_plot_order = [1 3 2 4];
        ha = tight_subplot(2,2,[.08 .08],[.1 .1],[.1 .1])
        strTitle =['PLV of ',dataBipolar.label(nChannelPairs(i,1)),' and strip',strrep(dataBipolar.label(nChannelPairs(i,2)),'1','') 'ss ',num2str(Set_Sizes{iSS_to_Plot})]
        suptitle(strTitle);
        if i~=1
            k = k+1;
        end
    end
    axes(ha(stripPL_plot_order(i-(k*4))));
    semilogx(freqAxis_fix,abs(PLV_Pairs_fix{i,iSS_to_Plot}),'k','LineWidth',3);
    hold on;
%     semilogx(freqAxis,abs(PLV_Pairs_encod{i,iSS_to_Plot}),'g','LineWidth',3);
    semilogx(freqAxis,abs(PLV_Pairs_maint{i,iSS_to_Plot}),'r','LineWidth',3);
    xlim([4 100]);
    ylim([0 1]);
    ylabel(dataBipolar.label(nChannelPairs(i,2)));
    
end

%% PLV heatmap
FreqBand = [9,18]; % same as granger for maint

[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
PLV_In_Band_SS = [];
iSS_to_Plot =4 %[1:5]
for iPair = 21:24%53:56%1:size(nChannelPairs,1) %193:256%1:384%
    for iSS = iSS_to_Plot
        PLV_temp = PLV_Pairs_maint{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

iSS_ToPLot = 4;
BChan = 44;

BipChan_to_Plot = [21:24];%find(nChannelPairs==BChan);
PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),1,4);

figure;
imagesc(PLV_In_Band_ToPlot,[0 0.6])

%     str = strcat(strrep(sprintf('%s',dataBipolar.label{BChan+2}),'_',' '), ' [ ', num2str(FreqBand(1)),{' '}, ....
%         num2str(FreqBand(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]')
%     title(str)
ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'

colorbar

%% Scalp PLV with Depth Electrodes during fixation
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_fix_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat'],'PLV_Pairs_Scalp_fix_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else%In case you forgot to declare the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat'],'PLV_Pairs_Scalp_fix_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end  


%% Scalp PLV with Depth Electrodes during encoding
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_enc_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end  


%% Scalp PLV with Depth Electrodes during maintenance
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_maint_depth,nChannelPairs_Sclp_depth,freqAxis_sclp,dataBipolar_appended] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100_new.mat'],'PLV_Pairs_Scalp_maint_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_f_1_1_30_5_100_new.mat'],'PLV_Pairs_Scalp_maint_depth','nChannelPairs_Sclp_depth','freqAxis_sclp','dataBipolar_appended','-v7.3');
end  


%% Load the PLV vars
strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_f_1_1_30_5_100_new.mat']);
load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat']);

%% Plot Scalp Depth PLV for one or more bipolar channels
Bipolar_Channels = unique(nChannelPairs_Sclp_depth(:,1));
Bip_chans_to_plot = Bipolar_Channels(1:end);
iSS = 4; % Set Size [6 8]
% Scalp topoplot related
FreqBand_Sclp = [6,7];  
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

for i = 1:numel(Bip_chans_to_plot)
    nFig = figure
    indBipChan = find(nChannelPairs_Sclp_depth==Bip_chans_to_plot(i));
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


%% Scalp PLV with Strip Electrodes during fixation
PLV_Pair_flag = 2; % Value 1 for Hippocampus/Depth Coupling, Value 2 for ECoG grid/strip coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_fix_strip,nChannelPairs_Sclp_strip,freqAxis_sclp_strip,dataBipolar_appended_sclp_strip] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];
if strcmp(coupling_type,'Scalp-Strip')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_1_30_5_100_all_scalp_chans_new.mat'],'PLV_Pairs_Scalp_fix_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Strip','_PLV_Pairs_fixation_f_1_1_30_5_100_all_scalp_chans_new.mat'],'PLV_Pairs_Scalp_fix_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
end  

%% Scalp PLV with Strip Electrodes during encoding
PLV_Pair_flag = 2; % Value 1 for Hippocampus/Depth Coupling, Value 2 for ECoG grid/strip coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_enc_strip,nChannelPairs_Sclp_strip,freqAxis_sclp_strip,dataBipolar_appended_sclp_strip] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];
if strcmp(coupling_type,'Scalp-Strip')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Strip','_PLV_Pairs_encoding_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
end  


%% Scalp PLV with Strip Electrodes during maintenance
PLV_Pair_flag = 2; % Value 1 for Hippocampus/Depth Coupling, Value 2 for ECoG grid/strip coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_maint_strip,nChannelPairs_Sclp_strip,freqAxis_sclp_strip,dataBipolar_appended_sclp_strip] = calculate_PLV_scalp_HL_P38(macro_data,data_Scalp_ses,TrialInformationTable_Scalp,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];
if strcmp(coupling_type,'Scalp-Strip')
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100_all_scalp_chans_new.mat'],'PLV_Pairs_Scalp_maint_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
else %In case you forgot to declare  the coupling type in the parameters and flags section
    save([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Strip','_PLV_Pairs_maintenance_f_1_1_30_5_100_all_scalp_chans_new.mat'],'PLV_Pairs_Scalp_maint_strip','nChannelPairs_Sclp_strip','freqAxis_sclp_strip','dataBipolar_appended_sclp_strip','-v7.3');
end 


%% Load the PLV vars
strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];
load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Strip','_PLV_Pairs_fixation_f_1_1_30_5_100_new.mat']);
load([strVariableFolder,'Patient_',num2str(38,'%.2d '),'_','Scalp-Strip','_PLV_Pairs_maintenance_f_1_1_30_5_100_all_scalp_chans_new.mat']);


%% Plot Scalp Strip PLV for one or more scalp channels
Scalp_Channels = unique(nChannelPairs_Sclp_strip(:,1));
Scalp_chans_to_plot = Scalp_Channels(1:end);
iSS = 4; % Set Size [6 8]
% Scalp topoplot related
FreqBand_Sclp = [7,9];  
[~,indFreq1] = min(abs(freqAxis_sclp_strip-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp_strip-FreqBand_Sclp(2)));

for i = 1:numel(Scalp_chans_to_plot)
    nFig = figure
    indScalpChan = find(nChannelPairs_Sclp_strip==Scalp_chans_to_plot(i));
    ha =  tight_subplot(3,4,[.05 .08],[.1 .08],[.08 .08])
    plot_order = [1:12]%[1 4 7 10 2 5 8 11 3 6 9 12];
    set(gcf,'color','white')
    for j = 1:numel(indScalpChan)
        datavector{j}(i) = median(abs(PLV_Pairs_Scalp_maint_strip{indScalpChan(j),iSS}(indFreq1:indFreq2)));
        Strip_channel = strcat({' '},dataBipolar_appended_sclp_strip{iSS}.label{nChannelPairs_Sclp_strip(indScalpChan(j),2)});
        axes(ha(plot_order(j))); 
        semilogx(freqAxis_sclp_strip,abs(PLV_Pairs_Scalp_fix_strip{indScalpChan(j),iSS}),'k','LineWidth',3);
        hold on;
        semilogx(freqAxis_sclp_strip,abs(PLV_Pairs_Scalp_maint_strip{indScalpChan(j),iSS}),'r','LineWidth',3);
        ylim([0 1]);
        xlim([4 100]);
        ylabel(Strip_channel);
        hold off;
    end
    suptitle(sprintf('Scalp Electrode %s - Strip PLV',dataBipolar_appended_sclp_strip{iSS}.label{Scalp_chans_to_plot(i)}))
    
end
for k = 1:12
    nFig = figure
    Scalp_topoplot(strPaths.Toolboxes,FreqBand_Sclp,freqAxis_sclp_strip,datavector{k},dataBipolarScalp,4,nChannelPairs_Sclp_strip,1)
end
