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
%% Load data for 2 sessions
Session1                   =      [Drive_Letter,'Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_01_Part_01'];
Session2                   =        [Drive_Letter,'Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_02_Part_01'];
data_s1                    =      load(Session1);
data_s2                    =      load(Session2);
TrialInformationTable      =      [data_s1.TrialInformationTable;data_s2.TrialInformationTable] ;
data_s1                    =      data_s1.dataMacro;
data_s2                    =      data_s2.dataMacro;


%% Load scalp data for 2 sessions
Sclp_Session1                    =      [Drive_Letter,'Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_01_Part_01'];
Sclp_Session2                    =      [Drive_Letter,'Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_02_Part_01'];
data_s1_scalp                    =      load(Sclp_Session1);
data_s2_scalp                    =      load(Sclp_Session2);
TrialInformationTable_Scalp      =      [data_s1_scalp.TrialInformationTable;data_s2_scalp.TrialInformationTable] ;
data_scalp_s1                    =      data_s1_scalp.dataScalp;
data_scalp_s2                    =      data_s2_scalp.dataScalp;

%% Merge sessions
cfg = [];
dataBipolar = ft_appenddata(cfg,data_s1,data_s2);
dataBipolarScalp = ft_appenddata(cfg,data_scalp_s1,data_scalp_s2);

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
GridChans = find(contains(strChannelNameList,'GL'));
dataRerefAhipp = reref_toSeparate_Chan(dataBipolar, 4, 'all','avg');
dataRerefGrid = reref_toSeparate_Chan(dataBipolar, 5, 'all','avg');
dataRerefPhipp = reref_toSeparate_Chan(dataBipolar, 4, 'all','avg');

%selection of channels
cfg = [];
cfg.channel = AhippChans;
dataRerefAhipp = ft_selectdata(cfg,dataRerefAhipp);
cfg = [];
cfg.channel = GridChans;
dataRerefGrid = ft_selectdata(cfg,dataRerefGrid);
cfg = [];
cfg.channel = PhippChans;
dataRerefPhipp = ft_selectdata(cfg,dataRerefPhipp);

% append them again
cfg = [];
dataBipolar = ft_appenddata(cfg,dataRerefAhipp,dataRerefGrid);
dataBipolar = ft_appenddata(cfg,dataBipolar,dataRerefPhipp);

%% Reref based on the average of the signals
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =  [3 4] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
dataBipolar       = ft_preprocessing(cfg,dataBipolar);
%% Laplace
dataBipolar.label = strrep(dataBipolar.label,'_','');

cfg = []
cfg.channel = 'all';
cfg.reref = 'yes';
cfg.refmethod = 'laplace';
cfg.refchannel = [9:72];
cfg.groupchans = 'yes';
dataBipolar = ft_preprocessing(cfg,dataBipolar);

%% Apply montage
clear montage;
montage.labelold = dataBipolar.label;
num_bipolar_chans=6;
num_reference_chans=size(montage.labelold,1);
%prepare the montage matrix
montage_matrix = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
chans_to_be_bipolar = [1 2 1 3 2 3 73 74 73 75 74 75 ] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end
%Append the bipolar channel labels to the reference channels labels
for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew = {dataBipolar.label{:},montage.labelnew{:}};
montage.tra = montage_matrix;
% dataBipolar = ft_apply_montage(dataBipolar,montage);

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod  = 'bipolar'
cfg.refchannel = [82 85]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));
%% Bipolar grid
clear montage
nBipolarChannelList = sort([9:72-1 10:72]);
a = [([8:8:56]*2)-1 ([8:8:56])*2];
chans_to_be_bipolar = nBipolarChannelList(setdiff(1:length(nBipolarChannelList),a))
montage.labelold = dataBipolar.label';
nGridChans = 64;
num_bipolar_chans = nGridChans/8*7;
num_reference_chans=size(montage.labelold,1);
montage_matrix = eye(num_reference_chans-nGridChans+num_bipolar_chans,num_reference_chans);
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end

for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew = [dataBipolar.label{1:8} montage.labelnew  dataBipolar.label{73:86}];
montage.tra = montage_matrix;

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod  = 'bipolar'
cfg.refchannel = [9 64]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);



%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500%80;%500 ;
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


%% Scalp EEG Preprocessing
if EEG_Preproc
    %% Visualize scalp raw data
    cfg =[];
    cfg.viewmode = 'vertical';
    cfg.preproc.lpfilter = 'yes'
    cfg.preproc.lpfreq = 45;
    ft_databrowser(cfg,dataBipolarScalp);
    
    
    %%
    cfg = [];
    cfg.channel = {'all' '-_A1','-_A2','-_T1','-_T2','-_C3','-_C4','-_Cz','-_O1','-_P3','-_Submm','-_Submp'}
    dataBipolarScalp_Rej_Chans = ft_preprocessing(cfg,dataBipolarScalp);
    
    %% Reject trial
    cfg          = [];
    cfg.method   = 'trial';
    % cfg.alim     = 1e-12;
    dummy        = ft_rejectvisual(cfg,dataBipolarScalp_Rej_Chans);
    
    cfg =[];
    cfg.viewmode = 'vertical';
    cfg.preproc.lpfilter = 'yes'
    cfg.preproc.lpfreq = 45;
    ft_databrowser(cfg,dataBipolarScalp);
    
    %% PSD of scalp eeg after removal of noisy channels
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:120;
    cfg.tapsmofrq   = 2;
    % cfg.channel = 1:23;
    fr_maint_sclp = ft_freqanalysis(cfg,dataBipolarScalp);
    
    fig = figure;
    plot(fr_maint_sclp.freq,10*log10(fr_maint_sclp.powspctrm),'LineWidth',2);
    fig = figure;
    plot(fr_maint_sclp.freq,fr_maint_sclp.powspctrm,'LineWidth',2);
    fr_maint_sclp.label = strrep(fr_maint_sclp.label,'_',' ');
    legend(fr_maint_sclp.label)
    
    
    %% ICA
    cfg = [];
    cfg.method       = 'runica'
    cfg.channel      = {'all' '-_Submm','-_Submp'};
    cfg.trials       = setdiff([1:77],[8 14 34]); % Reject trials: 8,14,34
    cfg.numcomponent = 'all';
    cfg.demean       = 'yes';
    cfg.feedback     =  'text'
    IC_components    = ft_componentanalysis(cfg,dataBipolarScalp)
    
    IC_components.topolabel = strrep(IC_components.topolabel,'_','');
    IC_components.cfg.channel = strrep(IC_components.cfg.channel,'_','');
    save([strPaths.Data,'IC_components_scalp_EEG'],'IC_components');
    
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
    
    
    %% Reject components after visual evaluation
    cfg = [];
    cfg.component = [1 2 3 4 5 8 10 21 23];
    cfg.demean = 'yes';
    dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);
    %%
    
    cfg = [];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg,dataBipolarScalp_Reconstructed);
    
    
    %% PSD after reconstruction - maintenance
    cfg = [];
    cfg.latency = [-2,-1/dataBipolarScalp_Reconstructed.fsample];
    dataBipolar_Scalp_SS = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.channel     =  {'all'}% '-C3','-C4','-Cz','-O1','-P3','-T1','-T2'}
    cfg.foi         = 0.5:0.5:120;
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    fr_maint_sclp_recons = ft_freqanalysis(cfg,dataBipolar_Scalp_SS);
    freqAxis = fr_maint_sclp_recons.freq;
    
    
    figure;
    
    plot(freqAxis,10*log10(fr_maint_sclp_recons.powspctrm))
    legend(fr_maint_sclp_recons.label)
    
    
    
    freqBand = [6 9]
    %     freqAxis = fr_maint_sclp_recons.freq;
    [~,indFreq1] = min(abs(freqAxis-freqBand(1)));
    [~,indFreq2] = min(abs(freqAxis-freqBand(2)));
    datavector = median(fr_maint_sclp_recons.powspctrm(:,indFreq1:indFreq2)')
    
    Scalp_topoplot(strPaths.Toolboxes,freqBand,freqAxis,datavector,fr_maint_sclp_recons,4,1,1)
    
    
    
    %% PSD after reconstruction - fixation,maintenance,encoding
    
    
    
    %% Select only correct trials
    [dataBipolar_Scalp_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolarScalp_Reconstructed,TrialInformationTable);
    
    %% Divide data into set sizes
    %   [4] -> Set Size 1
    %       [6] -> Set Size 2
    %           [8] -> Set Size 3
    %               [6 8] -> Set Size 4
    %                   [4 6 8] -> Set Size 5
    Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
    [dataBipolar_Scalp_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_Scalp_SS,TrialInformationTable);
    
    
    %%
    cfg = [];
    cfg.latency = [-2,-1/dataBipolarScalp_Reconstructed.fsample];
    dataBipolar_Scalp_SS_maint = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    dataBipolar_Scalp_SS_maint_hi = ft_selectdata(cfg,dataBipolar_Scalp_SS{5});
    dataBipolar_Scalp_SS_maint_low = ft_selectdata(cfg,dataBipolar_Scalp_SS{1});
    
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolarScalp_Reconstructed.fsample];
    dataBipolar_Scalp_SS_enc = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    dataBipolar_Scalp_SS_enc_hi = ft_selectdata(cfg,dataBipolar_Scalp_SS{5});
    dataBipolar_Scalp_SS_enc_low = ft_selectdata(cfg,dataBipolar_Scalp_SS{1});
    
    
    
    cfg = [];
    cfg.latency = [-6,-5-1/dataBipolarScalp_Reconstructed.fsample];
    dataBipolar_Scalp_SS_fix = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    dataBipolar_Scalp_SS_fix_hi = ft_selectdata(cfg,dataBipolar_Scalp_SS{5});
    dataBipolar_Scalp_SS_fix_low = ft_selectdata(cfg,dataBipolar_Scalp_SS{1});
    
    dataBipolarScalp_Reconstructed2 = dataBipolarScalp_Reconstructed;
    dataBipolarScalp_Reconstructed2.label = strrep(strcat('_',dataBipolarScalp_Reconstructed.label),'_',' ');
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.channel     =  {'all'}
    cfg.foi         = 4:0.5:120;
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    fr_maint_sclp_recons = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_maint);
    fr_enc_sclp_recons = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_enc);
    fr_fix_sclp_recons = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_fix);
    
    fr_maint_sclp_recons_maint_hi = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_maint_hi);
    fr_maint_sclp_recons_maint_low = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_maint_low);
    fr_enc_sclp_recons_hi = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_enc_hi);
    fr_fix_sclp_recons_hi = ft_freqanalysis(cfg,dataBipolar_Scalp_SS_fix_hi);
    
    freqAxis = fr_maint_sclp_recons.freq;
    %%
    figure
    iSS_ToPLot = 4;
    set(gcf,'color','white')
    LogF_power_maint = logspace(fr_maint_sclp_recons.freq(1),fr_maint_sclp_recons.freq(end),length(fr_maint_sclp_recons.freq));
    LogF_power_fix = logspace(fr_fix_sclp_recons.freq(1),fr_fix_sclp_recons.freq(end),length(fr_fix_sclp_recons.freq));
    %     LogF_power_enc = logspace(fr_maint_sclp_recons.powspctrm(15,1),fr_maint_sclp_recons.powspctrm(15,end),length(fr_maint_sclp_recons.powspctrm));
    
    
    %     ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
    figLayout = Get_Figure_Scalp_Layout_Information();
    nScalpChans = size(dataBipolarScalp_Reconstructed.label,1);
    for ii = 15%1:nScalpChans%15%
        Scalp_channel = dataBipolarScalp_Reconstructed2.label{ii};
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_psd(ii) = figLayout{indSubplot,2};
        %         axes(ha(nSubplot_psd(ii)));
        semilogx(fr_maint_sclp_recons.freq,10*log10(fr_maint_sclp_recons.powspctrm(ii,:)),'k','LineWidth',2);
        %         semilogx(LogF_power_maint,10*log10(fr_maint_sclp_recons.powspctrm(ii,:)),'k','LineWidth',2);
        
        hold on;
        semilogx(fr_enc_sclp_recons.freq,10*log10(fr_enc_sclp_recons.powspctrm(ii,:)),'g','LineWidth',2);
        %         semilogx(LogF_power_maint,10*log10(fr_enc_sclp_recons.powspctrm(ii,:)),'g','LineWidth',2);
        
        hold on;
        semilogx(fr_fix_sclp_recons.freq,10*log10(fr_fix_sclp_recons.powspctrm(ii,:)),'r','LineWidth',2);
        %         semilogx(LogF_power_fix,10*log10(fr_fix_sclp_recons.powspctrm(ii,:)),'r','LineWidth',2);
        
        
        ylim([-30 30]);
        xlim([4 100])
        elec_pair = fr_maint_sclp_recons.label{ii};
        ylabel(elec_pair);
    end
    %     set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    %     set(ha(1:size(ha,1)),'FontSize',18);
    suptitle(['Scalp Channels Power Spectrum Maintenance vs Encoding vs Fixation, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
    for i = 1:size(ha,1)
        
        if ~ismember(i,nSubplot_psd)
            set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    
    %% PSD Scalp different periods High vs Low Workload
    
    figure
    iSS_ToPLot = 4;
    set(gcf,'color','white')
    %     ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
    figLayout = Get_Figure_Scalp_Layout_Information();
    nScalpChans = size(dataBipolarScalp_Reconstructed.label,1);
    for ii = 15%1:nScalpChans%15%
        Scalp_channel = dataBipolarScalp_Reconstructed2.label{ii};
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_psd(ii) = figLayout{indSubplot,2};
        %         axes(ha(nSubplot_psd(ii)));
        semilogx(fr_maint_sclp_recons_maint_hi.freq,10*log10(fr_maint_sclp_recons_maint_hi.powspctrm(ii,:)),'k','LineWidth',2);
        hold on;
        semilogx(fr_maint_sclp_recons_maint_low.freq,10*log10(fr_maint_sclp_recons_maint_low.powspctrm(ii,:)),'b','LineWidth',2);
        hold on;
        semilogx(fr_enc_sclp_recons_hi.freq,10*log10(fr_enc_sclp_recons.powspctrm(ii,:)),'g','LineWidth',2);
        hold on;
        semilogx(fr_fix_sclp_recons_hi.freq,10*log10(fr_fix_sclp_recons.powspctrm(ii,:)),'r','LineWidth',2);
        
        ylim([-25 30]);
        %         xlim([3 100])
        elec_pair = fr_maint_sclp_recons.label{ii};
        ylabel(elec_pair);
    end
    %     set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    %     set(ha(1:size(ha,1)),'FontSize',18);
    suptitle(['Scalp Channels Power Spectrum Maintenance vs Encoding vs Fixation, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
    for i = 1:size(ha,1)
        
        if ~ismember(i,nSubplot_psd)
            set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    
    
    
    
    %% Clean EEG signal
    cfg = [];
    cfg.channel = {'all'} % '-C3','-C4','-Cz','-O1','-P3','-T1','-T2'}
    dataBipolarScalp_Reconstructed = ft_selectdata(cfg,dataBipolarScalp_Reconstructed);
    
    %%
    cfg = [];
    cfg.viewmode = 'vertical';
    ft_databrowser(cfg,dataBipolarScalp_Reconstructed);
    
end

%% Evoked Potential Visual Stimulus

if Visual_Evoked_Potential
    
    cfg = []
    cfg.trial = 'all';
    cfg.latency  = [-5.1 -4.1]
    timelock_dataBipolar = ft_timelockanalysis(cfg,dataBipolar);
    
    for i=1:size(timelock_dataBipolar.label,1)
        timelock_dataBipolar.label{i,1} = strrep(sprintf('%s',timelock_dataBipolar.label{i,1}),'_',' ');
        
    end
    % Bipolar Channels
    figure(20)
    plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(81:84,:));
    legend(timelock_dataBipolar.label{81:84})
    % ylim([-140 60])
    
    % Grid Channels A1 to D8
    figure(30)
    plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(9:40,:));
    legend(timelock_dataBipolar.label{9:40})
    
    %Grid Channels E1 to H8
    figure(40)
    plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(41:72,:));
    legend(timelock_dataBipolar.label{41:72})
    
    %% Plot average of every grid electrode [-5.1 4.1] sec
    figure;
    ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    tRange = [find(timelock_dataBipolar.time == -4.8) find(timelock_dataBipolar.time == -4.68)]
    plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
    for ii = 1:64
        
        axes(ha(plot_order(ii))); plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(ii+AHL_chans,:));
        
        elec_pair = timelock_dataBipolar.label{ii+AHL_chans};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([-50 50])
        set(ha(1:64),'XTickLabel',''); set(ha,'YTickLabel','')
        Min_Value(ii) = min(timelock_dataBipolar.avg(ii+AHL_chans,tRange(1):tRange(2)));
    end
    
    %Min Value of the Grid Channels VEP
    Min_Value = reshape(Min_Value,8,8)
    
    figure;
    im=imagesc(Min_Value);
    str = strcat('Minimum Value of Evoked Potential - Grid Channels')
    title(str)
    ax = gca;
    ax.XTick = 1:8;
    ax.YTick = 1:8;
    ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
    ax.YTickLabel = 1:8;
    ax.YDir = 'normal'
    colorbar
    colormap(flipud(jet))
    
    %Hippocampus channels
    figure;
    ha = tight_subplot(2,8,[.01 .03],[.1 .01],[.01 .01])
    for ii = 1:16
        if ii<=8
            axes(ha(ii)); plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(ii,:));
            elec_pair = timelock_dataBipolar.label{ii};
            ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
            Min_ValueH(ii) = min(timelock_dataBipolar.avg(ii,tRange(1):tRange(2)));
            
        else
            axes(ha(ii)); plot(timelock_dataBipolar.time,timelock_dataBipolar.avg(ii+nGridChans,:));
            elec_pair = timelock_dataBipolar.label{ii+nGridChans};
            ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
            Min_ValueH(ii) = min(timelock_dataBipolar.avg(ii+nGridChans,tRange(1):tRange(2)));
            
        end
        ylim([-50 50])
        set(ha(1:16),'XTickLabel',''); set(ha,'YTickLabel','')
    end
    
    %Min Value of the Hippocampus Channels' EP
    Min_ValueH = reshape(Min_ValueH',8,2)
    
    figure;
    im=imagesc(Min_ValueH);
    str = strcat('Minimum Value of Evoked Potential - Hippocampus Channels')
    title(str)
    ax = gca;
    ax.YTick = 1:8;
    ax.XTick = 1:2;
    ax.YTickLabel = 1:8;
    ax.XTickLabel = {'AHL', 'PHL'};
    ax.YDir = 'normal'
    colorbar
    colormap(flipud(jet))
end



%% Extract 2 seconds of maintenance
if strcmp(TaskPeriod,'maint')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
%     %%
%     %Power spectrum of maintenance
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = [1:40]%100];
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{5});
    fr_maint_low =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    plotColor = {'r','g','b','c','m','k'}
    % for hippocampal channels
    for i=1:size(fr_maint.label,2) fr_maint.label{i} = strrep(sprintf('%s',fr_maint.label{i}),'_',' '); end
    figure;
    for i = Bip_chans(1):Bip_chans(end)
        plot(fr_maint.freq,10*log10(fr_maint.powspctrm(i,:)),plotColor{i-80},'LineWidth',2);
        hold on;
        title('Power Spectrum of Hippocampal Channels during Maintenance Set Size [6 8]');
        ylabel('Power 10*log10(\muV^2/Hz)')
        xlabel('Frequency(Hz)')
    end
    legend(fr_maint.label{Bip_chans(1):Bip_chans(end)})
%     
elseif strcmp(TaskPeriod,'encod')
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-4,-3-1/dataBipolar.fsample];
        %         cfg.latency = [-5,-4-1/dataBipolar.fsample];
        
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 1:30;
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_enc =ft_freqanalysis(cfg,dataBipolar_Ret_SS{5});
    plotColor = {'r','g','b','c','m','k'}
    
    % for hippocampal channels
    for i=1:size(fr_enc.label,2) fr_enc.label{i} = strrep(sprintf('%s',fr_enc.label{i}),'_',' '); end
    figure;
    for i = Bip_chans(1):Bip_chans(end)
        plot(fr_enc.freq,10*log10(fr_enc.powspctrm(i,:)),plotColor{i-80},'LineWidth',2);
        hold on;
        title('Power Spectrum of Hippocampal Channels during Encoding Set Size [6 8]');
        ylabel('Power 10*log10(\muV^2/Hz)')
        xlabel('Frequency(Hz)')
    end
    legend(fr_enc.label{Bip_chans(1):Bip_chans(end)})
    
    
    
elseif strcmp(TaskPeriod,'fix')
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
    cfg.foi         = [1:40];
    %     cfg.foi         = 80:0.5:200;
    
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_fix =ft_freqanalysis(cfg,dataBipolar_Ret_SS{5});
    
    
else
    %% Extract 2 seconds of retrieval latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [0,2-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
end



%% PSD Maint vs Encoding
if PSD_calc
    %Grid Channels
    figure;
    ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
    Maint_Grid_spctm = fr_maint.powspctrm(9:72,:); %grid chans
    Enc_Grid_spctm = fr_enc.powspctrm(9:72,:); %grid chans
    
    for i = 1:64
        axes(ha(plot_order(i)));
        %     m=plot(fr_maint.freq,10*log10(Maint_Grid_spctm(i,:)),'r');
        m=semilogx(fr_maint.freq,10*log10(Maint_Grid_spctm(i,:)),'k');
        hold on;
        %     e=plot(fr_enc.freq,10*log10(Enc_Grid_spctm(i,:)),'g');
        e=semilogx(fr_enc.freq,10*log10(Enc_Grid_spctm(i,:)),'g');
        
        f =semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(i,:)),'r');
        %     title('Power Spectrum of Channels Maintenance vs Encoding Set Size [6 8]');
        %     ylabel('Power 10*log10(\muV^2/Hz)')
        %     xlabel('Frequency(Hz)')
        elec_pair = fr_maint.label{i+8};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([-30,0])
        %         ylim([-0.1 0.4])
        %         ylim([-10,500])
        %     legend([m(1),e(1)],'Maintenance','Encoding')
    end
    %     set(ha(1:64),'XTickLabel',''); set(ha,'YTickLabel','')
    suptitle('Power Spectrum of Grid Channels Maintenance vs Encoding vs Fixation Set Size [6 8]');
    
    
    %% Plot
    freqAxis = [0.5:0.5:30]
    FreqBand = [8,10];
    
    [~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
    [~,indFreq2] = min(abs(freqAxis-FreqBand(2)));
    clear PSD_maint
    for i = 1:64
        PSD_temp = Maint_Grid_spctm(i,indFreq1:indFreq2);
        PSD_maint(i,:) = PSD_temp;
        
    end
    Median_PSD_maint = median(PSD_maint');
    Median_PSD_maint = reshape(Median_PSD_maint,8,8);
    figure;imagesc(Median_PSD_maint,[50,300])
    str = strcat('Median Linear PSD,' , ' Set Size',{' ['}, num2str(Set_Sizes{4}), ' ]', ', Maintenance, Frequency Band [ ',num2str(FreqBand(1)),'-', num2str(FreqBand(2)),' ] Hz' )
    title(str)
    ax = gca;
    ax.XTick = 1:8;
    ax.YTick = 1:8;
    ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
    ax.YTickLabel = 1:8;
    ax.YDir = 'normal'
    
    colorbar
    colormap jet
    %%
    %Hippocampus Channels
    figure;
    ha = tight_subplot(6,4,[.01 .03],[.1 .01],[.01 .01])
    HC_Channels = [1:8 73:86]
    Maint_HC_spctm = fr_maint.powspctrm;
    Enc_HC_spctm = fr_enc.powspctrm;
    Fix_HC_spctm = fr_fix.powspctrm;
    
    for i = 1:size(HC_Channels,2)
        axes(ha(i));
        m=plot(fr_maint.freq,10*log10(Maint_HC_spctm(HC_Channels(i),:)),'r');
        hold on;
        e=plot(fr_enc.freq,10*log10(Enc_HC_spctm(HC_Channels(i),:)),'g');
        %     title('Power Spectrum of Channels Maintenance vs Encoding Set Size [6 8]');
        %     ylabel('Power 10*log10(\muV^2/Hz)')
        %     xlabel('Frequency(Hz)')
        elec_pair = fr_maint.label{HC_Channels(i)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([-10,35])
        %     legend([m(1),e(1)],'Maintenance','Encoding')
    end
    set(ha(1:size(HC_Channels,2)),'XTickLabel',''); set(ha,'YTickLabel','')
    suptitle('Power Spectrum of Hippocampus Channels Maintenance vs Encoding Set Size [6 8]');
    %
    figure;
    plotColor = {'r','g','b','k'}
    plotColor2 = {'--r','--g','--b','--k'}
    plotColor3 = {'-.r*','-.g*','-.b*','-.k*',}
    chans_to_plot = [82 83 85 86]
    for i = 2%1:size(chans_to_plot,2)
        %        P(1) = semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(chans_to_plot(i),:)),plotColor{i},'LineWidth',2);
        PL(i) = semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(chans_to_plot(i),:)),plotColor{i},'LineWidth',2);
        hold on;
        %         semilogx(fr_maint_low.freq,10*log10(fr_maint_low.powspctrm(chans_to_plot(i),:)),plotColor{i},'LineWidth',2);
        title('Power Spectrum of Hippocampus Channels Maintenance vs Encoding vs Fixation Set Size [6 8]');
        ylabel('Power 10*log10(\muV^2/Hz)');
        xlabel('Frequency(Hz)');
    end
    for i =2%1:size(chans_to_plot,2)
        P(2) = semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(chans_to_plot(i),:)),plotColor2{i},'LineWidth',2);
        
        hold on;
        P(3) = plot(fr_fix.freq,10*log10(fr_fix.powspctrm(chans_to_plot(i),:)),plotColor3{i},'LineWidth',2);
        xlim([4 100])
    end
    
    for i= Bip_chans(1):Bip_chans(end) fr_maint.label{i}= strrep(sprintf('%s',fr_maint.label{i}),'_',' '); end
    L{1} = 'maintenance'
    L{2} = 'encoding'
    L{3} = 'fixation'
    %     legend(P,L);
    %     legend(PL,fr_maint.label{chans_to_plot})
    
end


%% Channel pairs to run
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;
chans_to_be_bipolar = [1:3 73:75]; %% Comment this if you have bipolar hipp chans
if PLV_depth_electrodes == 1
    
    nChannelPairs = [81 83;81 84; 82 83; 82 84];
    
else
    nGridChans = 64;
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
    %% Comment this if you have bipolar hipp chans
    nChannelPairs = []
    for i = 1:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs; [[repelem(chans_to_be_bipolar(i),length(GridChans))] ; [GridChans]']']
    end
    %%
    nChannelPairs = [[ones(nGridChans,1),(1:nGridChans)'+AHL_chans];[nGridChans+AHL_chans+ones(nGridChans,1),(1:nGridChans)'+AHL_chans]];
    
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChans,1),(1:nGridChans)'+ AHL_chans]];
    end

end
strChannelNameList = dataBipolar.label;
% nPairs_att1 = [ones(5,1)+83,[10 11 13 18 19]']
%% Loop over pairs
% PLV_Pairs = {};
% Coh_Pairs = {};
tStart = tic;
for iPair = [266 267 275 274]%257:320%1:size(nChannelPairs,1)
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
    %     tic
    %     clear PLV_SS Coh_SS
    for iSS = 1:nSet_Size
        %% Frequency analysis
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'fourier';
%         cfg.pad= 20;
        cfg.foi         = [1:40]%[1:29 30:5:100];
            
%         end
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
    %     toc
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
%% Save PLV Pairs
if 0
    if Grid_Hippocampus_coupling
        strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\']
        mkdir(strVariableFolder);
        
        if ~GammaFreqPLV
            if strcmp(TaskPeriod,'maint')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_test_f_30_5_100.mat'],'PLV_Pairs','nChannelPairs','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_retrieval.mat'],'PLV_Pairs','nChannelPairs','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_test_f_30_5_100.mat'],'PLV_Pairs','nChannelPairs','-v7.3')
                
            end
        else
            if strcmp(TaskPeriod,'maint')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_gamma_smooth.mat'],'PLV_Pairs','nChannelPairs','freqAxis','TaskPeriod','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_gamma_smooth.mat'],'PLV_Pairs','nChannelPairs','freqAxis','TaskPeriod','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_retrieval_gamma_smooth.mat'],'PLV_Pairs','nChannelPairs','freqAxis','TaskPeriod','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_gamma_smooth.mat'],'PLV_Pairs','nChannelPairs','freqAxis','TaskPeriod','-v7.3')
                
            end
            
        end
    end
end
%% Load PLV_Pairs for a given period on Grid Hippocampus coupling

if Grid_Hippocampus_coupling && strcmp(TaskPeriod,'maint')
    if ~exist('freqAxis')
        %         freqAxis = 0.5:0.5:30;
    end
    if ~exist('PLV_Pairs')
        
        PLV_Pairs_struct =  load('F:/Vasileios/Task Analysis/Data/Analysis Data/Grid_Hippocampus coupling/Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance.mat');
        PLV_Pairs = PLV_Pairs_struct.PLV_Pairs;
        nChannelPairs = PLV_Pairs_struct.nChannelPairs;
        
    end
elseif Grid_Hippocampus_coupling && strcmp(TaskPeriod,'encod')
    PLV_Pairs_struct =  load('F:/Vasileios/Task Analysis/Data/Analysis Data/Grid_Hippocampus coupling/Patient_42_Grid-Hippocampus_PLV_Pairs_encoding.mat');
    PLV_Pairs = PLV_Pairs_struct.PLV_Pairs;
    nChannelPairs = PLV_Pairs_struct.nChannelPairs;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualizations
%% Plot
iPair = 82%130;
figure;
for iSS = 4%1:nSet_Size
    strChannelNames = dataBipolar_Ret_SS{iSS}.label(nChannelPairs(iPair,:));%label(nChannelPairs(iPair,:));
    semilogx(freqAxis,abs(PLV_Pairs{iPair,iSS}),'Color','r','LineWidth',3)
    hold on%   
    xlim([4,100])
    
    ylim([0,1])
    %     set(gca,'XTick',[1,5:5:30])
%     title(sprintf('%s %s %s %s',strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'),', Encoding,', 'Set Size,',num2str(Set_Sizes{iSS})))
end

%%
FreqBand = [4,8];
% FreqBand = [18,30];

[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));

%% Plot
PLV_In_Band_SS = [];
iSS_to_Plot =4; %[1:5]
% taskperiod = 'PLV_Pairs_maint';
% cell elements PLV_Pairs_enc
for iPair = 257:320%129:size(nChannelPairs,1) %193:256%1:384%
    for iSS = iSS_to_Plot
        %         switch taskperiod
        %             case 'PLV_Pairs_encode'
        %                 PLV_temp = PLV_Pairs_enc{iPair,iSS}(indFreq1:indFreq2);
        %             case 'PLV_Pairs_maint'
        PLV_temp = PLV_Pairs{iPair,iSS}(indFreq1:indFreq2);
        %         end
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

%% Plot the heatmap of the PLV for a certain bipolar channel and for a certain set size
%% Important
iSS_ToPLot = 4;
BChan = 83; %Select the Bipolar Channel for every pair of which with the grid chans you want to plot the PLV values
if PLV_depth_electrodes
    BipChan_to_Plot = find(nChannelPairs==BChan);
    PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
    PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),1,2);
    
    figure;
    imagesc(PLV_In_Band_ToPlot)
    
    str = strcat(strrep(sprintf('%s',dataBipolar.label{BChan}),'_',' '), ' [ ', num2str(FreqBand(1)),{' '}, ....
        num2str(FreqBand(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]')
    title(str)
    ax = gca;
    ax.XTick = 1:2;
    ax.YTick = 1;
    ax.YTickLabel = {timelock_dataBipolar.label{BChan}};
    ax.XTickLabel = {'PHL1-PHL2','PHL2-PHL3'};
    ax.YDir = 'normal'
    
    colorbar
    colormap jet
else
    BipChan_to_Plot = find(nChannelPairs==BChan);
    PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
    PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),8,8);
    
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
%     colormap jet
end
%% Plot the PLV spectrum of the bipolar channels for every combination with grid channels for a certain set size
figure
iSS_ToPLot = 4; %Set Size to Calculate the Median

ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';

for ii = 1:64
    for j = 2%1:num_bipolar_chans
        axes(ha(plot_order(ii))); plot(freqAxis,abs(PLV_Pairs_enc{ii+((j+1)*64),iSS_ToPLot}),strPlotColors{j});
        temp = abs(PLV_Pairs_enc{ii+((j+1)*64),iSS_ToPLot}(indFreq1:indFreq2));
        Data_FOI(ii,j)=median(temp);
        elec_pair = dataBipolar.label{nChannelPairs(ii+((j+1)*64),2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([0 1])
        hold on;
    end
    set(ha(1:64),'XTickLabel',''); set(ha,'YTickLabel','')
end
ttl = ['PLV spectrum of Bipolar Hippocampus channels on ECoG Grid during Encoding Set Size ',num2str(Set_Sizes{iSS_to_Plot})]
suptitle(ttl)
Median_Data_FOI = reshape(median(Data_FOI'),8,8);
MedianF{iSS_ToPLot} = Median_Data_FOI;

%% Plot the median of Bipolar channels' PLV with the grid channels for a certain set size
iSS_ToPLot = 4; %Set Size to Plot
figure;

im=imagesc(MedianF{iSS_ToPLot},[0 1]);
str = strcat('Median PLV,' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ],', ' n=', num2str(size(dataBipolar_Ret_SS{iSS_ToPLot}.trial,2)),', Encoding')
title(str)
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

%% Hippocampus - Grid PLV maintenance vs encoding
iSS_to_Plot = [4];
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for j = 385:448%1:size(nChannelPairs,1)
    if j == 385%1||mod(j,64) == 1
        figure;
        ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01]);
        tlt = [strrep(sprintf('%s',dataBipolar_Ret_SS{iSS_to_Plot}.label{nChannelPairs(j,1)}),'_',' '),' Hippocampus-ECoG coupling ','Maintenance vs Encoding ','Set Size ',num2str(Set_Sizes{iSS_to_Plot})];
        suptitle(tlt)
    end
    for iSS = iSS_to_Plot
        if j<=64
            axes(ha(plot_order(j)));
        else
            axes(ha(plot_order(j-floor((j-1)/64)*64)))
        end
        plot(freqAxis,abs(PLV_Pairs_maint{j,iSS}),'r','LineWidth',2);  ylim([0,1]);
        hold on;
        plot(freqAxis,abs(PLV_Pairs_enc{j,iSS}),'g','LineWidth',2);    ylim([0,1]);
        elec_pair = dataBipolar_Ret_SS{iSS_to_Plot}.label{nChannelPairs(j,2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        
    end
    set(ha(1:nGridChans),'XTickLabel',''); set(ha,'YTickLabel','')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Scalp-/ECoG/Hippocampus Coupling Maintenance/Encoding
%Scalp data
cfg               = []
cfg.channel       = {'all', '-_Submm', '-_Submp'}
cfg.reref         =  'yes'
cfg.refchannel    =   {'A1' 'A2'}
cfg.refmethod     =  'median'
dataBipolar_Sclp       = ft_preprocessing(cfg,dataBipolarScalp);


for i=1:size(dataBipolar_Sclp.label,1)
    dataBipolar_Sclp.label{i} = strrep(sprintf('%s',dataBipolar_Sclp.label{i} ),'_',' ')
end

cfg = [];
cfg.keepsampleinfo = 'no';
dataBipolar_sclp_grid = ft_appenddata(cfg,dataBipolar_Sclp,dataBipolar);

for i=1:size(dataBipolar_sclp_grid.label,1)
    dataBipolar_sclp_grid.label{i} = strrep(sprintf('%s',dataBipolar_sclp_grid.label{i} ),'_',' ')
end

%% Select only correct trials
[dataBipolar_sclp_grid_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar_sclp_grid,TrialInformationTable_Scalp);

%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_sclp_grid_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_sclp_grid_SS,TrialInformationTable);



%% Extract 2 seconds of maintenance
if strcmp(TaskPeriod,'maint')
    nSet_Size = size(dataBipolar_sclp_grid_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_sclp_grid_SS{iSS}.fsample];
        dataBipolar_scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_sclp_grid_SS{iSS});
    end
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:30;
    cfg.tapsmofrq   = 2;
    cfg.channel = 1:23;
    fr_maint_sclp = ft_freqanalysis(cfg,dataBipolar_scalp_Ret_SS{5});
elseif strcmp(TaskPeriod,'encod')
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_sclp_grid_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar.fsample];
        dataBipolar_scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_sclp_grid_SS{iSS});
    end
    %-------------------------------------------------------------------------
    %Power spectrum of scalp channels
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:30;
    cfg.tapsmofrq   = 2;
    cfg.channel = 1:23;
    fr_enc_sclp =ft_freqanalysis(cfg,dataBipolar_scalp_Ret_SS{5});
end
%% PSD scalp maint vs enc
if ~exist('fr_enc_sclp')
    error('Power Spectrum of Scalp Channels is not calculated. Please calculate it !')
    
else
    figure
    set(gcf,'color','white')
    ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
    figLayout = Get_Figure_Scalp_Layout_Information();
    nScalpChans = size(dataBipolar_Sclp.label,1);
    for ii = 1:nScalpChans
        Scalp_channel = dataBipolar_sclp_grid.label{ii};
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_psd(ii) = figLayout{indSubplot,2};
        axes(ha(nSubplot_psd(ii)));
        plot(fr_maint_sclp.freq,10*log10(fr_maint_sclp.powspctrm(ii,:)),'r','LineWidth',2);
        hold on;
        plot(fr_enc_sclp.freq,10*log10(fr_enc_sclp.powspctrm(ii,:)),'g','LineWidth',2);
        ylim([-15 30]);
        elec_pair = fr_enc_sclp.label{ii};
        ylabel(elec_pair);
    end
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    set(ha(1:size(ha,1)),'FontSize',18);
    % suptitle(['Scalp Channels Power Spectrum Maintenance vs Encoding, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
    for i = 1:size(ha,1)
        
        if ~ismember(i,nSubplot_psd)
            set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    
    % Plot PSD of maintenance in subplot
    freqBand = [8 10];
    freqAxis_sclp = fr_maint_sclp.freq;
    [~,indFreq1] = min(abs(freqAxis_sclp-freqBand(1)));
    [~,indFreq2] = min(abs(freqAxis_sclp-freqBand(2)));
    
    figure;
    set(gcf,'color','white')
    ha = tight_subplot(6,4,[.01 .03],[.1 .01],[.01 .01]);
    nScalpChans = size(dataBipolar_Sclp.label,1);
    for ii = 1:nScalpChans
        axes(ha(ii));
        plot(fr_maint_sclp.freq,10*log10(fr_maint_sclp.powspctrm(ii,:)),'r','LineWidth',2);
        PSD_band = mean(10*log10(fr_maint_sclp.powspctrm(ii,indFreq1:indFreq2)));
        ylim([-15 30]);
        elec_pair = fr_enc_sclp.label{ii};
        strLbl = sprintf('%s %2.1f',elec_pair, PSD_band);
        ylabel(strLbl);
    end
    strTitle = ['Power Spectral Density during Maintenance']
    suptitle(strTitle)
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    set(ha(1:size(ha,1)),'FontSize',9.5);
    set(ha(end),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1]);
end
%Calculation of PLV in Scalp-Hippocampus coupling either during maintenance or
%encoding and comparison between the two periods
%% Scalp PLV with  Hippocampus during maintenance
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_maint_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_cleanEEG_f_30:5_100.mat'],'PLV_Pairs_Scalp_maint_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_cleanEEG_f_30_5_100.mat'],'PLV_Pairs_Scalp_maint_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end



%% Prepare scalp topoplot
FreqBand_Sclp = [6,7];
nScalpChans = size(data_s1_scalp.dataScalp.label,1)-2;
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));
iSS_High = 4;
for j = 3%1:6
    for i = 1:nScalpChans
        datavector{j}(i) = median(abs(PLV_Pairs_Scalp_maint_HC{(i+(j-1)*nScalpChans),iSS_High}(indFreq1:indFreq2)));
        
    end
end
dataBipolarScalp_Reconstructed.label = strrep(strcat('_',dataBipolarScalp_Reconstructed.label),'_',' ');
for j = 3%1:6
    Scalp_topoplot(strPaths.Toolboxes,FreqBand_Sclp,freqAxis_sclp,datavector{j},dataBipolarScalp_Reconstructed,4,nChannelPairs_Sclp_HC,1)
    title(['PLV of Hippocampus channel ',strrep(dataBipolar.label{Bip_chans(j)},'_',' '),' to Scalp freq Band', '[', num2str(FreqBand_Sclp(1)),':',num2str(FreqBand_Sclp(2)),']'])
end




%% Scalp HippoCampus
% Plot the spectrum of a pair for every set size

iPair =17; %1:nChannelPairs

figure
for iSS = 1:5
    strChannelNames = dataBipolar_sclp_grid.label(nChannelPairs_Sclp_HC(iPair,:));
    plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{iPair,iSS}),strPlotColors{iSS},'LineWidth',2)
    hold on    xlim([1,30])
    
    ylim([0,0.8])
    set(gca,'XTick',[1,5:5:30])
    title(strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'))
end
%%
% Choose the frequency range you want to focus
FreqBand_Sclp = [8,10];
% FreqBand_Sclp = [15,20];

[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_maint_HC,1)
    for iSS = 1:size(PLV_Pairs_Scalp_maint_HC,2)
        PLV_temp = PLV_Pairs_Scalp_maint_HC{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

% Plot the PLV spectrum of the bipolar channels for every combination with scalp channels for a certain set size
figure
set(gcf,'color','white')
iSS_ToPLot = 5; %Set Size to Calculate the Median
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();
nScalpChans = size(dataBipolar_Sclp.label,1);
for ii = 1:nScalpChans
    for j = 1:num_bipolar_chans
        Scalp_channel = dataBipolar_sclp_grid.label{ii};
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot(ii) = figLayout{indSubplot,2};
        axes(ha(nSubplot(ii))); plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{ii+((j-1)*nScalpChans),iSS_ToPLot}),strPlotColors{j});
        temp = abs(PLV_Pairs_Scalp_maint_HC{ii+((j-1)*nScalpChans),iSS_ToPLot}(indFreq1:indFreq2));
        Data_FOI(ii,j)=median(temp);
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp_HC((ii+((j-1)*nScalpChans)),2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([0 1])
        hold on;
    end
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
end
suptitle(['Scalp-Hippocampus Coupling during Maintenance, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
    end
    
end


Median_Data_FOI = median(Data_FOI');
MedianF{iSS_ToPLot} = Median_Data_FOI;
%%
% Plot the median of Bipolar channels' PLV with the scalp channels for a certain set size
iSS_ToPLot = 5; %Set Size to Plot
figure;

im=imagesc(MedianF{iSS_ToPLot},[0 0.25]);
str = strcat('Median PLV,' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ],', ' n=', num2str(size(dataBipolar_Ret_SS{iSS_ToPLot}.trial,2)),', Maintenance')
title(str)
ax = gca;
ax.XTick = 1:size(dataBipolar_Sclp.label,1);
ax.YTick = 1;
ax.XTickLabel = dataBipolar_Sclp.label;
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

%% Scalp - Hippocampus coupling during encoding
%% Scalp PLV with Hippocampus
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_enc_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_cleanEEG_reref_auricular.mat'],'PLV_Pairs_Scalp_enc_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_encoding_cleanEEG_reref_auricular.mat'],'PLV_Pairs_Scalp_enc_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end




%% Scalp Hippocampus coupling during fixation
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_fix_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\'];
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_cleanEEG_f_30_5_100.mat'],'PLV_Pairs_Scalp_fix_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_cleanEEG_f_30_5_100.mat'],'PLV_Pairs_Scalp_fix_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end


%% Visualize the PLV spectrum for every Hippocampal channel to scalp maint - vs encoding vs fixation
%% Very Important Figures!!

strPlotColors_Hipp = {'r','g','b','k','m','c'};

for j = 1:6
    figure
    set(gcf,'color','white')
    iSS_ToPLot = 4; %Set Size to Calculate the Median
    ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
    figLayout = Get_Figure_Scalp_Layout_Information();
    for i = 1:nScalpChans
        Scalp_Chan = strrep(strcat('_',dataBipolarScalp_Reconstructed.label{i}),'_','');
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_Chan));
        nSubplot(i) = figLayout{indSubplot,2};
        axes(ha(nSubplot(i)));
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{(i+(j-1)*nScalpChans),iSS_ToPLot}),strPlotColors_Hipp{4})
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{(i+(j-1)*nScalpChans),iSS_ToPLot}),strPlotColors_Hipp{2})
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_HC{(i+(j-1)*nScalpChans),iSS_ToPLot}),strPlotColors_Hipp{1})
        Sig1 = abs(PLV_Pairs_Scalp_maint_HC{(i+(j-1)*nScalpChans),iSS_ToPLot});
        Sig2 = abs(PLV_Pairs_Scalp_fix_HC{(i+(j-1)*nScalpChans),iSS_ToPLot});
        indMaxbandSclp_Hc = Difference_Bar(Sig1,Sig2,freqAxis_sclp,[0 1],0.05,'k')
        ylim([0 1]);
        xlim([3 100]);
        elec_pair = strrep(dataBipolarScalp_Reconstructed.label{i},'_',' ');
        ylabel(elec_pair);
    end
    %     set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    suptitle(['PLV of Hippocampus channel ',strrep(dataBipolar.label{Bip_chans(j)},'_',' '),' to Scalp' ])
    for i = 1:size(ha,1)
        
        if ~ismember(i,nSubplot)
            set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
iPair =12; %1:nChannelPairs

figure
for iSS = 1:5
    strChannelNames = dataBipolar_sclp_grid.label(nChannelPairs_Sclp_HC(iPair,:));
    plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{iPair,iSS}),strPlotColors{iSS},'LineWidth',2)
    hold on    xlim([1,30])
    
    ylim([0,0.8])
    set(gca,'XTick',[1,5:5:30])
    title(strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'))
end
clear Data_FOI;
%%
% Choose the frequency range you want to focus
FreqBand_Sclp = [8,10];
% FreqBand_Sclp = [15,20];

[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_enc_HC,1)
    for iSS = 1:size(PLV_Pairs_Scalp_enc_HC,2)
        PLV_temp = PLV_Pairs_Scalp_enc_HC{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

% Plot the PLV spectrum of the bipolar channels for every combination with scalp channels for a certain set size
figure
set(gcf,'color','white')
iSS_ToPLot = 5; %Set Size to Calculate the Median
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();
nScalpChans = size(dataBipolar_Sclp.label,1);

for ii = 1:nScalpChans
    for j = 1:num_bipolar_chans
        Scalp_channel = dataBipolar_sclp_grid.label{ii};
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot(ii) = figLayout{indSubplot,2};
        axes(ha(nSubplot(ii))); plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{ii+((j-1)*nScalpChans),iSS_ToPLot}),strPlotColors{j});
        temp = abs(PLV_Pairs_Scalp_enc_HC{ii+((j-1)*nScalpChans),iSS_ToPLot}(indFreq1:indFreq2));
        Data_FOI(ii,j)=median(temp);
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp_HC((ii+((j-1)*nScalpChans)),2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        ylim([0 1])
        hold on;
    end
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
end
suptitle(['Scalp-Hippocampus Coupling during Encoding, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
    end
    
end

Median_Data_FOI = median(Data_FOI');
MedianF{iSS_ToPLot} = Median_Data_FOI;

% Plot the median of Bipolar channels' PLV with the scalp channels for a certain set size
% iSS_ToPLot = 5; %Set Size to Plot
figure;

im=imagesc(MedianF{iSS_ToPLot});
str = strcat('Median PLV,' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ],', ' n=', num2str(size(dataBipolar_Ret_SS{iSS_ToPLot}.trial,2)),', Encoding')
title(str)
ax = gca;
ax.XTick = 1:size(dataBipolar_Sclp.label,1);
ax.YTick = 1;
ax.XTickLabel = dataBipolar_Sclp.label;
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

%Scalp Hippocampus Maintenance vs Encoding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% iPair =1; %1:nChannelPairs
iSS_to_Plot = [4];
figLayout = Get_Figure_Scalp_Layout_Information();
FreqBand_Sclp = [8,10];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));
clear PLV_AHL23_scalp
iPLV=1;

for j = 1:size(nChannelPairs_Sclp_HC,1)
    Scalp_channel = dataBipolar_sclp_grid.label{nChannelPairs_Sclp_HC(j,2)};
    i = find(strcmpi(Scalp_channel,dataBipolar_sclp_grid.label));
    indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
    nSubplot(i) = figLayout{indSubplot,2};
    if nChannelPairs_Sclp_HC(j,2) == 1
        figure;
        set(gcf,'color','white')
        ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
        lbl = strrep(sprintf('%s',dataBipolar_sclp_grid.label{nChannelPairs_Sclp_HC(j,1)}),'_',' ');
        tlt = [lbl,' Scalp-Hippocampus coupling ','Maintenance vs Encoding ','Set Size ',num2str(Set_Sizes{iSS_to_Plot})];
        suptitle(tlt)
        
    end
    for iSS = iSS_to_Plot
        %         if j<=23
        %             axes(ha(j));
        %         else
        %             axes(ha(j-(floor((j-1)/23)*23)))
        %         end
        axes(ha(nSubplot(i))); plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{j,iSS}),'r','LineWidth',2);
        PLV_band = PLV_Pairs_Scalp_maint_HC{j,iSS}(indFreq1:indFreq2);
        
        ylim([0,1]);
        hold on;
        %         plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{j,iSS}),'g','LineWidth',2);    ylim([0,1]);
        %         hold on;
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp_HC(j,2)};
        elec_pair = strrep(sprintf('%s',elec_pair),'_',' ');
        
        meanPLV_band=mean(abs(PLV_band));
        ylabel(sprintf('%s %2.1f',elec_pair,meanPLV_band));
        if nChannelPairs_Sclp_HC(j,1) == 106
            
            PLV_AHL23_scalp{iPLV,1}=meanPLV_band;
            PLV_AHL23_scalp{iPLV,2} = elec_pair;
            iPLV = iPLV +1;
        end
        
        if j>=47&&j<=69 %% Only for AHL2-3
            if ~isempty(PLV_significance{floor((j-1)/23),i})
                plot(PLV_significance{floor((j-1)/23),i}.plv_sign_band,0.05*ones(length(PLV_significance{floor((j-1)/23),i}.plv_sign_band),1)','Color','magenta','LineWidth',2);
            end
        end
        %     title(sprintf('%s %s %s',strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'),' Encoding,','For every set size'))
    end
    
    set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
    if mod(j,23) == 0
        for i = 1:size(ha,1)
            
            if ~ismember(i,nSubplot)
                set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Calculation of PLV in Scalp-ECoG coupling either during maintenance or
%encoding and comparison between the two periods
%% Scalp PLV with ECoG Grid
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_maint,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];


if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_clean_EEG_f_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_maint','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_maintenance_clean_EEG_f_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_maint','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end


%% Scalp ECoG related
iPair =64; %1:nChannelPairs

figure
for iSS = 1:5
    strChannelNames = dataBipolar_sclp_grid.label(nChannelPairs_Sclp(iPair,:));
    plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{iPair,iSS}),strPlotColors{iSS},'LineWidth',2)
    hold on    xlim([1,30])
    
    ylim([0,0.8])
    set(gca,'XTick',[1,5:5:30])
    title(strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'))
end

% Choose the frequency range you want to focus
FreqBand_Sclp = [6,9];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_maint,1)
    for iSS = 1:size(PLV_Pairs_Scalp_maint,2)
        PLV_temp = PLV_Pairs_Scalp_maint{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end


% Plot the PLV heatmap of  of scalp's Pz or Fz with the grid channels for a certain set size

iSS_ToPLot = 4;
BChan = 15; %Select the Bipolar Channel for every pair of which with the grid chans you want to plot the PLV values

BipChan_to_Plot = find(nChannelPairs_Sclp==BChan);
PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),8,8);

figure;
imagesc(PLV_In_Band_ToPlot,[0 0.7])

% str = strcat(strrep(sprintf('%s',dataBipolar_sclp_grid.label{BChan}),'_',' '), ' [ ', num2str(FreqBand_Sclp(1)),{' '}, ....
%     num2str(FreqBand_Sclp(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]', ' Maintenance')
% title(str)
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

% Plot the PLV spectrum of the Pz,Fz scalp chans for every combination with grid chans for a certain set size
%%
figure
iSS_ToPLot = 5; %Set Size to Calculate the Median

ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
num_chans = 2;
nScalpChans = size(dataBipolar_Sclp.label,1);
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for ii = 1:nGridChans
    for j = 1:num_chans
        axes(ha(plot_order(ii))); plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{ii+((j-1)*nGridChans),iSS_ToPLot}),strPlotColors{j});
        temp = abs(PLV_Pairs_Scalp_maint{ii+((j-1)*nGridChans),iSS_ToPLot}(indFreq1:indFreq2));
        %         Data_FOI(ii,j)=median(temp);
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp((ii+((j-1)*nGridChans)),2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        
        ylim([0 1])
        hold on;
    end
    set(ha(1:nGridChans),'XTickLabel',''); set(ha,'YTickLabel','')
    
    
end
ttl = [sprintf('%s %s','Pz/Fz - ECoG coupling during Maintenance, Set Size  ',num2str(Set_Sizes{iSS_ToPLot}))];
suptitle(ttl)
% Median_Data_FOI = reshape(median(Data_FOI'),8,8);
% MedianF{iSS_ToPLot} = Median_Data_FOI;



%% Scalp - ECoG coupling during encoding
%% Scalp PLV with  ECoG Grid
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_enc,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);
strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];


if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_clean_EEG_f_1_30_all_scalp_chans.mat'],'PLV_Pairs_Scalp_enc','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_encoding_clean_EEG_f_1_30_all_scalp_chans.mat'],'PLV_Pairs_Scalp_enc','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end


%% Scalp - ECoG coupling during fixation
%% Scalp PLV with  ECoG Grid
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_fix,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,TrialInformationTable,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);
strVariableFolder = [Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\'];


if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_clean_EEG_f_30_5_100_all_scalp_chans.mat'],'PLV_Pairs_Scalp_fix','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_fixation_clean_EEG_f_30_5_100_all_scalp_chans.mat'],'PLV_Pairs_Scalp_fix','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end


%% Visualize PLV spectrum for Scalp Channel P3 to grid maint vs encod
%% Important Figures!!

plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for j = 15%[12 13 15 17]%1:size(nChannelPairs_Sclp,1)/nGridChans
    figure
    iSS_ToPLot = 4; %Set Size to Calculate the Median
    strPlotColors_Hipp = {'r','g','b','k','m','c'};
    
    %     ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    nScalpChans = size(dataBipolarScalp_Reconstructed.label,1);
    for i = 31%1:nGridChans
        %         axes(ha(plot_order(i)));
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{(i+(j-1)*nGridChans),iSS_ToPLot}),strPlotColors_Hipp{4})
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{(i+(j-1)*nGridChans),iSS_ToPLot}),strPlotColors_Hipp{2})
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix{(i+(2-1)*nGridChans),iSS_ToPLot}),strPlotColors_Hipp{1})
        
        elec_pair = strrep(dataBipolar_Ret_SS{iSS_ToPLot}.label{i+AHL_chans},'_',' ');
        ylabel(elec_pair);
        ylim([0 1]);
        xlim([3 100])
    end
    ttl = ['PLV Spectrum of channel ',dataBipolarScalp_Reconstructed.label{nChannelPairs_Sclp((i+(j-1)*nGridChans),1)},' to Grid during Maintenance']
    suptitle(ttl)
end

% Plot PLV of  H3 channel to Scalp chans during fix,enc,maint freq 1:1:30
%% Very Important Figure

FreqBand_Sclp = [10,14];
FreqBand_Sclp = [6 9];
FreqBand_Sclp = [10 15];
strPlotColors_Hipp = {'r','g','b','k','m','c'};

[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));
iSS_High = 4;

% Load Percentile for PLV significance
PLV_Rand_Enc_Maint_Variables_Pair_SS = load([Drive_Letter,'Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\PLV_significance_percentile_Scalp_Hippocampus_Maint_Fixation.mat']);
PLV_Rand_Enc_Maint_Variables_Pair_SS = PLV_Rand_Enc_Maint_Variables_Pair_SS.PLV_Rand_Enc_Maint_Variables_Pair_SS;

figure
iSS_ToPLot = 4;
set(gcf,'color','white')
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();
nScalpChans = size(dataBipolarScalp_Reconstructed.label,1);
H3_loc = find(strcmp(dataBipolar_Ret_SS{iSS_ToPLot}.label,'_GL_H2')==1);
H3_loc = H3_loc - AHL_chans;
nPairs = find((nChannelPairs_Sclp(:,2) == H3_loc+AHL_chans+nScalpChans) == 1);
for i = H3_loc
    for j = 1:nScalpChans
        Scalp_channel = strrep(strcat('_',dataBipolarScalp_Reconstructed.label{j}),'_',' ');
        indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
        nSubplot_plv(j) = figLayout{indSubplot,2};
        Vars = PLV_Rand_Enc_Maint_Variables_Pair_SS{iSS_ToPLot}{j};
        Maint_Fix_Diff = abs(PLV_Pairs_Scalp_maint{nPairs(j),iSS_ToPLot})- abs(PLV_Pairs_Scalp_fix{nPairs(j),iSS_ToPLot});
        axes(ha(nSubplot_plv(j)));
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{nPairs(j),iSS_ToPLot}),strPlotColors_Hipp{4});
        datavector_H3_scalp(j) = median(abs(PLV_Pairs_Scalp_maint{nPairs(j),iSS_ToPLot}(indFreq1:indFreq2)));
        
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{nPairs(j),iSS_ToPLot}),strPlotColors_Hipp{2});
        hold on;
        semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix{nPairs(j),iSS_ToPLot}),strPlotColors_Hipp{1});
        indMaxBand_H3_sclp = Difference_Bar(Maint_Fix_Diff,Vars.RandPrc{95},freqAxis_sclp,[0 1],0.05,'k');
        ylabel(Scalp_channel);
        ylim([0 1]);
        xlim([3 100]);
        
    end
end
channel = strrep(dataBipolar_Ret_SS{iSS_ToPLot}.label{H3_loc+AHL_chans},'_',' ');
% suptitle(['PLV spectra - significance for different task periods of channel: ',channel, ' to Scalp'])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot_plv)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
    end
    
end

% Topoplot

Scalp_topoplot(strPaths.Toolboxes,FreqBand_Sclp,freqAxis_sclp,datavector_H3_scalp,dataBipolarScalp_Reconstructed,4,nChannelPairs_Sclp,1)


%%

% PLV map Encoding - P3

FreqBand_Sclp = [8,10];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_enc,1)
    for iSS = 1:size(PLV_Pairs_Scalp_enc,2)
        PLV_temp = PLV_Pairs_Scalp_enc{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

iSS_ToPLot = 4;
BChan = 15; %Select the Channel for every pair of which with the grid chans you want to plot the PLV values

BipChan_to_Plot = [1:64];%find(nChannelPairs_Sclp==BChan);
PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),8,8);

figure;
imagesc(PLV_In_Band_ToPlot,[0.4 1])
colorbar
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

str = strcat(strrep(sprintf('%s',dataBipolarScalp_Reconstructed.label{BChan}),'_',' '), ' [ ', num2str(FreqBand_Sclp(1)),{' '}, ....
    num2str(FreqBand_Sclp(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]', ' Encoding')
title(str)

% PLV map Maint - P3

FreqBand_Sclp = [5,9];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_maint,1)
    for iSS = 1:size(PLV_Pairs_Scalp_maint,2)
        PLV_temp = PLV_Pairs_Scalp_maint{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

iSS_ToPLot = 4;
BChan = 15; %Select the Channel for every pair of which with the grid chans you want to plot the PLV values

BipChan_to_Plot = [1:64]%find(nChannelPairs_Sclp==BChan);
PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),8,8);

figure;
imagesc(PLV_In_Band_ToPlot,[0.4 1])
colorbar
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

str = strcat(strrep(sprintf('%s',dataBipolarScalp_Reconstructed.label{BChan}),'_',' '), ' [ ', num2str(FreqBand_Sclp(1)),{' '}, ....
    num2str(FreqBand_Sclp(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]', ' Maintenance')
title(str)

%%
% Visualizations
iPair =1; %1:nChannelPairs

figure
for iSS = 1:5
    strChannelNames = dataBipolar_sclp_grid.label(nChannelPairs_Sclp(iPair,:));
    plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{iPair,iSS}),strPlotColors{iSS},'LineWidth',2)
    hold on    xlim([1,30])
    
    ylim([0,1])
    set(gca,'XTick',[1,5:5:30])
    title(sprintf('%s %s %s',strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'),' Encoding,','For every set size'))
end
%%
% Choose the frequency range you want to focus
FreqBand_Sclp = [5,9];
[~,indFreq1] = min(abs(freqAxis_sclp-FreqBand_Sclp(1)));
[~,indFreq2] = min(abs(freqAxis_sclp-FreqBand_Sclp(2)));

% Plot
PLV_In_Band_SS = [];
for iPair = 1:size(PLV_Pairs_Scalp_maint,1)
    for iSS = 1:size(PLV_Pairs_Scalp_maint,2)
        PLV_temp = PLV_Pairs_Scalp_maint{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end


% Plot the PLV heatmap of  of scalp's Pz or Fz with the grid channels for a certain set size

iSS_ToPLot = 4;
BChan = 15; %Select the Channel for every pair of which with the grid chans you want to plot the PLV values

BipChan_to_Plot = find(nChannelPairs_Sclp==BChan);
PLV_In_Band_ToPlot = PLV_In_Band_SS(:,iSS_ToPLot);
PLV_In_Band_ToPlot = reshape(PLV_In_Band_ToPlot(BipChan_to_Plot),8,8);

figure;
imagesc(PLV_In_Band_ToPlot,[0 0.7])

% str = strcat(strrep(sprintf('%s',dataBipolar_sclp_grid.label{BChan}),'_',' '), ' [ ', num2str(FreqBand_Sclp(1)),{' '}, ....
%     num2str(FreqBand_Sclp(2)) ,' ] Hz ' , ' Set Size',{' ['}, num2str(Set_Sizes{iSS_ToPLot}), ' ]', ' Encoding')
title(str)
ax = gca;
ax.XTick = 1:8;
ax.YTick = 1:8;
ax.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax.YTickLabel = 1:8;
ax.YDir = 'normal'

colorbar
colormap jet

% Plot the PLV spectrum of the Pz,Fz scalp chans for every combination with grid chans for a certain set size
%%
figure
iSS_ToPLot = 5; %Set Size to Calculate the Median

ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
num_chans = 2;
nScalpChans = size(dataBipolar_Sclp.label,1);
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for ii = 1:nGridChans
    for j = 1:num_chans
        axes(ha(plot_order(ii))); plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{ii+((j-1)*nGridChans),iSS_ToPLot}),strPlotColors{j});
        temp = abs(PLV_Pairs_Scalp_enc{ii+((j-1)*nGridChans),iSS_ToPLot}(indFreq1:indFreq2));
        %         Data_FOI(ii,j)=median(temp);
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp((ii+((j-1)*nGridChans)),2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        
        ylim([0 1])
        hold on;
    end
    set(ha(1:nGridChans),'XTickLabel',''); set(ha,'YTickLabel','')
    
    
end
ttl = [sprintf('%s %s','Pz/Fz - ECoG coupling during Encoding, Set Size  ',num2str(Set_Sizes{iSS_ToPLot}))];
suptitle(ttl)
%% Scalp-ECoG coupling maintenance vs encoding
iPair =1; %1:nChannelPairs
iSS_to_Plot = [4];
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for j = 1:size(nChannelPairs_Sclp,1)
    if j == 1||j == size(nChannelPairs_Sclp,1)-nGridChans+1
        figure;
        ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01]);
        tlt = [dataBipolar_Sclp.label{nChannelPairs_Sclp(j,1)},' Scalp-ECoG coupling ','Maintenance vs Encoding ','Set Size ',num2str(Set_Sizes{iSS_to_Plot})];
        suptitle(tlt)
    end
    for iSS = iSS_to_Plot
        strChannelNames = dataBipolar_sclp_grid.label(nChannelPairs_Sclp(iPair,:));
        if j<=64
            axes(ha(plot_order(j)));
        else
            axes(ha(plot_order(j-64)))
        end
        plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{j,iSS}),'r','LineWidth',2);  ylim([0,1]);
        hold on;
        plot(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{j,iSS}),'g','LineWidth',2);    ylim([0,1]);
        elec_pair = dataBipolar_sclp_grid.label{nChannelPairs_Sclp(j,2)};
        ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
        
        %     title(sprintf('%s %s %s',strrep(sprintf('%s - %s',strChannelNames{1},strChannelNames{2}),'_','\_'),' Encoding,','For every set size'))
    end
    set(ha(1:nGridChans),'XTickLabel',''); set(ha,'YTickLabel','')
    
end

%% Statistics
% Check the script: PLV_Significance_calculation
% [ stat ] = Get_PLV_significance_P42DS(PLV_Pairs_retention, [0.5:0.5:30], freqAxis, [8 12], [0.5 30], 4 )

