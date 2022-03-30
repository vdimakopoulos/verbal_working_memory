
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
Session1                   =      'F:\Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_01_Part_01';
Session2                   =      'F:\Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_02_Part_01';
data_s1                    =      load(Session1);
data_s2                    =      load(Session2);
TrialInformationTable      =      [data_s1.TrialInformationTable;data_s2.TrialInformationTable] ;
data_s1                    =      data_s1.dataMacro;
data_s2                    =      data_s2.dataMacro;


%% Load scalp data for 2 sessions
Sclp_Session1                    =      'F:\Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_01_Part_01';
Sclp_Session2                    =      'F:\Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_02_Part_01';
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
    strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\';
elseif Grid_Hippocampus_coupling
     coupling_type = 'Grid-Hippocampus';
     strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\';
elseif Scalp_Grid_coupling
    coupling_type = 'Scalp-Grid';
     strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp-Grid coupling\';
end
    
%%
IncorrectTrials = TrialInformationTable(find(~TrialInformationTable.Correct),:);
cfg = [];
cfg.trials = find(~TrialInformationTable.Correct);
dataBipolar = ft_selectdata(cfg,dataBipolar);
dataBipolarScalp = ft_selectdata(cfg,dataBipolarScalp);
TrialInformationTable = TrialInformationTable(find(~TrialInformationTable.Correct),:);
%% Reref based on the average of the signals
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
cfg               =   [];
cfg.reref         =  'yes' 
cfg.refchannel    =  [3 4] %white matter contacts referencing
cfg.refmethod     =  'avg'%'avg'
% cfg.channel       = [1:8] % NEW LINES
dataBipolar       = ft_preprocessing(cfg,dataBipolar);

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

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);

%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar,IncorrectTrials);
[dataBipolar_ScalpSS,TrialInformationTableScalp_SS,nTrialList_TT_SS]= Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolarScalp,IncorrectTrials);

%% Extract 2 seconds of maintenance
if strcmp(TaskPeriod,'maint')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];

        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    %%
    %Power spectrum of maintenance
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;

    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
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
    
elseif strcmp(TaskPeriod,'encod')
    %% Extract 2 seconds of encoding latency
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;

    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    fr_enc =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
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
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;

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
%% PSD

figure;
semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(67,:)),'r','LineWidth',3);
hold on;
semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(67,:)),'g','LineWidth',3);
semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(67,:)),'k','LineWidth',3);
xlim([4 100]);
ylim([-10 40]);
title(strrep(fr_maint.label{67},'_',' '));

figure;
semilogx(fr_maint.freq,10*log10(fr_maint.powspctrm(83,:)),'r','LineWidth',3);
hold on;
semilogx(fr_enc.freq,10*log10(fr_enc.powspctrm(83,:)),'g','LineWidth',3);
semilogx(fr_fix.freq,10*log10(fr_fix.powspctrm(83,:)),'k','LineWidth',3);
xlim([4 100]);
ylim([-20 20]);
title(strrep(fr_maint.label{83},'_',' '));
%% Channel pairs to run
%To be configured for every combination of electrode pairs you have
PLV_depth_electrodes = 0; % Boolean Flag for calculating the PLV between the bipolar channels in the depth electrodes
clear nChannelPairs;

if PLV_depth_electrodes == 1
    
    nChannelPairs = [81 83;81 84; 82 83; 82 84];
    
else
    nGridChans = 64;
    AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
    PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
    nChannelPairs = [[ones(nGridChans,1),(1:nGridChans)'+AHL_chans];[nGridChans+AHL_chans+ones(nGridChans,1),(1:nGridChans)'+AHL_chans]];
    for i=1:2:length(chans_to_be_bipolar)
        nChannelPairs = [nChannelPairs;[round(i/2)+size(montage.labelold,1)+zeros(nGridChans,1),(1:nGridChans)'+ AHL_chans]];
    end

end
%% Time frequency granger for all hipp-grid pairs
indBipChans = find(contains(strChannelNameList(nChannelPairs(:,1)),'-'));
hippGrid_ChanPairs = nChannelPairs(indBipChans,:);
freq = [1:100];
channel_cmb = {'-','_'};
strChannelNameList = dataBipolar_SS{1}.label;
[Granger_TFR_ss,Granger_struct] = getTFR_granger_all_hippGridPairs(freq,strChannelNameList,hippGrid_ChanPairs,channel_cmb,dataBipolar_SS)


%% Visualize the results from TFR granger
Bip_chans_to_plot = 83%unique(hippGrid_ChanPairs(:,1));
grangerTimeAxis = [-6:0.25:2];
grangerFreqAxis = [1:100];
figure;
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
iSS =5;
clim = [-0.2 0.2];%[-0.05 0.05];
for i = 1:numel(Bip_chans_to_plot)
    indBipChan = find(hippGrid_ChanPairs==Bip_chans_to_plot(i))
    figure('units','normalized','outerposition',[0 0 1 1])
    set(gcf,'Color','white')
    ha = tight_subplot(8,8,[.06 .03],[.1 .01],[.04 .03])
    for j = 1:numel(indBipChan)
        chanPair = indBipChan(j);
        axes(ha(plot_order(j)));
        contourf(grangerTimeAxis,grangerFreqAxis,-Granger_TFR_ss{chanPair,iSS},100,'LineColor','none');
        set(gca,'clim',clim,'yscale','log');
        set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid
        colormap(bluewhitered(128));
        set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
        ylim([4 100]);
        ylim([5 30])
        colorbar;
        ylabel('Frequency (Hz)');
        xlabel('Time (s)');
        set(gca,'FontSize',16,'box','off');
        strTitle = strrep(strChannelNameList(hippGrid_ChanPairs(chanPair,:)),'_','');
        title(sprintf('%s',strTitle{2}));
        
    end
    suptitle(sprintf('%s - Grid Granger',strTitle{1}));
end
%%
tStart = tic;
for iPair = 257:320%314%129:size(nChannelPairs,1)
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
%         cfg.foi         = 0.5:0.5:30;
        if GammaFreqPLV
            cfg.foi         = 0.5:10:120;
        else
%             cfg.foi         = 0.5:0.5:30;
            cfg.foi      = [1:1:30 30:5:100];
%             cfg.foi         = 4:1:100;
%               cfg.foi         = 30:5:100;

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

 if Grid_Hippocampus_coupling
        strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Grid_Hippocampus coupling\'
        mkdir(strVariableFolder);
        
        if ~GammaFreqPLV
            if strcmp(TaskPeriod,'maint')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'encod')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_f_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'retr')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_retrieval.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
            elseif strcmp(TaskPeriod,'fix')
                save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_f_1_30_5_100.mat'],'PLV_Pairs','nChannelPairs','freqAxis','-v7.3')
                
            end
        end
 end
 
 %% Plot PLV results on a single channel pair
 strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Grid_Hippocampus coupling\'

 strFilenameMaint = [strVariableFolder,'Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_f_1_1_30_5_100.mat'];
 strFilenameEnc   = [strVariableFolder,'Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_f_1_30_5_100.mat'];
 strFilenameFix   = [strVariableFolder,'Patient_42_Grid-Hippocampus_PLV_Pairs_fixation_f_1_30_5_100.mat'];
 
 % Pairs to run analysis for
 strChannelNameList = dataBipolar.label;
 nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
 ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','_AHL2-_AHL3')));
 ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_GL_C2')));
 nPairs_ToRun = intersect(nPairs_ToRun,ind);
 nPairs_ToRun = intersect(nPairs_ToRun,ind2);
 
 iSS = 4;
 
 figure;
 load(strFilenameFix);
 semilogx(freqAxis,abs(PLV_Pairs{nPairs_ToRun,iSS}),'k','LineWidth',3);
 hold on;
 load(strFilenameEnc);
 semilogx(freqAxis,abs(PLV_Pairs{nPairs_ToRun,iSS}),'g','LineWidth',3);
 load(strFilenameMaint);
 semilogx(freqAxis,abs(PLV_Pairs{nPairs_ToRun,iSS}),'r','LineWidth',3);
 strTitle = ['Hipp - ',strrep(strChannelNameList{nChannelPairs(nPairs_ToRun,2)},'_',' ')];
 title(strTitle);

 ylim([0 1]);
 xlim([4 100]);
%% Plot the heatmap of the PLV for a certain bipolar channel and for a certain set size

FreqBand = [16,29];  

[~,indFreq1] = min(abs(freqAxis-FreqBand(1)));
[~,indFreq2] = min(abs(freqAxis-FreqBand(2)));

% Plot
PLV_In_Band_SS = [];
iSS_to_Plot =4; %[1:5]
for iPair = 257:320%129:size(nChannelPairs,1) %193:256%1:384%
    for iSS = iSS_to_Plot
        
        PLV_temp = PLV_Pairs{iPair,iSS}(indFreq1:indFreq2);
        PLV_temp = abs(PLV_temp);
        PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
    end
end

%Important
iSS_ToPLot = 4;
BChan = 83; %Select the Bipolar Channel for every pair of which with the grid chans you want to plot the PLV values
if PLV_depth_electrodes == 1
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
    imagesc(PLV_In_Band_ToPlot,[0 0.7])
    
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
    colormap jet
    strtitle = ['[',num2str(FreqBand(1)),num2str(FreqBand(2)),']',' Hz'];
    title(strtitle);
end


%% Scalp Hippocampus

%Scalp data 
cfg               = []
cfg.channel       = {'all', '-_Submm', '-_Submp'}
cfg.reref         =  'yes' 
cfg.refchannel    =   {'_A1' '_A2'}
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

%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
[dataBipolar_sclp_grid_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_sclp_grid,IncorrectTrials);

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
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 1:23;
    fr_maint_sclp = ft_freqanalysis(cfg,dataBipolar_scalp_Ret_SS{4});
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
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 1:23;
    fr_enc_sclp =ft_freqanalysis(cfg,dataBipolar_scalp_Ret_SS{4});
elseif strcmp(TaskPeriod,'fix')
    nSet_Size = size(dataBipolar_sclp_grid_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar.fsample];
        dataBipolar_scalp_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_sclp_grid_SS{iSS});
    end
    %-------------------------------------------------------------------------
    %Power spectrum of scalp channels
    cfg = []; cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
    cfg.foi         = 0.5:0.5:100;
    cfg.tapsmofrq   = 2;
    cfg.channel = 1:23;
    fr_fix_sclp =ft_freqanalysis(cfg,dataBipolar_scalp_Ret_SS{4});
    
end

%% Plot PSD of a scalp channel for every task period
freqAxis = 0.5:0.5:100;
figure;
semilogx(fr_fix.freq,10*log10(fr_fix_sclp.powspctrm(15,:)),'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(fr_enc_sclp.powspctrm(15,:)),'g','LineWidth',3);
semilogx(freqAxis,10*log10(fr_maint_sclp.powspctrm(15,:)),'r','LineWidth',3);
xlim([4 100]);
ylim([-10 50]);
title(strrep(fr_maint_sclp.label{15},'_',' '))
 
%For every scalp channel

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
    semilogx(fr_fix_sclp.freq,10*log10(fr_fix_sclp.powspctrm(ii,:)),'k','LineWidth',3);
    hold on;
    semilogx(fr_enc_sclp.freq,10*log10(fr_enc_sclp.powspctrm(ii,:)),'g','LineWidth',3);
    semilogx(fr_maint_sclp.freq,10*log10(fr_maint_sclp.powspctrm(ii,:)),'r','LineWidth',3);
    ylim([-20 50]);
    xlim([4 100])
    elec_pair = fr_enc_sclp.label{ii};
    ylabel(elec_pair);
end
set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% set(ha(1:size(ha,1)),'FontSize',12);
% suptitle(['Scalp Channels Power Spectrum Maintenance vs Encoding, Set Size ', num2str(Set_Sizes{iSS_ToPLot})])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot_psd)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
    end
    
end


%% Scalp PLV with  Hippocampus
%maintenance
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_maint_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp_Hippocampus coupling\';
mkdir(strVariableFolder)
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_maint_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_maintenance_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_maint_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end    

% encoding
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_enc_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp_Hippocampus coupling\';
mkdir(strVariableFolder)
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_encoding_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_enc_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end  


% fixation
PLV_Pair_flag = 1; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: fixation
[PLV_Pairs_Scalp_fix_HC,nChannelPairs_Sclp_HC,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);


strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp_Hippocampus coupling\';
mkdir(strVariableFolder)
if strcmp(coupling_type,'Scalp-Hippocampus')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_fix_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Hippocampus','_PLV_Pairs_fixation_cleanEEG_f_1_1_30_5_100.mat'],'PLV_Pairs_Scalp_fix_HC','nChannelPairs_Sclp_HC','freqAxis_sclp','-v7.3');
end  


%% Plot results for a single hipp-scalp channel pair 

 strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp_Hippocampus coupling\';
 strFilenameMaint = [strVariableFolder,'Patient_42_Scalp-Hippocampus_PLV_Pairs_maintenance_cleanEEG_f_1_1_30_5_100.mat'];
 strFilenameEnc   = [strVariableFolder,'Patient_42_Scalp-Hippocampus_PLV_Pairs_encoding_cleanEEG_f_1_1_30_5_100.mat'];
 strFilenameFix   = [strVariableFolder,'Patient_42_Scalp-Hippocampus_PLV_Pairs_fixation_cleanEEG_f_1_1_30_5_100.mat'];
 
 % Pairs to run analysis for
 strChannelNameList = [dataBipolar_Sclp.label; dataBipolar.label'];
 nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
 ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Sclp_HC(:,1))','_AHL2-_AHL3')));
 ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Sclp_HC(:,2))',' P3')));
 nPairs_ToRun = intersect(nPairs_ToRun,ind);
 nPairs_ToRun = intersect(nPairs_ToRun,ind2);
 
 iSS = 4;
 
 figure;
 load(strFilenameFix);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix_HC{nPairs_ToRun,iSS}),'k','LineWidth',3);
 hold on;
 load(strFilenameEnc);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc_HC{nPairs_ToRun,iSS}),'g','LineWidth',3);
 load(strFilenameMaint);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint_HC{nPairs_ToRun,iSS}),'r','LineWidth',3);
 strTitle = ['Hipp - ',strrep(strChannelNameList{nChannelPairs_Sclp_HC(nPairs_ToRun,2)},'_',' ')];
 title(strTitle);

 ylim([0 1]);
 xlim([4 100]);
 
 
 %% Scalp-Grid
%maintenance
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 2; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_maint,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp-Grid coupling\';
mkdir(strVariableFolder)

if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_maint','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_maintenance_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_maint','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end  

%encoding
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 1; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_enc,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp-Grid coupling\';
mkdir(strVariableFolder)

if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_enc','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_encoding_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_enc','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end  

%fixation
PLV_Pair_flag = 2; % Value 1 for Hippocampus Coupling, Value 2 for ECoG grid coupling
Task_Period_Flag = 3; %Value for choosing the latency between different task periods, 1:encoding, 2:maintenance, 3: retrieval
[PLV_Pairs_Scalp_fix,nChannelPairs_Sclp,freqAxis_sclp] = calculate_PLV_scalp_HL(dataBipolar,data_scalp_s1,data_scalp_s2,IncorrectTrials,montage,chans_to_be_bipolar,PLV_Pair_flag,Task_Period_Flag);

strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp-Grid coupling\';
mkdir(strVariableFolder)

if strcmp(coupling_type,'Scalp-Grid')
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_fix','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
else
    save([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_','Scalp-Grid','_PLV_Pairs_fixation_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref.mat'],'PLV_Pairs_Scalp_fix','nChannelPairs_Sclp','freqAxis_sclp','-v7.3');
end  


%% Plot results for a single grid-scalp channel pair 

strVariableFolder = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Incorrect_Trials\Scalp-Grid coupling\';
strFilenameMaint = [strVariableFolder,'Patient_42_Scalp-Grid_PLV_Pairs_maintenance_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref'];
strFilenameEnc   = [strVariableFolder,'Patient_42_Scalp-Grid_PLV_Pairs_encoding_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref'];
strFilenameFix   = [strVariableFolder,'Patient_42_Scalp-Grid_PLV_Pairs_fixation_clean_EEG_f_1_1_30_5_100_all_scalp_chans_reref'];
 
 % Pairs to run analysis for
 strChannelNameList = [dataBipolar_Sclp.label; dataBipolar.label'];
 nPairs_ToRun = 1:size(nChannelPairs_Sclp,1); % all pairs for real values
 ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Sclp(:,1))',' P3')));
 ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs_Sclp(:,2))','_GL_D6')));
 nPairs_ToRun = intersect(nPairs_ToRun,ind);
 nPairs_ToRun = intersect(nPairs_ToRun,ind2);
 
 iSS = 4;
 
 figure;
 load(strFilenameFix);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_fix{nPairs_ToRun,iSS}),'k','LineWidth',3);
 hold on;
 load(strFilenameEnc);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_enc{nPairs_ToRun,iSS}),'g','LineWidth',3);
 load(strFilenameMaint);
 semilogx(freqAxis_sclp,abs(PLV_Pairs_Scalp_maint{nPairs_ToRun,iSS}),'r','LineWidth',3);
 strTitle = ['EEG',strrep(strChannelNameList{nChannelPairs_Sclp(nPairs_ToRun,1)},'_',' '), ' - ',strrep(strChannelNameList{nChannelPairs_Sclp(nPairs_ToRun,2)},'_',' ')];
 title(strTitle);

 ylim([0 1]);
 xlim([4 100]);
 