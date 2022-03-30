%% Paths
% Folder that contains the folder 'Sternberg Task'
strPaths.Main = 'E:\MATLAB Codes\';
strPaths.GeneralFunctions = [strPaths.Main,'General Functions\'];

% Main folder of task and subfolders
strPaths.Sternberg                                        = [strPaths.Main,'Sternberg Task\'];
strPaths.EEGLABDataClean.Main                             = [strPaths.Sternberg,'EEGLAB Data Clean 1\'];
strPaths.EEGLABDataClean.Patients.Main                    = [strPaths.EEGLABDataClean.Main,'Patients\'];
strPaths.EEGLABDataClean.Patients.Macro                   = [strPaths.EEGLABDataClean.Patients.Main,'Macro Data\'];
strPaths.EEGLABDataClean.Patients.Macro_Bipolar           = [strPaths.EEGLABDataClean.Patients.Main,'Macro Bipolar Data\'];
strPaths.EEGLABDataClean.Patients.Scalp                   = [strPaths.EEGLABDataClean.Patients.Main,'Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Macro_Bipolar_Scalp     = [strPaths.EEGLABDataClean.Patients.Main,'Macro Bipolar Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Micro                   = [strPaths.EEGLABDataClean.Patients.Main,'Micro Data\'];
strPaths.EEGLABDataClean.Patients.First_3_Bipolar_Scalp   = [strPaths.EEGLABDataClean.Patients.Main,'Macro First 3 Bipolar Scalp Data\'];
strPaths.EEGLABDataClean.Patients.Params                  = [strPaths.EEGLABDataClean.Patients.Main,'Artifact Rejection Params\'];

% Recording information
strPaths.PatientRecordingInformation = [strPaths.Sternberg,'Patient Information\'];

% Results
strPaths.Results = [strPaths.Sternberg,'Results\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'E:\MATLAB Codes\Toolboxes\fieldtrip-20170925\';
% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = 'E:\MATLAB Codes\Toolboxes\eeglab14_1_1b\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
%  addpath(genpath(strPaths.Sternberg))
% addpath(genpath(strPaths.GeneralFunctions))
addpath(genpath(strPaths.Toolboxes.FieldTrip))
% rmpath(genpath(strPaths.GeneralFunctions))

% Remove EEGLAB from path
rmpath(genpath(strPaths.Toolboxes.EEGLAB))


% Plot colors
strPlotColors = {'b','g','r','c','k','m'};
ft_defaults

% % Resolve legend errors
%  restoredefaultpath;
% rehash toolboxcache

set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

%% Load data for 2 sessions
Session1                   =      'E:\MATLAB Codes\Sternberg Task\NLX Data Extraction\Sessions\Macro Data\Macro_Data_Sessions_Patient_42_Session_01_Part_01';
Session2                   =      'E:\MATLAB Codes\Sternberg Task\NLX Data Extraction\Sessions\Macro Data\Macro_Data_Sessions_Patient_42_Session_02_Part_01';
data_s1                    =      load(Session1);
data_s2                    =      load(Session2);
TrialInformationTable      =      [data_s1.TrialInformationTable;data_s2.TrialInformationTable] ;
data_s1                    =      data_s1.dataMacro;
data_s2                    =      data_s2.dataMacro;


%% Load scalp data for 2 sessions
Sclp_Session1                    =      'E:\MATLAB Codes\Sternberg Task\NLX Data Extraction\Sessions\Scalp Data\Scalp_Data_Sessions_Patient_42_Session_01_Part_01';
Sclp_Session2                    =      'E:\MATLAB Codes\Sternberg Task\NLX Data Extraction\Sessions\Scalp Data\Scalp_Data_Sessions_Patient_42_Session_02_Part_01';
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
    
%% Reref based on the average of the signals
% [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
cfg               =   [];
cfg.reref         =  'yes' 
cfg.refchannel    =  [3,4] %white matter contacts referencing
cfg.refmethod     =  'avg'
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
cfg.refmethod = 'bipolar'
cfg.refchannel = [82 85]
cfg.montage = montage;
dataBipolar = ft_preprocessing(cfg,dataBipolar);

Bip_chans = (find(contains(dataBipolar.label, '-')==1));

%% Downsample to 500 Hz
cfg = [];
cfg.resamplefs = 500;
dataBipolar = ft_resampledata(cfg,dataBipolar);
dataBipolarScalp = ft_resampledata(cfg,dataBipolarScalp);

%% Append scalp with ECoG and hippocampus
if Scalp_Grid_coupling || Scalp_Hippocampus_coupling

    cfg = [];
    cfg.keepsampleinfo = 'no' 
    dataBipolar = ft_appenddata(cfg,dataBipolarSclp,dataBipolar);
end



%% Select only correct trials
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
TaskPeriod = 'retr';
if strcmp(TaskPeriod,'fix')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-6,5-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
     %Power spectrum fixation
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
%     cfg.foi         = 0.5:0.5:30;
    cfg.foi         = 0.5:0.5:120;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    
    if cfg.foi(end)== 120
        fr_fix_high_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
        fr_fix_low_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
        
    else
       fr_fix_high =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
       fr_fix_low =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    
    end
    
end


if strcmp(TaskPeriod,'encod')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
     %Power spectrum encoding
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
%     cfg.foi         = 0.5:0.5:30;
    cfg.foi         = 0.5:0.5:120;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    
    if cfg.foi(end)== 120
       fr_enc_high_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
       fr_enc_low_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    else
       fr_enc_high =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
       fr_enc_low =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    end
    
end


if strcmp(TaskPeriod,'maint')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    %Power spectrum maintenance
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
%     cfg.foi         = 0.5:0.5:30;
    cfg.foi         = 0.5:0.5:120;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    if cfg.foi(end)== 120
        fr_maint_high_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
        fr_maint_low_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    else
        fr_maint_high =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
        fr_maint_low =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    end
end
if strcmp(TaskPeriod,'retr')
    nSet_Size = size(dataBipolar_SS,2);
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [0,2-1/dataBipolar_SS{iSS}.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
     %Power spectrum retrieval
    cfg = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'pow';
%     cfg.foi         = 0.5:0.5:30;
    cfg.foi         = 0.5:0.5:120;
    cfg.tapsmofrq   = 2;
    cfg.channel = 'all';
    if cfg.foi(end)== 120
        fr_retr_high_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
        fr_retr_low_gamma =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    else
        fr_retr_high =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
        fr_retr_low =ft_freqanalysis(cfg,dataBipolar_Ret_SS{1});
    end
end



%% Plot the power spectrum for high and low workload
%Scalp-channels:

nGridChans = 64;
AHL_chans = length(find(contains(montage.labelold, 'AHL')==1));
PHL_chans = length(find(contains(montage.labelold, 'PHL')==1));
nScalpChans = size(dataBipolarSclp.label,1);

if length(dataBipolar_Ret_SS{iSS_Low}.label)<= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end

strPlotColors = {'r','g','b','k',};
boolLinearScale = 0;

% Power spectrum for Fixation,encoding,maintenance high, maintenance low
fig = figure
set(gcf,'color','white')
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();
for i = 1:nScalpChans
    Scalp_channel = dataBipolar.label{i};
    indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
    nSubplot_psd(i) = figLayout{indSubplot,2};
    axes(ha(nSubplot_psd(i)));
    if ~boolLinearScale
        plot(fr_fix_high.freq,10*log10(fr_fix_high.powspctrm(i,:)),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high.freq,10*log10(fr_enc_high.powspctrm(i,:)),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low.freq,10*log10(fr_maint_low.powspctrm(i,:)),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high.freq,10*log10(fr_maint_high.powspctrm(i,:)),strPlotColors{4},'LineWidth',2);
       
        ylim([-15 30]);
    else
        plot(fr_fix_high.freq,fr_fix_high.powspctrm(i,:),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high.freq,fr_enc_high.powspctrm(i,:),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low.freq,fr_maint_low.powspctrm(i,:),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high.freq,fr_maint_high.powspctrm(i,:),strPlotColors{4},'LineWidth',2);
        ylim([-1000 2000]);
    end
    elec_pair = fr_fix_high.label{i};
    ylabel(elec_pair);
end

set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% set(ha(1:size(ha,1)),'FontSize',18);
suptitle(['Power spectral density in different task periods Scalp Channels'])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot_psd)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
        
    end
    
end
% legend('Fixation','Encoding','Maintenance High Workload','Maintenance Low Workload')


%% Power spectrum for high frequencies (gamma)
% Grid channels
Grid_Ch_Range = [AHL_chans+1:nGridChans+AHL_chans];
figure
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
if length(dataBipolar_Ret_SS{iSS_Low}.label)>= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end 
clear indMaxBand indMaxBand2;
boolLinearScale = 1;
for i = Grid_Ch_Range
    pos = find(Grid_Ch_Range == i);
    axes(ha(plot_order(pos)));
    if ~boolLinearScale
        plot(fr_fix_high_gamma.freq,10*log10(fr_fix_high_gamma.powspctrm(i,:)),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high_gamma.freq,10*log10(fr_enc_high_gamma.powspctrm(i,:)),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low_gamma.freq,10*log10(fr_maint_low_gamma.powspctrm(i,:)),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high_gamma.freq,10*log10(fr_maint_high_gamma.powspctrm(i,:)),strPlotColors{4},'LineWidth',2);
        ylim([-15 40]);
        xlim([2 100]);

    else
        plot(fr_fix_high_gamma.freq,fr_fix_high_gamma.powspctrm(i,:),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high_gamma.freq,fr_enc_high_gamma.powspctrm(i,:),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low_gamma.freq,fr_maint_low_gamma.powspctrm(i,:),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high_gamma.freq,fr_maint_high_gamma.powspctrm(i,:),strPlotColors{4},'LineWidth',2);
%          ylim([0 8]);
        ylim([0 200])
        xlim([2 100]);
    end
    if ~boolLinearScale
        Signal_d1 = fr_maint_high_gamma.powspctrm(i,:);
        Signal_d2 = fr_enc_high_gamma.powspctrm(i,:);
        yShift1 = -15;
        yShift2 = -14;
        color1  = 'm';
        color2  = 'c';
        indMaxBand = Difference_Bar(Signal_d1,Signal_d2,fr_maint_high_gamma.freq,[-15 40],yShift1,color1);
        indMaxBand2 = Difference_Bar(Signal_d2,Signal_d1,fr_maint_high_gamma.freq,[-15 40],yShift2,color2);
                
    else
        Signal_d1 = fr_maint_high_gamma.powspctrm(i,:);
        Signal_d2 = fr_enc_high_gamma.powspctrm(i,:);
        yShift1 = 0.15;
        yShift2 = 0.24;
        color1  = 'm';
        color2  = 'c';
        indMaxBand = Difference_Bar(Signal_d1,Signal_d2,fr_maint_high_gamma.freq,[0 200],yShift1,color1);
        indMaxBand2 = Difference_Bar(Signal_d2,Signal_d1,fr_maint_high_gamma.freq,[0 200],yShift2,color2);
              
    end
    
    elec_pair = fr_retr_high.label{i};
    ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
end
set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% set(ha(1:size(ha,1)),'FontSize',12);



suptitle(['Power spectral density in different task periods -  Grid Channels - Gamma Frequency'])

%% Power spectrum for low frequencies 
% Grid channels
Grid_Ch_Range = [AHL_chans+1:nGridChans+AHL_chans];
figure
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
if length(dataBipolar_Ret_SS{iSS_Low}.label)>= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end

boolLinearScale = 1;
for i = Grid_Ch_Range
    pos = find(Grid_Ch_Range == i);
    axes(ha(plot_order(pos)));
    if ~boolLinearScale
        plot(fr_fix_high.freq,10*log10(fr_fix_high.powspctrm(i,:)),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high.freq,10*log10(fr_enc_high.powspctrm(i,:)),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low.freq,10*log10(fr_maint_low.powspctrm(i,:)),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high.freq,10*log10(fr_maint_high.powspctrm(i,:)),strPlotColors{4},'LineWidth',2);
        ylim([-10 40]);
        
    else
        plot(fr_fix_high.freq,fr_fix_high.powspctrm(i,:),strPlotColors{1},'LineWidth',2);
        hold on;
        plot(fr_enc_high.freq,fr_enc_high.powspctrm(i,:),strPlotColors{2},'LineWidth',2);
        hold on;
        plot(fr_maint_low.freq,fr_maint_low.powspctrm(i,:),strPlotColors{3},'LineWidth',2);
        hold on;
        plot(fr_maint_high.freq,fr_maint_high.powspctrm(i,:),strPlotColors{4},'LineWidth',2);
        ylim([-1000 2000]);
    end
    if ~boolLinearScale
        Signal_d1 = 10*log10(fr_maint_high.powspctrm(i,:));
        Signal_d2 = 10*log10(fr_enc_high.powspctrm(i,:));
        yShift1 = -10;
        yShift2 = -8;
        color1  = 'm';
        color2  = 'c';
        indMaxBand = Difference_Bar(Signal_d1,Signal_d2,fr_maint_high.freq,[-10 40],yShift1,color1);
        indMaxBand2 = Difference_Bar(Signal_d2,Signal_d1,fr_maint_high.freq,[-10 40],yShift2,color2);
                
    else
        Signal_d1 = fr_maint_high.powspctrm(i,:);
        Signal_d2 = fr_enc_high.powspctrm(i,:);
        yShift1 = -1000;
        yShift2 = -800;
        color1  = 'm';
        color2  = 'c';
        indMaxBand = Difference_Bar(Signal_d1,Signal_d2,fr_maint_high.freq,[-1000 2000],yShift1,color1);
        indMaxBand2 = Difference_Bar(Signal_d2,Signal_d1,fr_maint_high.freq,[-1000 2000],yShift2,color2);
    end
    
    elec_pair = fr_retr_high.label{i};
    ylabel(strrep(sprintf('%s',elec_pair),'_',' '));
end
set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% set(ha(1:size(ha,1)),'FontSize',12);



suptitle(['Power spectral density in different task periods -  Grid Channels '])




%% PLV spectra for different task periods - Scalp Hippocampus


PLV_Pairs_Enc      = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_encoding.mat');
nChannelPairs      = PLV_Pairs_Enc.nChannelPairs_Sclp_HC
PLV_Pairs_Enc      = PLV_Pairs_Enc.PLV_Pairs_Scalp_enc_HC;
PLV_Pairs_Maint    = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_maintenance.mat');
PLV_Pairs_Maint    = PLV_Pairs_Maint.PLV_Pairs_Scalp_maint_HC;
PLV_Pairs_Fixation = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Scalp_Hippocampus coupling\Patient_42_Scalp-Hippocampus_PLV_Pairs_fixation.mat');
freqAxis_fix       = PLV_Pairs_Fixation.freqAxis_sclp
PLV_Pairs_Fixation = PLV_Pairs_Fixation.PLV_Pairs_Scalp_fix_HC;

PLV_Pairs_Grid_enc = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding.mat');
PLV_Pairs_Grid_maint = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance.mat');

%
fig = figure
set(gcf,'color','white')
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();
iSS_High = 4; % Set Size [6 8]
iSS_Low = 1;  % Set Size [4]
freqAxis = fr_maint_high.freq;
if length(dataBipolar_Ret_SS{iSS_Low}.label)<= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end
j = 3; % 3rd Bipolar channel corresponding to AHl2-3
for i = 1:nScalpChans
    Scalp_channel = dataBipolar.label{i};
    indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
    nSubplot_plv(i) = figLayout{indSubplot,2};
    axes(ha(nSubplot_plv(i)));
    plot(freqAxis_fix,abs(PLV_Pairs_Fixation{i+((j-1)*nScalpChans),iSS_High}),strPlotColors{1},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Enc{i+((j-1)*nScalpChans),iSS_High}),strPlotColors{2},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Maint{i+((j-1)*nScalpChans),iSS_Low}),strPlotColors{3},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Maint{i+((j-1)*nScalpChans),iSS_High}),strPlotColors{4},'LineWidth',2);
    ylim([0,1])
  
    
    elec_pair = fr_fix_high.label{i};
    ylabel(elec_pair);
end

set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')
% set(ha(1:size(ha,1)),'FontSize',18);
channel = strrep(dataBipolar_Ret_SS{iSS_Low}.label{nChannelPairs(i+((j-1)*nScalpChans))},'_',' ');
suptitle(['PLV spectra for different task periods on Scalp-',channel])
for i = 1:size(ha,1)
    
    if ~ismember(i,nSubplot_plv)
        set(ha(i),'XColor',[1,1,1],'YColor',[1,1,1],'Color',[1,1,1])
    end
    
end

%% PLV spectra for different task periods - Grid Hippocampus

PLV_Pairs_Grid_enc = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding.mat');
nChannelPairs_Grid = PLV_Pairs_Grid_enc.nChannelPairs;
PLV_Pairs_Grid_enc = PLV_Pairs_Grid_enc.PLV_Pairs;
PLV_Pairs_Grid_maint = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance.mat');
PLV_Pairs_Grid_maint = PLV_Pairs_Grid_maint.PLV_Pairs;

PLV_Pairs_Grid_fix = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_fixation.mat');
PLV_Pairs_Grid_fix = PLV_Pairs_Grid_fix.PLV_Pairs;

if length(dataBipolar_Ret_SS{iSS_Low}.label)>= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end

fig = figure;
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';

iSS_High = 4; % Set Size [6 8]
iSS_Low = 1;  % Set Size [4]
freqAxis = fr_maint_high.freq;

j = 4; % 4th Bipolar channel corresponding to AHl1-3
for i = 1:nGridChans
    axes(ha(plot_order(i)));
    plot(freqAxis_fix,abs(PLV_Pairs_Grid_fix{i+((j-1)*nGridChans),iSS_High}),strPlotColors{1},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High}),strPlotColors{2},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_Low}),strPlotColors{3},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High}),strPlotColors{4},'LineWidth',2);
    hold on;
    
    Signal_d1 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High});
    Signal_d2 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_Low});
    Signal_d3 = abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High});
    yShift1 = 0.05;
    color1 = 'm';
    yShift2 = 0.01;
    color2 = 'c';
    indMaxBand = Difference_Bar(Signal_d1,Signal_d3,freqAxis,[0 1],yShift1,color1);
    indMaxBand2 = Difference_Bar(Signal_d3,Signal_d1,freqAxis,[0 1],yShift2,color2);

    ylim([0,1])

    elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
    ylabel(elec_pair);
end
set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')


channel = strrep(dataBipolar_Ret_SS{iSS_Low}.label{nChannelPairs_Grid(i+((j-1)*nGridChans),1)},'_',' ');
suptitle(['PLV spectra for different task periods on Grid-',channel])


%% PLV spectra for different task periods - Grid Hippocampus-High Gamma -FOI- 0.5:0.5:120

PLV_Pairs_Grid_enc_gamma = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_gamma');
freqAxis = PLV_Pairs_Grid_enc_gamma.freqAxis;
nChannelPairs_Grid_gamma = PLV_Pairs_Grid_enc_gamma.nChannelPairs;
PLV_Pairs_Grid_enc_gamma = PLV_Pairs_Grid_enc_gamma.PLV_Pairs;
PLV_Pairs_Grid_maint_gamma = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_gamma.mat');
PLV_Pairs_Grid_maint_gamma = PLV_Pairs_Grid_maint_gamma.PLV_Pairs;
PLV_Pairs_Grid_fix_gamma = load('F:\Vasileios\Task Analysis\Data\Analysis Data\Grid_Hippocampus coupling\Patient_42_Grid-Hippocampus_PLV_Pairs_fixation_gamma.mat');
PLV_Pairs_Grid_fix_gamma = PLV_Pairs_Grid_fix_gamma.PLV_Pairs;



if length(dataBipolar_Ret_SS{iSS_Low}.label)>= nGridChans + AHL_chans + PHL_chans + nScalpChans
   error('Please ensure you adjusted the coupling parameters in the beginning of the script') 
end

fig = figure;
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';

iSS_High = 4; % Set Size [6 8]
iSS_Low = 1;  % Set Size [4]

j = 4; % 4th Bipolar channel corresponding to AHl1-3
for i = 1:nGridChans
    axes(ha(plot_order(i)));
%     plot(freqAxis,abs(PLV_Pairs_fix_gamma{i+((j-1)*nGridChans),iSS_High}),strPlotColors{1},'LineWidth',2);
    plot(freqAxis_fix,abs(PLV_Pairs_fix_gamma{i}),strPlotColors{1},'LineWidth',1); %temporary

    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_enc_gamma{i+((j-1)*nGridChans),iSS_High}),strPlotColors{2},'LineWidth',1);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint_gamma{i+((j-1)*nGridChans),iSS_Low}),strPlotColors{3},'LineWidth',1);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint_gamma{i+((j-1)*nGridChans),iSS_High}),strPlotColors{4},'LineWidth',1);
    hold on;
    
    Signal_d1 = abs(PLV_Pairs_Grid_maint_gamma{i+((j-1)*nGridChans),iSS_High});
    Signal_d2 = abs(PLV_Pairs_Grid_enc_gamma{i+((j-1)*nGridChans),iSS_High});
    yShift1 = 0.05;
    color1 = 'm';
    yShift2 = 0.01;
    color2 = 'c';
    indMaxBand = Difference_Bar(Signal_d1,Signal_d2,freqAxis,[0 1],yShift1,color1);
    indMaxBand2 = Difference_Bar(Signal_d2,Signal_d1,freqAxis,[0 1],yShift2,color2);

    ylim([0,1])

    elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
    ylabel(elec_pair);
end
set(ha(1:size(ha,1)),'XTickLabel',''); set(ha,'YTickLabel','')


channel = strrep(dataBipolar_Ret_SS{iSS_Low}.label{nChannelPairs_Grid(i+((j-1)*nGridChans),1)},'_',' ');
suptitle(['PLV spectra for different task periods in High Frequency Range on Grid-',channel])



%% PLV for Grid Hippocampus - FOI 0.5:10:120
PLV_Pairs_Grid_maint = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_gamma_smooth.mat'])
freqAxis = PLV_Pairs_Grid_maint.freqAxis;
PLV_Pairs_Grid_maint = PLV_Pairs_Grid_maint.PLV_Pairs;
PLV_Pairs_Grid_enc = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_gamma_smooth.mat'])
freqAxis = PLV_Pairs_Grid_enc.freqAxis;
PLV_Pairs_Grid_enc = PLV_Pairs_Grid_enc.PLV_Pairs;
strPlotColors = {'r','g','b','k',};
iSS_High = 4;
iSS_Low =1;


PLV_Pairs_Fixation = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_gamma_smooth.mat'])
PLV_Pairs_Fixation = PLV_Pairs_Fixation.PLV_Pairs;

figure;
ha =  tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
j = 4; % 4th Bipolar channel corresponding to AHl1-3
for i = 1:nGridChans
    axes(ha(plot_order(i)));
    plot(freqAxis,abs(PLV_Pairs_Fixation{i+((j-1)*nGridChans),iSS_High}),strPlotColors{1},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High}),strPlotColors{2},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High}),strPlotColors{4},'LineWidth',2);
    hold on;
      
    Signal_d1 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High});
    Signal_d2 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_Low});
    Signal_d3 = abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High});
    yShift1 = 0.05;
    color1 = 'm';
    yShift2 = 0.01;
    color2 = 'c';
    indMaxBand = Difference_Bar(Signal_d1,Signal_d3,freqAxis,[0 1],yShift1,color1);
    indMaxBand2 = Difference_Bar(Signal_d3,Signal_d1,freqAxis,[0 1],yShift2,color2);

    ylim([0,1])

    elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
    ylabel(elec_pair);
end

ttl = ['PLV Spectrum for Fixation-Encoding-Maintenance, Channel ', strrep(dataBipolar.label{Bip_chans(j-2)},'_','')]
suptitle(ttl)



%% PLV for Grid Hippocampus - FOI 1:1:25
PLV_Pairs_Grid_maint = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_maintenance_test.mat'])
freqAxis = 1:1:25;
PLV_Pairs_Grid_maint = PLV_Pairs_Grid_maint.PLV_Pairs;
PLV_Pairs_Grid_enc = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_encoding_test.mat'])
PLV_Pairs_Grid_enc = PLV_Pairs_Grid_enc.PLV_Pairs;
PLV_Pairs_Fixation = load([strVariableFolder,'Patient_',num2str(42,'%.2d '),'_',coupling_type,'_PLV_Pairs_fixation_test.mat'])
PLV_Pairs_Fixation = PLV_Pairs_Fixation.PLV_Pairs;

strPlotColors = {'r','g','b','k',};
iSS_High = 4;
iSS_Low =1;

figure;
ha =  tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
j = 4; % 4th Bipolar channel corresponding to AHl1-3
for i = 1:nGridChans
    axes(ha(plot_order(i)));
    plot(freqAxis,abs(PLV_Pairs_Fixation{i+((j-1)*nGridChans),iSS_High}),strPlotColors{1},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High}),strPlotColors{2},'LineWidth',2);
    hold on;
    plot(freqAxis,abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High}),strPlotColors{4},'LineWidth',2);
    hold on;
      
    Signal_d1 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High});
    Signal_d2 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_Low});
    Signal_d3 = abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High});
    Signal_d4 = abs(PLV_Pairs_Fixation{i+((j-1)*nGridChans),iSS_High})
    yShift1 = 0.05;
    color1 = 'm';
    yShift2 = 0.01;
    color2 = 'c';
    indMaxBand = Difference_Bar(Signal_d1,Signal_d3,freqAxis,[0 1],yShift1,color1);
    indMaxBand2 = Difference_Bar(Signal_d3,Signal_d1,freqAxis,[0 1],yShift2,color2);
    indMaxBand = Difference_Bar(Signal_d1,Signal_d4,freqAxis,[0 1],2*yShift1,'r');

    ylim([0,1])

    elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
    ylabel(elec_pair);
end

ttl = ['PLV Spectrum for Fixation-Encoding-Maintenance, Channel ', strrep(dataBipolar.label{Bip_chans(j-2)},'_','')]
suptitle(ttl)





%% PLV for Grid Hippocampus - FOI 30:5:100

PLV_Pairs_Grid_maint = load([strVariableFolder,'Patient_42_Grid-Hippocampus_PLV_Pairs_maintenance_f_30_5_100.mat'])
freqAxis = PLV_Pairs_Grid_maint.freqAxis; 
PLV_Pairs_Grid_maint = PLV_Pairs_Grid_maint.PLV_Pairs;

PLV_Pairs_Grid_enc = load([strVariableFolder, 'Patient_42_Grid-Hippocampus_PLV_Pairs_encoding_test_f_30_5_100.mat'])
PLV_Pairs_Grid_enc = PLV_Pairs_Grid_enc.PLV_Pairs;
PLV_Pairs_Fixation = load([strVariableFolder,'Patient_42_Grid-Hippocampus_PLV_Pairs_fixation_test_f_30_5_100.mat'])
PLV_Pairs_Fixation = PLV_Pairs_Fixation.PLV_Pairs;

strPlotColors = {'r','g','b','k',};
iSS_High = 4;
iSS_Low =1;


plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
for j = 3:8
    figure;
    ha =  tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])
    for i = 1:nGridChans
        axes(ha(plot_order(i)));
        semilogx(freqAxis,abs(PLV_Pairs_Fixation{i+((j-1)*nGridChans),iSS_High}),strPlotColors{1},'LineWidth',1);
        hold on;
        semilogx(freqAxis,abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High}),strPlotColors{2},'LineWidth',1);
        hold on;
        semilogx(freqAxis,abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High}),strPlotColors{4},'LineWidth',1);
        hold on;
        
        Signal_d1 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_High});
        Signal_d2 = abs(PLV_Pairs_Grid_maint{i+((j-1)*nGridChans),iSS_Low});
        Signal_d3 = abs(PLV_Pairs_Grid_enc{i+((j-1)*nGridChans),iSS_High});
        Signal_d4 = abs(PLV_Pairs_Fixation{i+((j-1)*nGridChans),iSS_High});
        yShift1 = 0.05;
        color1 = 'k';
        yShift2 = 0.05;
        color2 = 'r';
        indMaxBand = Difference_Bar(Signal_d1,Signal_d4,freqAxis,[0 1],yShift1,color1);
        indMaxBand2 = Difference_Bar(Signal_d4,Signal_d1,freqAxis,[0 1],yShift2,color2);
%         indMaxBand = Difference_Bar(Signal_d1,Signal_d4,freqAxis,[0 1],2*yShift1,'r');
        
        ylim([0,1])
        
        elec_pair = strrep(dataBipolar.label{i+AHL_chans},'_',' ');
        ylabel(elec_pair);
    end
    
    ttl = ['PLV Spectrum for Fixation-Encoding-Maintenance, Channel ', strrep(dataBipolar.label{Bip_chans(j-2)},'_','')]
    suptitle(ttl)
end