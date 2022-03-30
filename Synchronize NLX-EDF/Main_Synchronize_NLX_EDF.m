
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


%% EDF
pID = 37;
patientInitials = 'PN';
strPaths.EDFdata = ['F:\Vasileios\Task Analysis\Data\Sternberg Task\EDF Data sternberg recordings\'];
strPaths.PatientEDFdata = sprintf('%s%d %s\\',strPaths.EDFdata,pID,patientInitials);
cd(strPaths.PatientEDFdata);
filesEDF = dir('*.edf');

%% NLX Paths
strPaths.NLXData = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\NLX Data for synchronization with EDF\';
strPaths.PatientNLXData = sprintf('%s%d %s\\',strPaths.NLXData,pID,patientInitials);
cd(strPaths.PatientNLXData);
filesNLX = dir('*.mat')
cd(strPaths.Main)
%% Synchronize EDF with NLX
for session = 4

cfg            = [];
cfg.dataset    = [strPaths.PatientEDFdata, filesEDF(session).name]

cfg.continuous = 'yes';
cfg.channel    = 'all';
dataEDF        = ft_preprocessing(cfg);
 
%% load and resample NLX stuff

NLX_dataset = filesNLX(session).name;
load([strPaths.PatientNLXData,NLX_dataset]);
loaded_data_NLX = data;

cfg=[];
cfg.resamplefs  = 256;
cfg.detrend     = 'no';
dataNLX         = ft_resampledata(cfg,loaded_data_NLX);
dataNLX.RejectedTrials = [];    
    
strSavePath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\EDF Data\';
mkdir(strSavePath);
%% synchronization
 
clear posr_sample dataEDFsynch
%  EDF channel for synch
ch_synch_edf = dataEDF.trial{1}(2,:) - dataEDF.trial{1}(1,:);
EEGch = 49:67;

for tr = 1:size(dataNLX.trial,2)
    
    ch_synch_nlx =   dataNLX.trial{1,tr}(15,:);
    fs = dataEDF.fsample;
    dataEDFsynch.time{tr}(1,:) = -6:1/fs:2-1/fs;% 2048
    dataEDFsynch.trial{tr}    = ones(length(EEGch)+length(dataNLX.label),length(ch_synch_nlx));
    
    if isempty(intersect(tr, dataNLX.RejectedTrials))
        
        ch_synch_nlx =   dataNLX.trial{1,tr}(15,:);
        [b,a] = butter(2,2*[5 30]/fs);
        L = min(length(ch_synch_edf), length(ch_synch_nlx));

        [r,lag] = xcorr(filtfilt(b,a,ch_synch_edf), filtfilt(b,a,ch_synch_nlx));
        [max_r,posr] = max(abs(r));
        posr = (lag(posr))/dataNLX.fsample;
        posr_sample(tr) = posr;
        
        % cut the NLX trial from dataEDFsynch
        good_int        = posr_sample(tr)*dataNLX.fsample+1:posr_sample(tr)*dataNLX.fsample+length(ch_synch_nlx);
        dataEDFsynch.trial{1,tr} = [dataEDF.trial{1,1}(EEGch,good_int); dataNLX.trial{1,tr}]; % trial "tr"
        dataEDFsynch.time{1,tr}   = dataNLX.time{1,tr};
        
%         figure, subplot(311),plot(lag/dataNLX.fsample,r)
%         subplot(3,1,[2 3]),
%         plot( dataEDF.time{1,1} - dataEDF.time{1,1}(1)-posr, filtfilt(b,a,ch_synch_edf),'b')
%         hold on,
%         plot( dataNLX.time{1,1} - dataNLX.time{1,1}(1)  , filtfilt(b,a,ch_synch_nlx) ,'r')
%         hold on,
%         plot( dataEDF.time{1,1} - dataEDF.time{1,1}(1)-posr,filtfilt(b,a,ch_synch_edf)+1000,'b')
%         grid on
     end
end

dataEDFsynch.label          = [dataEDF.label(EEGch);   dataNLX.label];
dataEDFsynch.fsample        = dataNLX.fsample;
dataEDFsynch.RejectedTrials = dataNLX.RejectedTrials;

dataEDF = dataEDFsynch;

dataEDF.setsize = TrialInformationTable.SetSize;
dataEDF.cfg     = dataNLX.cfg;

dataEDF.label = strrep(dataEDF.label,'EEG','')
dataEDF.label = strrep(dataEDF.label,'-REF','')

%%
fileName = sprintf('Sternberg_EDF_Patient %d %s session %d',pID, patientInitials,6);
save([strSavePath, fileName],'dataEDF','dataNLX') 
 
end