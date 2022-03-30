% load data file ('dataf') preprocessed with fieldtrip
% and show in eeglab viewer
%
% This function is provided as is. It only works for some specific type of
% data. This is a simple function to help the developer and by no mean
% an all purpose function.
%%%%%%%%%%%%%%%%%
% 171031 borec changed to include electrode labels by Ece.Boran@usz.ch
% eeglab cannot be in a function

% function [EEG] = fieldtrip2eeglab_171101(data)

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% if exist(dataf,'file')
%   load(dataf)
% end

% load chanlocs.mat
% EEG.chanlocs = chanlocs;
EEG.chanlocs = [];

for i=1:size(data.trial,2)
    EEG.data(:,:,i) = single(data.trial{i});
    %     EEG.event(i).latency = (i-1)*size(data.trial{i},2);
    %     EEG.event(i).type = 'ss 4';
    %     EEG.event(i).epoch = i;
    EEG.event((i-1)*4+1).latency = (i-1)*size(data.trial{i},2);
    EEG.event((i-1)*4+1).type = 'Fix';
    EEG.event((i-1)*4+1).epoch = i;
    
    EEG.event((i-1)*4+2).latency = (i-1)*size(data.trial{i},2)+1*data.fsample;
    EEG.event((i-1)*4+2).type = ['Stim',num2str(2*iSS+2)];
    EEG.event((i-1)*4+2).epoch = i;
    
    EEG.event((i-1)*4+3).latency = (i-1)*size(data.trial{i},2)+3*data.fsample;
    EEG.event((i-1)*4+3).type = 'Ret';
    EEG.event((i-1)*4+3).epoch = i;
    
    EEG.event((i-1)*4+4).latency = (i-1)*size(data.trial{i},2)+6*data.fsample;
    EEG.event((i-1)*4+4).type = 'Probe';
    EEG.event((i-1)*4+4).epoch = i;
end

EEG.setname    = 'sternberg'; %data.cfg.dataset;
EEG.filename   = '';
EEG.filepath   = '';
EEG.subject    = '';
EEG.group      = '';
EEG.condition  = '';
EEG.session    = [];
EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};
EEG.ref        = []; %'common';
%EEG.event      = [];
EEG.epoch      = [];
EEG.icawinv    = [];
EEG.icasphere  = [];
EEG.icaweights = [];
EEG.icaact     = [];
EEG.saved      = 'no';

%% Add channel names
for iChannel = 1:length(data.label)
    EEG.chanlocs(iChannel).labels = data.label{iChannel};
end

%%
eeglab redraw
% [EEG ALLEEG CURRENTSET] = eeg_retrieve(ALLEEG,1);
EEG = eeg_checkset( EEG );
% pop_eegplot( EEG, 1, 1, 1);
eeglab redraw;
