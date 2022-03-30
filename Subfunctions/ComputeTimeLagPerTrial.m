function [xcorrelogram,timeLags] = ComputeTimeLagPerTrial(dataBipolar,nChannelPairs,CxHipp_pair,foi,saveflag, pID,TaskPeriod)

%%
% Patient 42 -> maint foi [11 20]
% Patient 42 -> encod foi [13 18]

if nargin<5
    saveflag = 0;
elseif nargin ==5
    if saveflag
        error('Please provide the patient ID and the Task Period strings in order to save the figure');
    end
end
%% Find the hippocampal and the cortical channel
Hipp_chan = nChannelPairs(CxHipp_pair,1);
Cx_chan = nChannelPairs(CxHipp_pair,2);
nTrials = length(dataBipolar.trial);

%% Bandpass filtering
%
cfg = [];
cfg.demean = 'yes';
dataBipolar = ft_preprocessing(cfg,dataBipolar);

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq =  foi;
dataBipolar_Filtered = ft_preprocessing(cfg,dataBipolar);
% dataBipolar_Filtered = dataBipolar
%% Calculate cross correlation between hipp time series and  cortex time series
for i = 1:nTrials
    ts_Hipp = dataBipolar_Filtered.trial{i}(Hipp_chan,:);
    ts_Cortex = dataBipolar_Filtered.trial{i}(Cx_chan,:);
    [cor{i},lag{i}]=xcorr(ts_Hipp,ts_Cortex,'normalized');
end
%% Calculate the average envelope of the xcorrelogram across trials
cor_mat = cell2mat(cor);
cor_mat = reshape(cor_mat,nTrials,length(cor{1}));
std_cor = std(cor_mat);
avgEnvelope = abs(hilbert(cor{1}));
for i = 2:nTrials
    envelopeXcorr{i} = abs(hilbert(cor{i}));
    avgEnvelope = avgEnvelope+envelopeXcorr{i};
end
avgEnvelope = avgEnvelope/nTrials;

%% Visualization
timeLength  = length(dataBipolar.time{1});
fs = dataBipolar.fsample;
Lags = [-(timeLength-1)*1/fs:1/fs:(timeLength-1)*1/fs];
figure;
% for i = 1:nTrials
%     plot(Lags,cor{i},'k')
%     hold on;
% end
plot(Lags,avgEnvelope);
xlabel('Lags (s)');
title('Amplitude envelope cross correlation');
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlim([-0.4 0.4]);

set(gca,'FontSize',13,'TickDir','Out','box','off');
set(gcf,'Color','white');

hold on;
inbetween = [avgEnvelope(1:end/2+1)-std(avgEnvelope(1:end/2+1)), fliplr(avgEnvelope(1:end/2+1)+std(avgEnvelope(1:end/2+1)))];
% inbetween = [avgEnvelope(1:end/2+1)-std_cor((1:end/2+1)), fliplr(avgEnvelope(1:end/2+1)+std_cor((1:end/2+1)))];
x2 = [Lags(1:end/2+1), fliplr(Lags(1:end/2+1))];
hipp_std = fill(x2,inbetween,'r');
set(hipp_std,'FaceAlpha',0.3,'EdgeColor','none');

hold on;
inbetween = [avgEnvelope(end/2:end)-std(avgEnvelope(end/2:end)), fliplr(avgEnvelope(end/2:end)+std(avgEnvelope(end/2:end)))];
% inbetween = [avgEnvelope(end/2:end)-std_cor((end/2:end)), fliplr(avgEnvelope(end/2:end)+std_cor((end/2:end)))];
x2 = [Lags(end/2:end), fliplr(Lags(end/2:end))];
hipp_std = fill(x2,inbetween,'b');
set(hipp_std,'FaceAlpha',0.3,'EdgeColor','none');

%% Textboxes
annotation(gcf,'textbox',...
    [0.135714285714286 0.151380953942027 0.260714278476579 0.0761904746294022],...
    'Color',[1 0 0],...
    'String',{'Hipp to Cortex'},...
    'FontWeight','bold',...
    'FontSize',13,...
    'EdgeColor','none');

% Create textbox
annotation(gcf,'textbox',...
    [0.674999999999999 0.151380953942027 0.260714278476579 0.0761904746294022],...
    'Color',[0 0 1],...
    'String',{'Cortex to Hipp'},...
    'FontWeight','bold',...
    'FontSize',13,...
    'EdgeColor','none');
%% Return the variables
xcorrelogram = avgEnvelope;
timeLags = Lags;

if saveflag
    strPath = ['F:\Vasileios\Task Analysis\Analysis Results\','TimeLag Results\'];
    mkdir(strPath);
    strFilename = ['S',int2str(pID),'_Amplitude_Envelope_CrossCorrelation_',TaskPeriod];
    saveas(gcf,[strPath,strFilename]);
    saveas(gcf,[strPath,strFilename],'png');
end