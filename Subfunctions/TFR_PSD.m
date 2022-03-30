 %% Parameters

Frequency = 4:1:100;
% Log_freq = 10*log10(Frequency);
Logspace_freq = logspace(log10(Frequency(1)),log10(Frequency(end)),length(Frequency));

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.foi =4:1:100
cfg.tapsmofrq = 0.2*cfg.foi;
cfg.toi = -6:0.25:2 
cfg.t_ftimwin = 10./cfg.foi;
cfg.keeptrials = 'no';
TFR_psd = ft_freqanalysis(cfg,dataBipolar_SS{5})%dataBipolar_Scalp_SS{4})%

TFR_psd.label = strrep(TFR_psd.label,'_','');
TFR2 = TFR_psd
TFR2.powspctrm = 10*log10(TFR2.powspctrm);

%% Linear and logarithmic frequency TFR PSD
figure;
set(gcf,'color','white')
ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01])

for i =9:72%23:27%1:20%81:88
Chan_to_plot = i%Bip_chans(1)
clim = [0 1];
subplot_flag = 0;
axes(ha(i-8))

cfg = [];
cfg.baseline = [TFR2.time(1) TFR2.time(5)]
cfg.baselinetype = 'absolute';
TFR_baselined = ft_freqbaseline(cfg,TFR2)
PS = squeeze(TFR_baselined.powspctrm(Chan_to_plot,:,:));

% figure;
% set(gcf,'color','white')

if subplot_flag
    subplot(121)
end
TFR_psd_base = (squeeze(TFR_baselined.powspctrm(Chan_to_plot,:,:)))%/2%-0.5)/1.9%/3-1%
contourf(TFR_psd.time,TFR_psd.freq,TFR_psd_base,100,'LineColor','none')

%Axes properties
set(gca,'clim',clim,'yscale','log')
%     set(gca,'yscale','log')

set(gca,'ytick',ceil(logspace(log10(4),log10(100),5)),'YTickLabel',[4 10 20 40 100],'color',[0.01 0.01 0.56]) %background color of grid
% set(gca,'FontSize',20)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap 'jet' % to achieve correct colorbar for negative red peaks use flipud(jet) and -TFR_psd_base
colorbar
strTitle = ['TFR PSD at ', strrep(TFR_baselined.label{Chan_to_plot},'_',' ')];
% title(strTitle)
end
  
  %%
Chan_to_plot = 18
clim = [0 1];
subplot_flag = 0;
if subplot_flag
  subplot(122)
else
    figure
end
  contourf(TFR_psd.time,TFR_psd.freq,squeeze(TFR_baselined.powspctrm(Chan_to_plot,:,:))/1.5,100,'LineColor','none')
  %Axes properties
  set(gca,'clim',clim,'yscale','log')
  colormap jet
  colorbar
  set(gca,'color',[0.01 0.01 0.56])
  
  set(gca,'FontSize',20)
  set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
  colormap jet
  colorbar
%   strTitle = ['TFR PSD - ', strrep(TFR_baselined.label{Chan_to_plot},'_',' '),' Linear'];
%   title(strTitle)

  
  %%
  
  figure;imagesc(TFR_psd.time,[],squeeze(TFR2.powspctrm(90,:,:)),[0 30])
  set(gca,'ytick',round(linspace(4,100,6)),'yticklabel',round(Logspace_freq(round(linspace(4,97,6)))*10)/10)
  axis xy
  set(gca,'clim',[0 30])
  colormap jet
  colorbar
  
%% Calculate covariance and correlation of power

cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
% cfg.taper = 'hanning';
cfg.foi = [4:1:100];%ceil(logspace(log10(4),log10(100),97))%4:1:100%[0.5:0.5:30];
% numfoi = length(cfg.foi);
cfg.tapsmofrq = 2*cfg.foi;
cfg.toi = -6:0.25:2%-6:0.25:2
cfg.t_ftimwin = 1./cfg.foi;
cfg.keeptrials = 'no';
TFR_psd = ft_freqanalysis(cfg,dataBipolar_SS{5})
  
cfg = [];
% cfg.covariance         = 'yes';
cfg.keeptrials         = 'no';
cfg.removemean         = 'yes';
timelock = ft_timelockanalysis(cfg,TFR_psd);
Covariance_TFR =cov(timelock.avg)
d = sqrt(diag(Covariance_TFR));
r = Covariance_TFR./(d*d');
figure;
imagesc(TFR_psd.freq,TFR_psd.freq,r)
axis xy;
colormap jet;
colorbar;


%% Cross Correlation on channel AHL2-3 of 10 Hz and 20 Hz peaks in the powerspectrum


% Maintenance
nSet_Size = size(dataBipolar_SS,2);
dataBipolar_Ret_SS = {};
for iSS = 1:nSet_Size
    cfg = [];
    cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
    %         cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
    %         cfg.latency = [-1,-1/dataBipolar_SS{iSS}.fsample];
    
    
    dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

cfg = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';
cfg.output      = 'pow';
cfg.foi         = [11 23]

cfg.tapsmofrq   = 2;
cfg.channel = 83;
cfg.keeptrials  = 'yes';
fr_maint =ft_freqanalysis(cfg,dataBipolar_Ret_SS{4});
fr_maint.powspctrm = squeeze(fr_maint.powspctrm(:,1,:));
figure;
subplot(3,1,1)
plot([1:length(fr_maint.powspctrm)],10*log10(fr_maint.powspctrm(:,1)))
title('10 Hz')
subplot(3,1,2)
plot([1:length(fr_maint.powspctrm)],10*log10(fr_maint.powspctrm(:,2)))
title('20 Hz')
[r,lags] = xcorr(fr_maint.powspctrm(:,1),fr_maint.powspctrm(:,2),'normalized')
subplot(3,1,3)
stem(lags,r)
title('Cross Correlation')
suptitle('1 second of maintenance [-2 -1]')




phase_10Hz = angle(10*log10(fr_maint.powspctrm(:,1)));
phase_20Hz = angle(10*log10(fr_maint.powspctrm(:,2)));
figure;
subplot(3,1,1)
plot([1:length(fr_maint.powspctrm)],phase_10Hz)
title('Phase @ 10 Hz')
subplot(3,1,2)
plot([1:length(fr_maint.powspctrm)],phase_20Hz)
title('Phase @ 20 Hz')
subplot(3,1,3)
plot([1:length(fr_maint.powspctrm)],phase_20Hz-phase_10Hz)
title('Phase difference')
suptitle('2 seconds of maintenance [-2 0]')


%% Power of Power 
cfg = [];
cfg.method    = 'mtmconvol';
cfg.output    = 'pow';
cfg.taper     = 'hanning';
cfg.foi       = 4:1:100;
cfg.toi       = dataBipolar_Ret_SS{4}.time{1};
cfg.t_ftimwin = 1./cfg.foi;
cfg.keeptrials = 'yes';
cfg.channel = 83;
freq1 = ft_freqanalysis(cfg,dataBipolar_SS{4});



cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'pow';
cfg.taper     = 'hanning';
cfg.foilim    = [4 100];

cfg.keeptrials = 'no';

freq2 = ft_freqanalysis(cfg,freq1); %FieldTrip automatically converts the freq1 data to raw data.
                 
figure; imagesc(freq2.freq, freq2.freq, freq2.powspctrm/10)
axis xy
colorbar;
colormap jet;
