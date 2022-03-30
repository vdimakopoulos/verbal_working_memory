Frequency = 4:1:100;
% Log_freq = 10*log10(Frequency);
Logspace_freq = logspace(log10(Frequency(1)),log10(Frequency(end)),length(Frequency));

cfg = [];
cfg.method = 'mtmconvol';
% cfg.taper  = 'dpss';
cfg.output = 'pow';%'pow';
cfg.foi = [4:1:100]%[0.5:0.5:30];
cfg.tapsmofrq = 0.2*cfg.foi;
cfg.toi = -6:0.25:2% dataBipolar_SS{5}.time{1};
cfg.t_ftimwin = 10./cfg.foi;
cfg.keeptrials = 'yes';
TFR_PSD = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{4})%   dataBipolar_SS{4})

freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.foi       = [4:1:100];
freq.freqcfg.output    = 'fourier';
freq.freqcfg.tapsmofrq = 2;
freqdata1           = ft_freqanalysis(freq.freqcfg, dataBipolar_SS{4});

grangercfg = [];
grangercfg.method  = 'granger';
grangercfg.channelcmb = {'_AHL2-_AHL3','_GL_C2'};
grangercfg.granger.conditional = 'no';
grangercfg.granger.sfmethod = 'bivariate';

gdata = [];
gdata.g1_bivar_reg      = ft_connectivityanalysis(grangercfg, freqdata1);

figure;
plot(gdata.g1_bivar_reg.freq,gdata.g1_bivar_reg.grangerspctrm,'LineWidth',2);
set(gca,'FontSize',20);
Pair1 = strrep(strcat(grangercfg.channelcmb{1},'-',grangercfg.channelcmb{2}),'_',' ');
Pair2 = strrep(strcat(grangercfg.channelcmb{2},'-',grangercfg.channelcmb{1}),'_',' ');
legend(Pair1,Pair2);
xlim([4 100])



%% Calculate granger causality per task period for one pair
nPid = 38; %patient ID
nSet_Size = size(dataBipolar_SS,2);
Maint_Data = {};
iSS_current = 4; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]

for iSS = iSS_current% 1:nSet_Size
    cfg = [];
    cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
    Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

nSet_Size = size(dataBipolar_SS,2);
Enc_Data = {};
for iSS = iSS_current%1:nSet_Size
    cfg = [];
    cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
    Enc_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end

nSet_Size = size(dataBipolar_SS,2);
Fix_Data = {};
for iSS = iSS_current%1:nSet_Size
    cfg = [];
    cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample];
    Fix_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
end
iSS = 4;
freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.output    = 'fourier';
freq.freqcfg.taper = 'hanning'
freq.freqcfg.pad = 20;
Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS});
Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS});
Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS});



grangercfg = [];
grangercfg.method  = 'granger';
if nPid == 42
    grangercfg.channelcmb = {'_PHL2-_PHL3','_GL_C2'};%{'_AHL2','_T4'};%%{'_AHL2-_AHL3','_GL_C2'};2
    %     grangercfg.channelcmb = {'_AHL2','_GL_C2'}%{'_GL_C2','_AHL2-_AHL3'} %Patient 42 DS
    %     grangercfg.channelcmb = {'_AHL2-_AHL3','_T'}; %{'_GL_C2','_AHL2-_AHL3'} %Patient 42 DS
elseif nPid == 37
%         grangercfg.channelcmb = {'AHR2','PL1'};
    
        grangercfg.channelcmb = {'AHL2','PL1'};
%     grangercfg.channelcmb = {'PHL2-PHL3','T4'};
elseif nPid == 44
    grangercfg.channelcmb = {'AHL2-AHL3','T3'}%{'AHL1-AHL2','T3'}
elseif nPid == 38
        grangercfg.channelcmb = {'mPHR2','mTLSL1'}%{'mPHR1-mPHR2','mTLSL1'};%{'mPHR1','mTLSL1'}; %Patient 38 NP
    %     grangercfg.channelcmb = {'PHR1-PHR3','TLSL1'}; %Patient 38 NP
%     grangercfg.channelcmb = {'PHL2-PHL3','T3'}; %Patient 38 NP {'PHR2-PHR3','T5'}
elseif nPid == 40
    grangercfg.channelcmb ={'AHR1-AHR2','T4'}
elseif nPid == 14
    grangercfg.channelcmb = {'TL1-TL8','GL22'};
elseif nPid == 45
    grangercfg.channelcmb = {'PHR1-PHR2','T3'}; %{'AHR2-AHR3','T5'}
end
grangercfg.granger.conditional = 'no';
grangercfg.granger.sfmethod = 'bivariate';

gdata = [];
gdata.Maint    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
gdata.Enc      = ft_connectivityanalysis(grangercfg, Encoding_freq);
gdata.Fix      = ft_connectivityanalysis(grangercfg, Fixation_freq);

% close all;
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 0.2];
% figure;ft_connectivityplot(cfg, gdata.Maint);
% figure;ft_connectivityplot(cfg, gdata.Enc);


light_blue = [0.30,0.75,0.93];
light_red       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',light_red};
freq_ax    = gdata.Fix.freq;
figure;
% subplot(1,2,1)
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(1,:)*2,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(gdata.Enc.freq,gdata.Enc.grangerspctrm(2,:)/6,'Color',Colors{2},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(1,:)*2,'Color',Colors{3},'LineWidth',3);
semilogx(gdata.Maint.freq,gdata.Maint.grangerspctrm(2,:)/6,'Color',Colors{4},'LineWidth',3);
xlim([4 20])
% ylim([0 0.2]*100)
set(gca,'FontSize',16);
xlabel('Frequency (Hz)');
ylabel('Granger')
set(gca,'box', 'off')



%% Calculate Granger for different number of trials - convergence
iSS_current = 5; % Set Size [1 2 3 4 5] correspond to [4],[6],[8],[6 8],[4 6 8]
k=1;
nIterations = 50;

% for iTrial = 2:2:length(dataBipolar_SS{iSS_current+1}.trial)
for iRand = 1:nIterations
    trialpool = 1:length(dataBipolar_SS{iSS_current}.trial);
    nTrials_toSelect = size(TrialInformationTable_Scalp,1)-size(TrialInformationTable,1);
    trialsSelection = randi([trialpool(1) trialpool(end)],1,nTrials_toSelect);
    nPid = 37; %patient ID
    nSet_Size = size(dataBipolar_SS,2);
    Maint_Data = {};
    for iSS = iSS_current% 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Maint_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    nSet_Size = size(dataBipolar_SS,2);
    Enc_Data = {};
    for iSS = iSS_current%1:nSet_Size
        cfg = [];
        cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Enc_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    
    nSet_Size = size(dataBipolar_SS,2);
    Fix_Data = {};
    for iSS = iSS_current%1:nSet_Size
        cfg = [];
        cfg.latency = [-6,-5-1/dataBipolar_SS{iSS}.fsample];
        cfg.trials = trialsSelection;
        Fix_Data{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
    end
    freq                   = [];
    freq.freqcfg           = [];
    freq.freqcfg.method    = 'mtmfft';
    freq.freqcfg.foi       = 4:1:100;
    freq.freqcfg.output    = 'fourier';
    freq.freqcfg.tapsmofrq = 2;
    Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS_current});
    Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS_current});
    Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS_current});
    
    
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    if nPid == 42
        %         grangercfg.channelcmb = {'_AHL2-_AHL3','_GL_C2'}; %{'_GL_C2','_AHL2-_AHL3'} %Patient 42 DS
        grangercfg.channelcmb = {'_AHL2-_AHL3','T4'}; %{'_GL_C2','_AHL2-_AHL3'} %Patient 42 DS
    elseif nPid == 37
        %         grangercfg.channelcmb = {'PHL1-PHL2','PL1'};
        grangercfg.channelcmb = {'PHL1-PHL2','T6'};
    elseif nPid == 44
        grangercfg.channelcmb = {'AHL1-AHL2','T3'}
    elseif nPid == 38
        grangercfg.channelcmb = {'PHR1-PHR3','TLSL1'}; %Patient 38 NP
        %         grangercfg.channelcmb = {'PHR1-PHR3','T5'}; %Patient 38 NP
    elseif nPid == 40
        grangercfg.channelcmb ={'AHL1-AHL2','T4'}
    elseif nPid == 14
        grangercfg.channelcmb = {'TL1-TL8','GL22'};
    elseif nPid == 45
        grangercfg.channelcmb = {'AHR2-AHR4','T5'};
    end
    grangercfg.granger.conditional = 'no';
    grangercfg.granger.sfmethod = 'bivariate';
    
    gdata = [];
    gdata.Maint    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata.Enc      = ft_connectivityanalysis(grangercfg, Encoding_freq);
    gdata.Fix      = ft_connectivityanalysis(grangercfg, Fixation_freq);
    
    % close all;
    cfg           = [];
    cfg.parameter = 'grangerspctrm';
    cfg.zlim      = [0 0.2];
    % figure;ft_connectivityplot(cfg, gdata.Maint);
    % figure;ft_connectivityplot(cfg, gdata.Enc);
    Granger_trial_cortex_hipp_enc(k,:) = gdata.Enc.grangerspctrm(1,:);
    Granger_trial_hipp_cortex_enc(k,:) =  gdata.Enc.grangerspctrm(2,:);
    Granger_trial_hipp_cortex_maint(k,:) = gdata.Maint.grangerspctrm(1,:);
    Granger_trial_cortex_hipp_maint(k,:) = gdata.Maint.grangerspctrm(2,:)
    k = k+1;
    
end

% Visualize for maintenance
figure;
for iRand = 1:nIterations
    light_blue = [0.30,0.75,0.93];
    light_red       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',light_red};
    freq_ax    = gdata.Fix.freq;
    semilogx(gdata.Maint.freq,Granger_trial_hipp_cortex_maint(iRand,:)*100,'Color','r','LineWidth',3);
    hold on;
    xlim([4 30])
    set(gca,'FontSize',16);
    xlabel('Frequency (Hz)');
    ylabel('Granger')
    set(gca,'box', 'off')
end

for j = 1:size(Granger_trial_hipp_cortex_maint,2)
    medianGranger_randTrials_maint_hc(j) = median(Granger_trial_hipp_cortex_maint(:,j));
end
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_hc*100,'k','LineWidth',3);
xlim([4 30])

for j = 1:size(Granger_trial_cortex_hipp_maint,2)
    medianGranger_randTrials_maint_ch(j) = median(Granger_trial_cortex_hipp_maint(:,j));
end

% Visualize for encoding
figure;
for iRand = 1:nIterations
    light_blue = [0.30,0.75,0.93];
    light_red       = [1 0.45 0.45];
    Colors     = {'b',light_blue,'r',light_red};
    freq_ax    = gdata.Fix.freq;
    semilogx(gdata.Enc.freq,Granger_trial_cortex_hipp_enc(iRand,:)*100,'Color','b','LineWidth',3);
    hold on;
    xlim([4 30])
    set(gca,'FontSize',16);
    xlabel('Frequency (Hz)');
    ylabel('Granger')
    set(gca,'box', 'off')
    k = k+1;
end

for j = 1:size(Granger_trial_cortex_hipp_enc,2)
    medianGranger_randTrials_enc_ch(j) = median(Granger_trial_cortex_hipp_enc(:,j));
end
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_ch*100,'k','LineWidth',3);
xlim([4 30])


for j = 1:size(Granger_trial_hipp_cortex_enc,2)
    medianGranger_randTrials_enc_hc(j) = median(Granger_trial_hipp_cortex_enc(:,j));
end

figure;
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_ch*100,'Color',Colors{1},'LineWidth',3);
hold on;
semilogx(gdata.Maint.freq,medianGranger_randTrials_enc_hc*100,'Color',Colors{2},'LineWidth',3);
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_hc*100,'Color',Colors{3},'LineWidth',3);
semilogx(gdata.Maint.freq,medianGranger_randTrials_maint_ch*100,'Color',Colors{4},'LineWidth',3);

xlim([4 30])

%% %% Calculate granger causality per task period for multiple pairs

nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','AHR2')));
ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','PL')));
nPairs_ToRun = intersect(nPairs_ToRun,ind);
nPairs_ToRun = intersect(nPairs_ToRun,ind2);
gdata_Pairs = [];
% strChannelNameList = strrep(strChannelNameList,'_','');
% Maintenance_freq.label = strChannelNameList;
% Encoding_freq.label = strChannelNameList;
% Fixation_freq.label = strChannelNameList;
for iPair = 1:length(nPairs_ToRun)
    nPair = nPairs_ToRun(iPair);
    
    % Channel numbers
    nChannel_1 = nChannelPairs(nPair,1);
    nChannel_2 = nChannelPairs(nPair,2);
    
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    grangercfg.channelcmb = {char(strChannelNameList(nChannel_1)), char(strChannelNameList(nChannel_2))};%{'_AHL2-_AHL3','_GL_C2'}; %{'_GL_C2','_AHL2-_AHL3'}
    grangercfg.granger.conditional = 'no';
    % grangercfg.complex = 'yes';
    grangercfg.granger.sfmethod = 'bivariate';
    
    gdata_Pairs.Maint{iPair}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata_Pairs.Enc{iPair}      = ft_connectivityanalysis(grangercfg, Encoding_freq);
    gdata_Pairs.Fix{iPair}      = ft_connectivityanalysis(grangercfg, Fixation_freq);
    
    DeltaGranger_maint(iPair,:) =  gdata_Pairs.Maint{iPair}.grangerspctrm(1,:) -  gdata_Pairs.Maint{iPair}.grangerspctrm(2,:);
    DeltaGranger_enc(iPair,:) =  gdata_Pairs.Enc{iPair}.grangerspctrm(1,:) -  gdata_Pairs.Enc{iPair}.grangerspctrm(2,:);
    
    freqAxis = gdata_Pairs.Maint{1}.freq;
    FreqBand = [9 18];
    [~,indFreq1] = min(abs(freq_ax-FreqBand(1)));
    [~,indFreq2] = min(abs(freq_ax-FreqBand(2)));
    Dgranger_maint_InBand(iPair) = max(DeltaGranger_maint(iPair,indFreq1:indFreq2));
    Dgranger_enc_InBand(iPair) = mean(DeltaGranger_enc(iPair,indFreq1:indFreq2));
end



%% Visualization for every pair

light_blue = [0.30,0.75,0.93];
light_red       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',light_red};
freq_ax    = gdata_Pairs.Fix{iPair}.freq;
% GridPanels_Motor_activity = [12,23,24,32]; %B4,C7,C8,D8
% GridPanels_Sensory_activity = [22,31];     %C6,D7
% Aud_Cortex = [17 18 25 26];                %C1,C2,D1,D2
% Vis_Cortex = [57:59];                      %H1, H2, H3

figure('units','normalized','outerposition',[0 0 1 1]);
% ha = tight_subplot(8,8,[.01 .03],[.1 .01],[.01 .01]);
ha = tight_subplot(5,9,[.01 .03],[.1 .01],[.01 .01])
figLayout = Get_Figure_Scalp_Layout_Information();

% plot_order = [57:-8:1;58:-8:2;59:-8:3;60:-8:4;61:-8:5;62:-8:6;63:-8:7;64:-8:8]';
figure
for i = 1:length(nPairs_ToRun)
    Scalp_channel = strrep(strChannelNameList(nChannelPairs(nPairs_ToRun(i),2)),'_','')
    indSubplot = find(strcmpi(figLayout(:,1),Scalp_channel));
    nSubplot(i) = figLayout{indSubplot,2};
    
    axes(ha(nSubplot(i)));
    semilogx(freq_ax,gdata_Pairs.Enc{i}.grangerspctrm(1,:),'Color',Colors{1},'LineWidth',2);
    hold on;
    semilogx(freq_ax,gdata_Pairs.Enc{i}.grangerspctrm(2,:),'Color',Colors{2},'LineWidth',2);
    semilogx(freq_ax,gdata_Pairs.Maint{i}.grangerspctrm(1,:),'Color',Colors{3},'LineWidth',2);
    semilogx(freq_ax,gdata_Pairs.Maint{i}.grangerspctrm(2,:),'Color',Colors{4},'LineWidth',2);
    %     channel_lbl = strrep(strChannelNameList(i+8),'_',' ');
    ylabel(Scalp_channel);
    xlim([4 30]);
    ylim([0 0.2]);
    
end

set(ha(1:size(ha,1)),'Fontsize',10)
% set(ha(plot_order(GridPanels_Motor_activity)),'XColor','b','YColor','b');
% set(ha(plot_order(GridPanels_Sensory_activity)),'XColor','r','YColor','r');

% set(ha(plot_order(Aud_Cortex)),'XColor',[0.06 1 1],'YColor',[0.06 1 1]);
% set(ha(plot_order(Vis_Cortex)),'XColor',[0 0 0.53],'YColor',[0 0 0.53]);


lg = legend('Grid to HC encoding','HC to Grid encoding','Grid to HC maintenance','HC to Grid maintenance')
Lgnd_pos = [0.604134137691255,0.901149476394398,0.116145830663542,0.09198355343947];
set(lg,'Location','none','Position',Lgnd_pos,'Fontsize',12);
suptitle('Granger causality spectrum AHL2-3 - Grid')
%%
% Granger heatmat
Granger_In_Band_SS = [];
FreqBand = [15,20];
[~,indFreq1] = min(abs(freq_ax-FreqBand(1)));
[~,indFreq2] = min(abs(freq_ax-FreqBand(2)));

for iPair = 1:length(nPairs_ToRun)
    Granger_temp = gdata_Pairs.Maint{iPair};
    Granger_temp = Granger_temp.grangerspctrm(1,indFreq1:indFreq2)-Granger_temp.grangerspctrm(2,indFreq1:indFreq2);
    Granger_In_Band_SS(iPair) = mean(Granger_temp);
end
GrangerInBandToPlot = reshape(Granger_In_Band_SS(:),1,4);%reshape(Granger_In_Band_SS(:),8,8)
figure;
imagesc(GrangerInBandToPlot*100,[0 0.1]*100)
ax1 = gca;
ax1.XTick = 1:8;
ax1.YTick = 1:8;
ax1.XTickLabel = {'A','B','C','D','E','F','G','H'};
ax1.YTickLabel = 1:8;
ax1.YDir = 'normal'

colorbar
colormap(bluewhitered_pos)
set(gca,'FontSize',16,'box','off','clim',[0,10])
%% Visualization
figure;
cfg = [];
cfg.parameter = 'grangerspctrm';
% cfg.xlim = [Logspace_freq(1) Logspace_freq(end)];
cfg.zlim = [0 0.2];
cfg.xlim = [4 40];
% figure;
ft_connectivityplot(cfg,gdata.Fix)
suptitle('Fixation');
figure;
ft_connectivityplot(cfg,gdata.Enc)
suptitle('Encoding');
figure;
ft_connectivityplot(cfg,gdata.Maint)
suptitle('Maintenance');


%%

figure;
plot(gdata.Maint.freq,gdata.Maint.grangerspctrm,'LineWidth',2);
set(gca,'FontSize',20);
title('Granger spectrum during maintenance');
Pair1 = strrep(strcat(grangercfg.channelcmb{1},'-',grangercfg.channelcmb{2}),'_',' ');
Pair2 = strrep(strcat(grangercfg.channelcmb{2},'-',grangercfg.channelcmb{1}),'_',' ');
h = legend(Pair1,Pair2);
set(h,'Location','northwest')
xlim([4 100]);
ylim([0 0.2]);


plot(gdata.Enc.freq,gdata.Enc.grangerspctrm,'LineWidth',2);
set(gca,'FontSize',20);
title('Granger spectrum during encoding');
Pair1 = strrep(strcat(grangercfg.channelcmb{1},'-',grangercfg.channelcmb{2}),'_',' ');
Pair2 = strrep(strcat(grangercfg.channelcmb{2},'-',grangercfg.channelcmb{1}),'_',' ');
h = legend(Pair1,Pair2);
set(h,'Location','northwest')
xlim([4 100]);
ylim([0 0.2]);

figure;
plot(gdata.Fix.freq,gdata.Fix.grangerspctrm,'LineWidth',2);
set(gca,'FontSize',20);
title('Granger spectrum during fixation');
Pair1 = strrep(strcat(grangercfg.channelcmb{1},'-',grangercfg.channelcmb{2}),'_',' ');
Pair2 = strrep(strcat(grangercfg.channelcmb{2},'-',grangercfg.channelcmb{1}),'_',' ');
h = legend(Pair1,Pair2);
set(h,'Location','northwest')
xlim([4 100]);
ylim([0 0.2]);

