%In accordance with the script PLV_significance_calculation
% To be run together


Frequency = [4:1:100];
% Log_freq = 10*log10(Frequency);
Logspace_freq = logspace(log10(Frequency(1)),log10(Frequency(end)),length(Frequency));


%% Run PLV using cross spectrum
    %% Select pairs for analysis
    nPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real values
    ind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))','AHL2')));
    ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))','_GL_C2')));
    nPairs_ToRun = intersect(nPairs_ToRun,ind);
    nPairs_ToRun = intersect(nPairs_ToRun,ind2);
   
    
    %% Analysis for each pair of electrodes - Real values
    %% Save variables
    strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\'
    mkdir(strPaths.Results)
    strVariableFolder = [strPaths.Results,'Statistics\TFR_PLV\'];
    mkdir(strVariableFolder)
    
    for iPair = 1:length(nPairs_ToRun)
        nPair = nPairs_ToRun(iPair);
       
        %% Channel numbers
        nChannel_1 = nChannelPairs(nPair,1);
        nChannel_2 = nChannelPairs(nPair,2);
       
   

        
        %% Create data structure with the selected channels for each set size
        data_SinglePair_SS = cell(1,5);
        for iSS = 1:5
            cfg = [];
            cfg.channel = [nChannel_1,nChannel_2];
            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_SS{iSS});
        end
       
        %% Create dataset
        data_SinglePair_Rand_SS = data_SinglePair_SS;
       
        %% Phase coherence
        for iSS = 1:5
            %% Time-frequency analysis
            cfg = [];
            cfg.output     = 'powandcsd'; % fourier pow
            cfg.method     = 'mtmconvol';
            cfg.foi        = Logspace_freq;%4:1:30; % 5:0.5:30
            cfg.t_ftimwin  = 10./cfg.foi;
            cfg.tapsmofrq  = 0.3*cfg.foi; 
            fs             = dataBipolar_SS{iSS}.fsample;
            cfg.toi        = -6:0.25:2;%dataBipolar_Ret_SS{iSS}.time{1}%dataBipolar_Ret_SS{iSS}.time{1};%[-6:1/fs:2-1/fs] % -6:0.1:2
            cfg.keeptrials = 'yes';
            cfg.pad        = 'nextpow2';
            CrossFreq = ft_freqanalysis(cfg,data_SinglePair_Rand_SS{iSS});
           
            cfgFreq = cfg;
            
            cfg             = [];
            cfg.method      = 'plv';
            cfg.complex     = 'complex';
            PLV_Rand_SS{iSS}= ft_connectivityanalysis(cfg,CrossFreq);
            
         
            %% PLV using cross spectrum
%             PLV_SS{iSS} = squeeze(PLV_Rand_SS{iSS}.plvspctrm(1,2,:,:));
            PLV_SS{iSS} = squeeze(CrossFreq.crsspctrm);
            PLV_SS{iSS} = squeeze(nanmean(PLV_SS{iSS}./abs(PLV_SS{iSS}),1));
           
        end
       
        %% Store results for the channel pair
        strPairLabels = CrossFreq.label';
        freqAxis = CrossFreq.freq;
        timeAxis = CrossFreq.time;
        TimeInterval = [-6,2];
       
        %% Save real values
        save([strVariableFolder,'Patient_',num2str(42,'%.2d'),'_Real_PLV_Pair_',num2str(nPair,'%.4d'),'_',dataBipolar_SS{iSS}.label{nChannel_1},'_',dataBipolar_SS{iSS}.label{nChannel_2},'_mt_tsf_0_2f.mat'],'PLV_SS','strPairLabels','freqAxis','timeAxis','cfgFreq','TimeInterval','-v7.3')
%         clear PLV_SS strPairLabels freqAxis timeAxis TimeInterval CrossFreq
       
        %% Channel pair complete
        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))
        
    end


%% Log Frequency
clim = [0 1];
iSS_To_Plot = 5;
figure;
% subplot(2,2,1)
contourf(CrossFreq.time,CrossFreq.freq,abs(PLV_SS{iSS_To_Plot}),100,'LineColor','none')
%Axes properties
set(gca,'clim',clim,'yscale','log')
set(gca,'ytick',ceil(logspace(log10(6),log10(length(Frequency)),5)),'YTickLabel',[4,10,20,40,100]) %background color of grid
set(gca,'color',[0.01 0.01 0.56]);
set(gcf,'color','w');
set(gca,'FontSize',20)
set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])
colormap jet
colorbar
strTitle = ['TFR ',strrep(dataBipolar_SS{iSS_To_Plot}.label{nChannelPairs(nPair,1)},'_',' '),' - ' strrep(dataBipolar_SS{iSS_To_Plot}.label{nChannelPairs(nPair,2)},'_',' ')];
title(strTitle)
% 
% count = 3;
% ax = subplot(2,2,count);
% 
% tEnc = [find(CrossFreq.time==-4.5):find(CrossFreq.time==-3)];
% tMaint = [find(CrossFreq.time==-2):find(CrossFreq.time==-1)];
% 
% PLV_enc = mean(PLV_SS{iSS_To_Plot}(:,tEnc),2);
% % angle_enc = reshape(rad2deg(angle(PLV_enc)),size(PLV_enc,1)*length(tEnc),1);
% angle_enc = rad2deg(angle(PLV_enc));
% 
% figure;CircHist(angle_enc,36)
% % 
% % count = 4;
% % ax = subplot(2,2,count);
% PLV_maint = mean(PLV_SS{iSS_To_Plot}(:,tMaint),2);
% % angle_maint = reshape(rad2deg(angle(PLV_maint)),size(PLV_maint,1)*length(tMaint),1);
% angle_maint = rad2deg(angle(PLV_maint));
% 
% figure;CircHist(angle_maint,36)

% mean(angle_maint)
