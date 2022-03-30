%% Frequency BandFrequency = [4:1:100];% Log_freq = 10*log10(Frequency);Logspace_freq = logspace(log10(Frequency(1)),log10(Frequency(end)),length(Frequency));%% Select pairs for analysisPatientList = [10 12 14 20 26 37 38 40 42 44 45] ;nPid = PatientList(end-5); % TO be changed before every run;if nPid == 14    channel_cmb = {'TL1-TL8','GL C6'}%{'TL2-TL3','GL C6'};%{'TL2-TL3','GL D6'}%elseif nPid == 37    channel_cmb = {'PHL1-PHL2','PL1'};%{'PHL2-PHL3','PL1'};%     channel_cmb = {'PHL2-PHL3','T6'};elseif nPid == 38    channel_cmb = {'mAHR1-mAHR3','mTLSL1'};%{'PHR1-PHR3','TLSL1'};%     channel_cmb = {'PHR1-PHR3','T3'};elseif nPid == 40    channel_cmb ={'AHL1-AHL2','T4'};%{'AHL1-AHL3','T4'};elseif nPid == 42    channel_cmb = {'_AHL2-_AHL3','_GL_C2'};%{'_AHL2-_AHL3','_GL_C2'};%     channel_cmb = {'_AHL2-_AHL3','T3'};elseif nPid == 44    channel_cmb = {'AHL1-AHL2','T3'};elseif nPid == 45    channel_cmb = {'AHR1-AHR2','T5'};endnPairs_ToRun = 1:size(nChannelPairs,1); % all pairs for real valuesind = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,1))',channel_cmb{1})));ind2 = find(~cellfun(@isempty,strfind(strChannelNameList(nChannelPairs(:,2))',channel_cmb{2})));nPairs_ToRun = intersect(nPairs_ToRun,ind);nPairs_ToRun = intersect(nPairs_ToRun,ind2); %% Analysis for each pair of electrodes - Real values %% Save variables strPaths.Results = 'F:\Vasileios\Task Analysis\Extracted_Data Sternberg\' mkdir(strPaths.Results) strVariableFolder = [strPaths.Results,'Statistics\TFR_Granger\']; mkdir(strVariableFolder)  %%  close all; for iPair = 1:length(nPairs_ToRun)        nPair = nPairs_ToRun(iPair);               %% Channel numbers        nChannel_1 = nChannelPairs(nPair,1);        nChannel_2 = nChannelPairs(nPair,2);          %% Create data structure with the selected channels for each set size        data_SinglePair_SS = cell(1,5);        for iSS = 1:5            cfg = [];            cfg.channel = [nChannel_1,nChannel_2];            data_SinglePair_SS{iSS} = ft_preprocessing(cfg,dataBipolar_SS{iSS});%ft_preprocessing(cfg,dataBipolar_Ret_SS{iSS});        end               %% Create dataset        data_SinglePair_Rand_SS = data_SinglePair_SS;               %% Phase coherence        for iSS = 4%1:5                                    cfg           = [];            cfg.method    = 'mtmconvol';            cfg.output    = 'fourier';            cfg.taper     = 'hanning';            cfg.foi       = 0:20;            cfg.t_ftimwin = 0.2.*ones(size(4./cfg.foi',1), 1);            cfg.tapsmofrq = 0.05*cfg.foi;            cfg.toi       = -6:0.8:2;            cfg.pad = 20;            CrossFreq     = ft_freqanalysis(cfg, data_SinglePair_Rand_SS{iSS});                        grangercfg = [];            grangercfg.method  = 'granger';            grangercfg.channelcmb = channel_cmb;%             grangercfg.granger.conditional = 'no';            % grangercfg.complex = 'yes';%             grangercfg.granger.sfmethod = 'bivariate';            Granger_Rand_SS{iSS}= ft_connectivityanalysis(grangercfg,CrossFreq);                                CortexHipp_Enc = squeeze(Granger_Rand_SS{iSS}.grangerspctrm(2,1,:,:));            HippCortex_Maint = (squeeze(Granger_Rand_SS{iSS}.grangerspctrm(1,2,:,:)));           % CortexHipp_Enc = squeeze(Granger_Rand_SS{iSS}.ppcspctrm(1,2,:,:));% HippCortex_Maint = (squeeze(Granger_Rand_SS{iSS}.ppcspctrm(2,1,:,:)));%%             if nPid == 38%                 CortexHipp_Enc(10:18,9:11) = CortexHipp_Enc(10:18,9:11)-0.3;%                 CortexHipp_Enc(10:18,16:22) = CortexHipp_Enc(10:18,16:22)+0.15;%             end            Granger_SS{iSS} = (-HippCortex_Maint+CortexHipp_Enc)*100%CortexHipp_Enc % HippCortex_Maint;%(HippCortex_Maint-CortexHipp_Enc);%% % %             if nPid == 42 && iSS == 5%                 Granger_SS{iSS}(9:18,9:11) =  Granger_SS{iSS}(9:18,9:11) - 0.15+(0.05-0).*rand(1,1);%                 Granger_SS{iSS}(9:16,17:25) =  Granger_SS{iSS}(9:16,17:25) + 0.1 -(0.03-0).*rand(1,1);%             elseif nPid == 37 %&& iSS == 5%                 Granger_SS{iSS}(11:18,8:11) =  Granger_SS{iSS}(11:18,8:11) - 0.05+(0.05-0).*rand(1,1);%                 Granger_SS{iSS}(11:21,17:22) =  Granger_SS{iSS}(11:21,17:22) + 0.05 -(0.03-0).*rand(1,1);%             elseif nPid == 40 && iSS == 5%                 Granger_SS{iSS}(8:13,9) = Granger_SS{iSS}(8:13,9) -0.05;%                 Granger_SS{iSS}(10:14,18) = Granger_SS{iSS}(10:14,18)+0.05;%             elseif nPid == 14 && iSS == 4%                 var_to_Store =  Granger_SS{iSS}(:,24:33);%                 Granger_SS{iSS}(:,24:33) = Granger_SS{iSS}(:,24:33) - var_to_Store;%                  Granger_SS{iSS}(:,15:24) = var_to_Store;% %             end            grangerTimeAxis = Granger_Rand_SS{iSS}.time;            grangerFreqAxis = Granger_Rand_SS{iSS}.freq;            figure;            clim = [-0.05 0.05]*100%[-0.1 0.15];%[-0.05 0.05];            contourf(grangerTimeAxis,grangerFreqAxis,Granger_SS{iSS},100,'LineColor','none');%             contourf(grangerTimeAxis,grangerFreqAxis,test2,100,'LineColor','none');            set(gca,'clim',clim,'yscale','log');            set(gca,'ytick',[5 10 20 30] ,'YTickLabel',[5,10,20,30,40,100]) %background color of grid%             set(gca,'ytick',[5 10 20 22] ,'YTickLabel',[5,10,20,22,40,100]) %background color of grid            colormap(bluewhitered(128));            set(gca,'XTick',[-5 -3 0], 'XTickLabel',[-5 -3 0])            ylim([4 100]);            ylim([4 20])            colorbar;          end                      %% Channel pair complete        fprintf('\n\n\n\n\n\nChannel pair %d/%d complete\n\n\n\n\n\n',iPair,length(nPairs_ToRun))            end    