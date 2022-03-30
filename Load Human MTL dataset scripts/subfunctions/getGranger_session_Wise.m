function [gdata_session DeltaGranger_Maint] = getGranger_session_Wise(dataBipolar_Ret_SS_session,hipp_chanPairs,strChannelNameList)

%% Calculate Granger - Power Correlation for all channel pairs
clear gdata_session;
tStart = tic;
nSessions = size(dataBipolar_Ret_SS_session{1},2);
Fix_Data = dataBipolar_Ret_SS_session{1};
Enc_Data = dataBipolar_Ret_SS_session{2};
Maint_Data = dataBipolar_Ret_SS_session{3};

for iPair = 1:size(hipp_chanPairs,1)
    fprintf('Calculating Granger for pair %d for all set sizes\n',iPair)
    for iSS = 1:nSessions
        freq                   = [];
        freq.freqcfg           = [];
        freq.freqcfg.method    = 'mtmfft';
        freq.freqcfg.foi       = [1:1:100];%[4:1:100];%
        freq.freqcfg.output    = 'fourier';
        freq.freqcfg.tapsmofrq = 2;
        Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data{iSS});
        Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data{iSS});
        Fixation_freq          = ft_freqanalysis(freq.freqcfg, Fix_Data{iSS});
        
        grangercfg                          = [];
        grangercfg.method                   = 'granger';
        grangercfg.granger.conditional      = 'no';
        grangercfg.granger.sfmethod         = 'bivariate';
        grngChan1                           = strChannelNameList{hipp_chanPairs(iPair,1)};
        grngChan2                           = strChannelNameList{hipp_chanPairs(iPair,2)};
        grangercfg.channelcmb               = {grngChan1,grngChan2};
        gdata_session.Maint{iPair,iSS}      = ft_connectivityanalysis(grangercfg, Maintenance_freq);
%         gdata_session.Enc{iPair,iSS}        = ft_connectivityanalysis(grangercfg, Encoding_freq);
%         gdata_session.Fix{iPair,iSS}        = ft_connectivityanalysis(grangercfg, Fixation_freq);
        
        
        DeltaGranger_Maint{iSS} = gdata_session.Maint{iPair,iSS}.grangerspctrm(2,:)-gdata_session.Maint{iPair,iSS}.grangerspctrm(1,:);
        
    end
    tStop = toc(tStart);
    fprintf('\n\n\n\n\n\n\n\n\n\n Pair %d %f seconds elapsed \n',iPair,tStop)
end