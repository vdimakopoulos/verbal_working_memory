function [Granger_random_trial_sampling] = calculateGranger_random_trial_sampling(pID,iSS,nTrials_toSelect,dataMaint,dataEnc,channelPair)

k=1;
nIterations = 20;
ids = [42 38 37 40 45 44; 10 11 12 13 14 15]';
for iRand = 1:nIterations
    trialpool = 1:length(dataMaint.trial);
    trialsSelection = randperm(numel(trialpool),nTrials_toSelect);%randi([trialpool(1):trialpool(end)],nTrials_toSelect,1);
    cfg = [];
    cfg.trials = trialsSelection;
    Maint_Data = ft_selectdata(cfg,dataMaint);
    Enc_Data = ft_selectdata(cfg,dataEnc);

    strChannelNameList = dataEnc.label;
    freq                   = [];
    freq.freqcfg           = [];
    freq.freqcfg.method    = 'mtmfft';
    freq.freqcfg.output    = 'fourier';
    freq.freqcfg.taper     = 'hanning';
    freq.freqcfg.pad       = 20;
    % freq.freqcfg.tapsmofrq = 20;
    Maintenance_freq       = ft_freqanalysis(freq.freqcfg, Maint_Data);
    Encoding_freq          = ft_freqanalysis(freq.freqcfg, Enc_Data);
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    grangercfg.granger.conditional = 'no';
    grangercfg.granger.sfmethod = 'bivariate';
    grngChan1 = channelPair{1};
    grngChan2 = channelPair{2};
    grangercfg.channelcmb = {grngChan1,grngChan2};
    gdata.Maint    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata.Enc       = ft_connectivityanalysis(grangercfg, Encoding_freq);
    Granger_trial_cortex_hipp_enc(k,:) = gdata.Enc.grangerspctrm(2,:);
    Granger_trial_hipp_cortex_enc(k,:) =  gdata.Enc.grangerspctrm(1,:);
    Granger_trial_hipp_cortex_maint(k,:) = gdata.Maint.grangerspctrm(2,:);
    Granger_trial_cortex_hipp_maint(k,:) = gdata.Maint.grangerspctrm(1,:)
    k = k+1;

end

for j = 1:size(Granger_trial_hipp_cortex_maint,2)
    medianGranger_randTrials_maint_hc(j) = median(Granger_trial_hipp_cortex_maint(:,j));
end


for j = 1:size(Granger_trial_cortex_hipp_maint,2)
    medianGranger_randTrials_maint_ch(j) = median(Granger_trial_cortex_hipp_maint(:,j));
end

for j = 1:size(Granger_trial_cortex_hipp_enc,2)
    medianGranger_randTrials_enc_ch(j) = median(Granger_trial_cortex_hipp_enc(:,j));
end



for j = 1:size(Granger_trial_hipp_cortex_enc,2)
    medianGranger_randTrials_enc_hc(j) = median(Granger_trial_hipp_cortex_enc(:,j));
end

% Granger_random_trial_sampling_scalp{id} =Granger_random_trial_sampling_scalp{id-1};
if pID<10
    id = pID;
else
    idx = find(ids(:,1) == pID);
    id = ids(idx,2);
end
Granger_random_trial_sampling{id}.pID = pID
Granger_random_trial_sampling{id}.Median_gdata.Maint{1} = medianGranger_randTrials_maint_ch*100;
Granger_random_trial_sampling{id}.Median_gdata.Maint{2} = medianGranger_randTrials_maint_hc*100;
Granger_random_trial_sampling{id}.gdata.Maint{1} =Granger_trial_cortex_hipp_maint;
Granger_random_trial_sampling{id}.gdata.Maint{2} =Granger_trial_hipp_cortex_maint;


Granger_random_trial_sampling{id}.Median_gdata.Enc{1} = medianGranger_randTrials_enc_hc*100;
Granger_random_trial_sampling{id}.Median_gdata.Enc{2} = medianGranger_randTrials_enc_ch*100;
Granger_random_trial_sampling{id}.gdata.Enc{1} =Granger_trial_hipp_cortex_enc;
Granger_random_trial_sampling{id}.gdata.Enc{2} =Granger_trial_cortex_hipp_enc;
end