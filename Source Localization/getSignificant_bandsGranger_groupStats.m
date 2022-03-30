function [sgnf_bands_enc sgnf_bands_maint stat_maint_freq stat_enc_freq] = getSignificant_bandsGranger_groupStats()
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerEnc_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\sourceGrangerMaint_group.mat')
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\freq.mat')
clear dataMaint dataEnc
dataMaint{1}.label = 'Hipp - Temporal Superior Lobe';
dataMaint{2}.label = 'Temporal Superior Lobe - Hipp';
dataMaint{1}.dimord = 'subj_freq';
dataMaint{2}.dimord = 'subj_freq';
dataMaint{1}.freq = freq;
dataMaint{2}.freq = freq;

dataEnc{1}.label = 'Hipp - Temporal Superior Lobe';
dataEnc{2}.label = 'Temporal Superior Lobe - Hipp';
dataEnc{1}.dimord = 'subj_freq';
dataEnc{2}.dimord = 'subj_freq';
dataEnc{1}.freq = freq;
dataEnc{2}.freq = freq;


temp_maint_1 = {sourceGrangerMaint_group{1,:}};
temp_maint_2 = {sourceGrangerMaint_group{2,:}};
temp_enc_1 = {sourceGrangerEnc_group{2,:}};
temp_enc_2 = {sourceGrangerEnc_group{1,:}};

dataAll_subj_maint_dir_1 = reshape(cell2mat(temp_maint_1),length(sourceGrangerEnc_group),length(freq));
dataAll_subj_maint_dir_2 = reshape(cell2mat(temp_maint_2),length(sourceGrangerEnc_group),length(freq));


dataAll_subj_enc_dir_1 = reshape(cell2mat(temp_enc_1),length(sourceGrangerEnc_group),length(freq));
dataAll_subj_enc_dir_2 = reshape(cell2mat(temp_enc_2),length(sourceGrangerEnc_group),length(freq));

dataMaint{1}.powspctrm = dataAll_subj_maint_dir_1;
dataMaint{2}.powspctrm = dataAll_subj_maint_dir_2;
dataEnc{1}.powspctrm = dataAll_subj_enc_dir_1;
dataEnc{2}.powspctrm = dataAll_subj_enc_dir_2;


cfg = [];
cfg.event_comparisons = {[1 2]}
stat_maint = statistics_fr(cfg,dataMaint);
stat_enc = statistics_fr(cfg,dataEnc);

sgnf_bands_maint =stat_maint{1}.hipp_cortex_maint.freq(stat_maint{1}.hipp_cortex_maint.stat>30);
sgnf_bands_enc = stat_enc{1}.hipp_cortex_maint.freq(stat_enc{1}.hipp_cortex_maint.stat>35);

sgnf_bands_enc = unique(ceil(sgnf_bands_enc));
sgnf_bands_enc = sgnf_bands_enc(sgnf_bands_enc<=15);
sgnf_bands_maint = unique(ceil(sgnf_bands_maint));
sgnf_bands_maint =sgnf_bands_maint(sgnf_bands_maint<=15);

stat_maint_freq = find(stat_maint{1}.hipp_cortex_maint.stat>30);
stat_enc_freq =find(stat_enc{1}.hipp_cortex_maint.stat>35);