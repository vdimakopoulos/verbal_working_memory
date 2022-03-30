function [statistics_tfrs] = statistics_fr(cfg, data_cell)
% This function calculates the statistics for the tfr events
%
% cfg must contain:
%
%   cfg.event_comparisons = a cell array containing the events numbers for
%   the events to be compared
%
% cfg can contain anything that ft_freqstatistics recognizes

event_comparisons = cfg.event_comparisons;
n_comparisons = length(event_comparisons);
tfrs = data_cell;
statistics_tfrs = [];
n_subjects = length(tfrs);
data_comparison = data_cell;
for comparison_index = 1:n_comparisons
   cfg                  = [];
    cfg.method           = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability
    cfg.statistic        = 'indepsamplesT'; % use the independent samples T-statistic as a measure to
    % evaluate the effect at the sample level
    cfg.correctm         = 'fdr';
    % will be used for thresholding
    cfg.tail             = 0;          % -1, 1 or 0 (default = 0); one-sided or two-sided test
    cfg.alpha            = 0.05;      % alpha level of the permutation test
    cfg.numrandomization = 2500;
%     cfg.neighbours = [1 1];
    cfg.minnbchan        = 2;          % minimum number of neighborhood channels that is
    
    subj = 15;
    design = zeros(2,2*subj);
    for i = 1:subj
        design(1,i) = i;
    end
    for i = 1:subj
        design(1,subj+i) = i;
    end
    design(2,1:subj)        = 1;
    design(2,subj+1:2*subj) = 2;
    cfg.design = design;
    cfg.ivar = [2];
%     cfg.uvar = 1;
    this_stat = ft_freqstatistics(cfg, data_comparison{1},data_comparison{2});
    this_stat = rmfield(this_stat, 'cfg');
    statistics_tfrs = setfield(statistics_tfrs, 'hipp_cortex_maint', this_stat);
end

statistics_tfrs = {statistics_tfrs}; %% return as cell