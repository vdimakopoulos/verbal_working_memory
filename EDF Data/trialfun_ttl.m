function [trl, event] = trialfun_ttl(cfg)

% read the header information
hdr           = ft_read_header(cfg.dataset);

% read the events from the data
chanindx      = 1;
detectflank   = 'up';
threshold     = '(1/2)*nanmedian'; % or, e.g., 1/2 times the median for down flanks
event         = ft_read_event(cfg.dataset, 'chanindx', chanindx, 'detectflank', detectflank, 'threshold', threshold);

% define trials around the events
trl           = [];
pretrig       = 2* hdr.Fs; % e.g., 1 sec before trigger
posttrig      = 6* hdr.Fs; % e.g., 2 sec after trigger
for i = 1:numel(event)
    offset    = -hdr.nSamplesPre;  % number of samples prior to the trigger
    trlbegin  = event(i).sample - pretrig;
    trlend    = event(i).sample + posttrig;
    newtrl    = [trlbegin trlend offset];
    trl       = [trl; newtrl]; % store in the trl matrix
end