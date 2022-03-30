function [Sgranger] = FieldTripSpectralCGC(x, FREQINTERVAL, Fs)
% Clean up old FieldTrip Granger causalist code - RL Jenison 4/16/2020

[Ntrials,Nchannels,Ntime] = size(x);

t = [1:Ntime]/Fs;

clear data
for i = 1:Ntrials
    
    data.trial{i} = squeeze(x(i,:,:));
    data.time{i} = t;
end
label = ['a' 'b' 'c' 'd' 'e' 'f' 'g' 'h' 'i' 'j' 'k' 'l' 'm' 'n' 'o' 'p' 'q' 'r' ...
    's' 't' 'u' 'v' 'w' 'x' 'y' 'z' 'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' ...
    'I' 'J' 'K' 'L' 'M' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];
% label = strLabel;
for i = 1:Nchannels
    data.label{i} = label(i);
end

data.fsample = Fs;

data.cfg.ntrials     = Ntrials;
data.cfg.triallength = Ntime/Fs;
data.cfg.fsample     = Fs;
data.cfg.nsignal     = Nchannels;

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foi    = FREQINTERVAL;
cfg.tapsmofrq = 4;
cfg.keeptrials = 'yes';
cfg.keeptapers = 'yes';
cfg.demean                    = 'yes';
ft_warning off  
freq          = ft_freqanalysis(cfg, data);

cfg           = [];
cfg.method    = 'granger';
cfg.granger.conditional = 'no';

Sgranger       = ft_connectivityanalysis(cfg, freq);


