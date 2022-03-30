
% Fill time-frequency time-bins with spectral Granger causality
% Input time-series matrix [TRIALS x CHANNELS x TIME]

Ntime = size(x,3);
Nchannels = size(x,2);
INTERVAL = 100;
STEP = 1;
Fs = 500;
FREQINTERVAL = [4:1:100];
label = dataBipolar.label;
indx = slideWindow(1, Ntime, INTERVAL, STEP)';
M = length(indx);
clear GStime
for i = 1:M
    
    xtf = x(:,:,indx(i,1):indx(i,2));
   
 
    [Sgranger] = FieldTripSpectralCGC(xtf,FREQINTERVAL,Fs);
 
    GStime(i,:,:,:) = Sgranger.grangerspctrm;


end