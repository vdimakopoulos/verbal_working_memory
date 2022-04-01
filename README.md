# Analysis of information flow during a verbal working memory task

## The task
In the task, sets of consonants are presented and had to be memorized. The set size (4, 6, or 8 letters) determines WM workload. In each trial, presentation of a letter string (encoding period, 2 s) is followed by a delay (maintenance period, 3 s) where participants activate their phonological loop. After the delay, a probe letter is shown, and participants indicate whether the probe was presented during the encoding period (Boran et al., 2019a).

## Analysis steps
##### 1. Data Preprocessing
├── Appy bipolar montage  in  LFP channels

├── Rereference  hippocampal LFP and temporal ECoG  using separate iEEG contacts

|──  Resample the data to 500 Hz

|──  Reject large unitary artifacts in scalp EEG recordings

##### 2. Power Spectral Density (PSD)
├── Multitaper time-frequency transformation

├── PSD computation during encoding and maintenance  in hippocampal LFP, temporal ECoG and scalp EEG channels

├── Baseline against fixation 

##### 3. Functional Connectivity
├──  Phase locking value (PLV) calculation between  all pairs of hippocampal LFP and ECoG grid channels.

├──── Multitaper frequency transformation with 2 tapers using FFT,  (frequency resolution = 1 Hz)

##### 4. Source reconstruction of EEG using beamforming
├──  Forward problem 

├──── > Precomputed BEM head-model with 3 compartments (scalp,skull,brain)

├──── > EEG electrodes alignment to the scalp compartment

├──── > Determine  source grid locations acoording to AAL parcel

├──── > Source model and leadfield computation

├── Inverse problem 

├──── > Linearly constrained minimum variance beamformers

├──── > Reconstruction separate for fixation/encoding /maintenance

├──── > Baseline encoding/maintenance sources against fixation sources

##### 5. Spectral Granger causality (GC)
Evaluation of GC in the range [4 30] Hz
├── > Downsample the signals to 2*Nyquist frequency = 60 Hz

├── > Frequency transformation using multitaper FFT method with 2 Hann tapers, 20 seconds of padding

├── >  Non parametric spectral factorization  method for GC

├── >  GC between  hippocampal LFP and beamforming EEG sources 

##### 6. Statistics
├── >  Cluster-based nonparametric permutation tests

├── >  GC differences comparison to null distribution of differences

├── >  To test the statistical significance of the  spatial spread of contacts with high PSD/PLV/GC we computed the spatial collinearity on grid contacts agains null distribution

├── > Significant  beamforming EEG sources during encoding and maintenance are tested against fixation  with non-parametric permutation t-test

├── >  Group statistics to assess the direction of information flow  using group cluster based permutation t-test from FieldTrip

### Research Papers
[1] Persistent hippocampal neural firing and hippocampal-cortical coupling predict verbal working memory load
Boran E, Klaver P, Hilfiker P, Stieglitz L, Grunwald T, Sarnthein J, Science Advances 2019, Vol.5, Issue 3, doi: https://doi.org/10.1126/sciadv.aav3687 
[2] Information flows from hippocampus to auditory cortex during replay of verbal working memory items, Vasileios Dimakopoulos, Pierre Mégevand, Lennart Stieglitz, Lukas Imbach, Johannes Sarnthein, bioRxiv 2021.03.11.434989; doi: https://doi.org/10.1101/2021.03.11.434989
#### This code has been written by 
Vasileios Dimakopoulos (Vasileios.Dimakopoulos@usz.ch)
#### Parts of the preprocessing and statistics have been implemented initially by
Ece Boran (Ece.Boran@usz.ch)

#### For questions regarding the task contact
Johannes Sarnthein (Johannes.Sarnthein@usz.ch)
