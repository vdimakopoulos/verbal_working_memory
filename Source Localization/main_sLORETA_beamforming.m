%% Paths
Drive_Letter = 'F:\';
strPaths.Main = [Drive_Letter,'Vasileios\'];
strPaths.Project = [strPaths.Main, 'Task Analysis\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = [Drive_Letter,'Vasileios\Toolboxes\fieldtrip-20200315\'];

% EEGLAB toolbox
strPaths.Toolboxes.EEGLAB               = [Drive_Letter,'Vasileios\Toolboxes\eeglab14_1_1b\'];

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(strPaths.Main)
addpath(genpath(strPaths.Project))
addpath(strPaths.Toolboxes.FieldTrip)
% Remove EEGLAB from path
rmpath(genpath(strPaths.Toolboxes.EEGLAB))

%Add figure tools on toolbar
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

% open meeg
setenv('PATH','F:\Vasileios\Toolboxes\OpenMEEG\bin\')
setenv('LD_LIBRARY_PATH','F:\Vasileios\Toolboxes\OpenMEEG\lib\')
setenv('TMP','F:/Vasileios/Temp/')

flag_patient_specific = 0;
%% load EEG data for a patient
pID = 4;
[dataScalp,TrialInformationTable,dataScalp_Task_SS] = LoadEEGData(pID,Drive_Letter)
% dataScalp = dataScalp{4};
cd(strPaths.Main)

%% Patient Specific headmodel
if flag_patient_specific
    %% Read MRI
    mri = loadPatientMRI(pID)
    
    %% align MRI to CTF
    cfg = [];
    cfg.coordsys = 'ctf';
    cfg.method = 'interactive';
    cfg.viewresult = 'yes';
    mri_realigned = ft_volumerealign(cfg,mri);
    
    %% segment head compartments
    cfg           = [];
    cfg.output    = {'brain','skull','scalp'};
    % cfg.scalpthreshold = 0.2
    segmentedmri  = ft_volumesegment(cfg, mri_realigned);
    
    %% meshing the compartments
    cfg=[];
    % cfg.method = 'iso2mesh';
    cfg.tissue={'brain','skull','scalp'};
    cfg.numvertices = [3000 2000 1000];
    mesh = ft_prepare_mesh(cfg,segmentedmri);
    
    %% volume conduction model
    cfg        = [];
    cfg.method = 'openmeeg';
    cfg.tissue= {'brain','skull','scalp'};
    vol        = ft_prepare_headmodel(cfg,segmentedmri);
    
    %% visualization
    figure;
    % brain
    b = ft_plot_mesh(vol.bnd(1), 'facecolor',[0.3 0.6 0.3], 'facealpha', 0.3, 'edgecolor', 'none', 'edgealpha', 0.05);
    hold on;
    % skull
    sk = ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4,'facecolor',[0.8 0.8 0.8]);
    hold on;
    % scalp
    sc =ft_plot_mesh(vol.bnd(3),'edgecolor','none','facecolor',[0.2 0.2 0.2],'facealpha',0.3);
    % sagitall view
    set(gca,'view',[-1.659371794480957e+02,20.967737453644254])
    legend('brain','skull','scalp')
else
    %% Load template headmodel
    load('F:\Vasileios\Toolboxes\fieldtrip-20200315\template\headmodel\headmodel_dipoli.mat')
    vol = headmodel_dipoli;
    vol = ft_convert_units(vol,'mm')
    %% visualization
    figure;
    % brain
    b = ft_plot_mesh(vol.bnd(3), 'facecolor',[0.3 0.6 0.3], 'facealpha', 0.3, 'edgecolor', 'none', 'edgealpha', 0.05);
    hold on;
    % skull
    sk = ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4,'facecolor',[0.8 0.8 0.8]);
    hold on;
    % scalp
    sc =ft_plot_mesh(vol.bnd(1),'edgecolor','none','facecolor',[0.2 0.2 0.2],'facealpha',0.3);
    % sagitall view
    set(gca,'view',[-1.659371794480957e+02,20.967737453644254])
    legend('brain','skull','scalp')

end

%% read electrodes
elec = ft_read_sens('standard_1020.elc');

if flag_patient_specific
    %% realign electrodes
    nas = mri_realigned.cfg.fiducial.nas;
    lpa = mri_realigned.cfg.fiducial.lpa;
    rpa = mri_realigned.cfg.fiducial.rpa;
    
    transm = mri_realigned.transform;
    
    nas = ft_warp_apply(transm,nas, 'homogenous');
    lpa = ft_warp_apply(transm,lpa, 'homogenous');
    rpa = ft_warp_apply(transm,rpa, 'homogenous');
    
    fid.elecpos       = [nas; lpa; rpa];       % ctf-coordinates of fiducials
    fid.label         = {'Nz','LPA','RPA'};    % same labels as in elec
    fid.unit          = 'mm';                  % same units as mri
    
    % alignment
    cfg               = [];
    cfg.method        = 'fiducial';
    cfg.target        = fid;                   % see above
    cfg.elec          = elec;
    cfg.fiducial      = {'Nz', 'LPA', 'RPA'};  % labels of fiducials in fid and in elec
    elec_aligned      = ft_electroderealign(cfg);
else
    elec_aligned = elec;
end
%% keep only the scalp channels that exist in the dataset
chansEEG = dataScalp.label;
ind_elc = find(ismember(elec_aligned.label,chansEEG));
elec_aligned.chanpos = elec_aligned.chanpos(ind_elc,:);
elec_aligned.chantype = elec_aligned.chantype(ind_elc,:);
elec_aligned.chanunit = elec_aligned.chanunit(ind_elc,:);
elec_aligned.elecpos = elec_aligned.elecpos(ind_elc,:);
elec_aligned.label = elec_aligned.label(ind_elc,:);
if flag_patient_specific
    elec_aligned = rmfield(elec_aligned,'tra')
end
%% visualize electrodes
figure;
% elec_aligned = ft_convert_units(elec_aligned,'m');
% head surface (scalp)
ft_plot_mesh(vol.bnd(1), 'edgecolor','none','facealpha',0.8,'facecolor',[0.5 0.5 0.5]);
hold on;
% electrodes
ft_plot_sens(elec_aligned,'style', 'k');

%% realign electrodes again
cfg           = [];
cfg.method    = 'interactive';
cfg.elec      = elec_aligned;
cfg.headshape = vol.bnd(1);
elec_aligned  = ft_electroderealign(cfg);


%% source model

cfg                  = [];
cfg.elec             = elec_aligned;
cfg.headmodel        = vol;

% use a 3-D grid with a 1 cm resolution
cfg.resolution       = 2;   
cfg.sourcemodel.unit = 'cm';
cfg.inwardshift  = -1;
sourcemodel = ft_prepare_leadfield(cfg);
sourcemodel = ft_convert_units(sourcemodel,'mm')

inside = sourcemodel;
outside = sourcemodel;

inside.pos = sourcemodel.pos(sourcemodel.inside, :);
inside.inside = inside.inside(sourcemodel.inside, :)
inside.leadfield = inside.leadfield([sourcemodel.inside])
outside.pos = sourcemodel.pos(~sourcemodel.inside, :);


figure
hold on
ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
ft_plot_mesh(outside, 'vertexsize', 20);
ft_plot_headmodel(vol, 'facealpha', 0.1)

view(125, 10)

%% source model 2
cfg = [];
cfg.elec = elec_aligned;
cfg.headmodel = vol;
cfg.tight = 'yes'
cfg.inwardshift  = -1.5;
cfg.resolution       = 3;   
cfg.sourcemodel.unit = 'cm';
template_grid    = ft_prepare_sourcemodel(cfg);
figure
hold on
ft_plot_mesh(template_grid, 'vertexsize', 20, 'vertexcolor', 'red');
ft_plot_headmodel(vol, 'facealpha', 0.1)

view(125, 10)
% sourcemodel = template_grid;
%% leadfield
%% convert all models in SI Units
sourcemodel = ft_convert_units(sourcemodel,'m');
elec_aligned = ft_convert_units(elec_aligned,'m');
vol = ft_convert_units(vol,'m');

cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = vol;      %% how do currents spread?
cfg.elec        = elec_aligned; %% where are the sensors?

% how do sources and sensors connect?
sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);
sourcemodel_and_leadfield = ft_convert_units(sourcemodel_and_leadfield,'m')
plot_leadfield(elec_aligned,sourcemodel_and_leadfield,vol,vol.bnd(3))

%% save headmodel source model leadfield and electrodes
save(['F:/Vasileios/Task Analysis/Data/Analysis Data/SourceLocalization/',date,'_patient_',int2str(28), '_SourceLocalization_vars'],'sourcemodel_and_leadfield','sourcemodel','vol','elec_aligned','mri_realigned','mesh','-v7.3')

save(['F:/Vasileios/Task Analysis/Data/Analysis Data/SourceLocalization/',date,'_patient_',int2str(30), '_SourceLocalization_vars'],'sourcemodel_and_leadfield','sourcemodel','vol','elec_aligned','-v7.3')



%% select task data
TaskData = dataScalp_Task_SS{2}{4};
%% reject trials with artifacts

cfg          = [];
cfg.method   = 'summary';
cleanEEGdata        = ft_rejectvisual(cfg,TaskData)%dataScalp_Task_SS{3}{4});
%% transform data

cfg                  = [];
cfg.covariance       = 'yes';
timelock             = ft_timelockanalysis(cfg, cleanEEGdata)%dataScalp_Task_SS{3}{4});



%% sLORETA
% for i = 1:length(timelock.label)
cfg              = [];
cfg.method       = 'sloreta';
cfg.sourcemodel  = sourcemodel_and_leadfield%   sourcemodel_and_leadfield;
cfg.headmodel    = vol;
cfg.sloreta.lambda = 25;
cfg.elec = elec_aligned;
source = ft_sourceanalysis(cfg,cleanEEGdata);


cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
source_Intrp  = ft_sourceinterpolate(cfg, source ,vol.bnd(3));

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'pow';
% cfg.atlas = 'BNA_MPM_thr25_1.25mm.nii'
ft_sourceplot(cfg,source_Intrp);
% end

%% another way of visualization
cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'pow';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.surffile       = 'surface_pial_both.mat'
cfg.surfdownsample = 10;  % downsample to speed up processing
ft_sourceplot(cfg, source);

%% data for beamforming
cfg = [];
cfg.resamplefs = 160;
dataScalp_Task_SS{1}{4} = ft_resampledata(cfg,dataScalp_Task_SS{1}{4});

cfg = [];
cfg.parameter = 'trial';
cfg.operation = '(x1-x2)'
dataLCMV = ft_math(cfg,dataScalp_Task_SS{2}{4},dataScalp_Task_SS{1}{4})
dataLCMV2 = ft_math(cfg,dataScalp_Task_SS{3}{4},dataScalp_Task_SS{1}{4})

%% PCC beamforming
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = [10];
freq           = ft_freqanalysis(cfg, dataLCMV);

%% beamforming
cfg              = [];
cfg.method       = 'lcmv';
cfg.sourcemodel  = sourcemodel_and_leadfield;
cfg.headmodel    = vol;
cfg.elec = elec_aligned;
cfg.keeepfilter = 'yes'
cfg.keepfilter       = 'yes';
cfg.lcmv.projectmom = 'yes'; % project dipole time series in direction of maximal power (see below)
cfg.lcmv.kurtosis = 'yes'; % compute kurtosis at each location
% cfg.channel = {'Cz'};
cfg.atlas = 'BNA_MPM_thr25_1.25mm.nii'
cfg.inputcoord = 'mni';
sourceLCMV = ft_sourceanalysis(cfg,dataScalp_Task_SS{2}{4});
sourceLCMV_bsln = ft_sourceanalysis(cfg,dataScalp_Task_SS{1}{4});
sourceLCMV2 = ft_sourceanalysis(cfg,dataLCMV2);

sourceLCMV_bsln.avg.pow(find(isnan(sourceLCMV_bsln.avg.pow))) = 0.1;
s1 = sourceLCMV;
s2 = sourceLCMV2;
cfg = [];
cfg.parameter = 'pow'
cfg.operation = '((x1-x2)./x2)*100'
s1 = ft_math(cfg,sourceLCMV,sourceLCMV_bsln);cfg = [];
cfg = [];
cfg.parameter = 'pow'
cfg.operation = '((x1-x2)./x2)*100'
s2 = ft_math(cfg,sourceLCMV2,sourceLCMV_bsln);
% cfg                   = [];
% cfg.frequency         = freq.freq;
% cfg.elec = elec_aligned;
% cfg.method            = 'pcc';
% cfg.grid              = sourcemodel_and_leadfield;
% cfg.headmodel         = vol;
% cfg.keeptrials        = 'yes';
% cfg.pcc.lambda        = '10%';
% cfg.pcc.projectnoise  = 'yes';
% cfg.pcc.keepfilter    = 'yes';
% cfg.pcc.fixedori      = 'yes';
% sourceFreq = ft_sourceanalysis(cfg, freq);
% sourceFreq  = ft_sourcedescriptives([], sourceFreq);


cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'kurtosis';
source_Intrp  = ft_sourceinterpolate(cfg, s1 ,vol.bnd(3));

cfg              = [];
cfg.method       = 'surface';
cfg.funparameter = 'kurtosis';
ft_sourceplot(cfg,source_Intrp);
% 
% % find temporal regions
% aal = ft_read_atlas('F:\Vasileios\Toolboxes\fieldtrip-20200315\template\atlas\aal\ROI_MNI_V4.nii');
% temp_idx = strmatch('Temporal', aal.tissuelabel);
% 
% cfg = [];
% cfg.inputcoord = 'mni';
% cfg.atlas = aal;
% cfg.roi = aal.tissuelabel(temp_idx);
% mask = ft_volumelookup(cfg, sourceLCMV);
%% another way of visualization
% atlas = ft_read_atlas('ROI_MNI_V4.nii');
% sourceLCMV = ft_determine_coordsys(sourceLCMV);

cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'kurtosis';
cfg.maskparameter  = cfg.funparameter;
cfg.funcolormap    = 'jet';
cfg.surffile       = 'surface_pial_both.mat'
cfg.surfdownsample = 10;  % downsample to speed up processing
% cfg.funcolorlim = [0 50]
% cfg.atlas = atlas;
% cfg.roi = atlas.tissuelabel(temp_idx);
% cfg.inputcoord = 'mni';

ft_sourceplot(cfg, s1);


%% construct a virtual channel in the scalp data containing only the source reconstructed
% Find the maximum source reconstructed from the beamformer
data = dataScalp_Task_SS{2}{4};
data2 = dataScalp_Task_SS{3}{4};
[maxval, maxindx] = max(sourceLCMV.avg.kurtosis);
maxpos = sourceLCMV.pos(maxindx,:);

% Find which spatial beamformer filter is in that position
x = find(sourceLCMV.pos(:,1) == maxpos(:,1));
y = find(sourceLCMV.pos(:,2) == maxpos(:,2));
z = find(sourceLCMV.pos(:,3) == maxpos(:,3));
filter_pos = intersect(intersect(x,y),z);

beamformer = sourceLCMV.avg.filter{filter_pos};

chansel = ft_channelselection('EEG', data.label); % find the names
chansel = match_str(data.label, chansel);         % find the indices

sourcedata = [];
sourcedata.label = {'x', 'y', 'z'};
sourcedata.time = data.time;
for i=1:length(data.trial)
  sourcedata.trial{i} = beamformer * data.trial{i}(chansel,:);
end


timeseries = cat(2, sourcedata.trial{:});

[u, s, v] = svd(timeseries, 'econ');

timeseriesmaxproj = u(:,1)' * timeseries;

virtualchanneldata = [];
virtualchanneldata.label = {'cortex'};
virtualchanneldata.time = data.time;
for i=1:length(data.trial)
  virtualchanneldata.trial{i} = 10^2*u(:,1)' * beamformer * data.trial{i}(chansel,:);
end


cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, data);

cfg = [];
cfg.viewmode = 'vertical'; 
ft_databrowser(cfg, combineddata);
%%
data2 = dataScalp_Task_SS{3}{4};
[maxval, maxindx] = max(sourceLCMV2.avg.kurtosis);
maxpos = sourceLCMV2.pos(maxindx,:);

% Find which spatial beamformer filter is in that position
x = find(sourceLCMV2.pos(:,1) == maxpos(:,1));
y = find(sourceLCMV2.pos(:,2) == maxpos(:,2));
z = find(sourceLCMV2.pos(:,3) == maxpos(:,3));
filter_pos = intersect(intersect(x,y),z);

beamformer = sourceLCMV2.avg.filter{filter_pos};

chansel = ft_channelselection('EEG', data.label); % find the names
chansel = match_str(data2.label, chansel);         % find the indices

sourcedata = [];
sourcedata.label = {'x', 'y', 'z'};
sourcedata.time = data2.time;
for i=1:length(data2.trial)
  sourcedata.trial{i} = beamformer * data2.trial{i}(chansel,:);
end


timeseries2 = cat(2, sourcedata.trial{:});

[u, s, v] = svd(timeseries2, 'econ');

timeseriesmaxproj = u(:,1)' * timeseries2;

virtualchanneldata = [];
virtualchanneldata.label = {'cortex'};
virtualchanneldata.time = data2.time;
for i=1:length(data2.trial)
  virtualchanneldata.trial{i} = 10^2*u(:,1)' * beamformer * data2.trial{i}(chansel,:);
end


cfg = [];
combineddata2 = ft_appenddata(cfg, virtualchanneldata, data2);


cfg = [];
cfg.viewmode = 'vertical'; 
ft_databrowser(cfg, combineddata2);
%% Parameters for loading the hipp data
analysis_type = 'correct trials'; %'incorrect trials'; %or 'correct trials'
% Granger Causality calculation between a hippocampal contact of a depth 
%electrode and a contact in the cortex of a depth electrode (outer contacts)
flag_Hipp_GC_Ctx_Depth = 1; 
%% load hippocampal data

strDataPath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Dataset_2020_MTL_iEEG_scalpEEG\';
strDataPath_Macro = [strDataPath,'Macro\'];
strDataPath_Scalp = [strDataPath,'Scalp\'];
nSubjects = [1:9];
strPat_IDs = {'WC','SW','IA','HJ','PL','AM','HF','MA','MS'};
nSubject_ID = [28 22 19 30 33 13 23 29 16];
% specify the patient's data directory
selected_subj =  pID;
strPatientDataPath_Macro = sprintf('%s%d %s\\',strDataPath_Macro,nSubject_ID(selected_subj),strPat_IDs{selected_subj});
strPatientDataPath_Scalp = sprintf('%s%d %s\\',strDataPath_Scalp,nSubject_ID(selected_subj),strPat_IDs{selected_subj});


[dataBipolar, dataBipolar_Scalp, TrialInformationTable] = ...
    Load_Dataset_from_human_MTL(strPatientDataPath_Macro,strPatientDataPath_Scalp);

cd(strPaths.Main)

%% Reject Artifactual trials
[dataBipolar_clean, TrialInformationTable_iEEG_clean ] = ...
    Get_Only_Clean_Trials_FieldTrip(dataBipolar, TrialInformationTable.Macro);
%% Rereference
% Load anatomical locations for the current patient
strAnatLocationsPath = 'F:\Vasileios\Task Analysis\Data\Sternberg Task\Dataset_2020_MTL_iEEG_scalpEEG\Anatomical Locations\';
strFilePrefix = 'Anatomical_Locations_Patient_';
strSubjAnatLocationFile = sprintf('%s%s%d.mat',strAnatLocationsPath,strFilePrefix,nSubject_ID(selected_subj));
load(strSubjAnatLocationFile);

% Find white matter contacts
nWM_contacts = find(strcmp(AnatomicalLocations.Anatomical_Location,'no_label_found'))
strWM_labels = AnatomicalLocations.Labels(nWM_contacts);
hipp_chans_WM_contacts = find(contains(strWM_labels,'H'));
refChannel_iEEG = nWM_contacts(hipp_chans_WM_contacts(1));
% Rereferencing iEEG
cfg               =   [];
cfg.reref         =  'yes'
cfg.refchannel    =   refChannel_iEEG;%white matter contacts referencing
cfg.refmethod     =  'avg';
dataBipolar_reref       = ft_preprocessing(cfg,dataBipolar_clean);
dataBipolar_reref.label = strrep(dataBipolar.label,'m','');

%% Apply montage
[montage,chans_to_be_bipolar,dataBip,Bip_chans] = GetDataInBipolar(dataBipolar_reref,selected_subj,flag_Hipp_GC_Ctx_Depth);

%% Downsample to 200 Hz
cfg = [];
cfg.resamplefs = 40;%200;
dataBipolarResampled = ft_resampledata(cfg,dataBip);
combineddata = ft_resampledata(cfg,combineddata);
combineddata2 = ft_resampledata(cfg,combineddata2);
%% macro data
% Select only correct trials
switch analysis_type
    case 'correct trials'
        [dataBipolar_SS,TrialInformationTableCorrect] = Get_Only_Correct_Trials_FieldTrip(dataBipolarResampled,TrialInformationTable_iEEG_clean);
    case 'incorrect trials'
        [dataBipolar_SS, TrialInformationTableCorrect] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolarResampled, TrialInformationTable_iEEG_clean )
end

%% Divide data into set sizes
%   [4] -> Set Size 1
%       [6] -> Set Size 2
%           [8] -> Set Size 3
%               [6 8] -> Set Size 4
%                   [4 6 8] -> Set Size 5
Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
TrialInformationTableCorrect.SetSize = TrialInformationTableCorrect.Setsize;
[dataBipolar_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTableCorrect);

%% Select the latencies for every task period for each task period
[dataBipolar_Ret_SS{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
[dataBipolar_Ret_SS{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
[dataBipolar_Ret_SS{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance

%% Append the hippocampal with source data
combineddata.time = dataBipolar_Ret_SS{2}{4}.time;
combineddata2.time = dataBipolar_Ret_SS{3}{4}.time;
cfg = [];
dataEnc = ft_appenddata(cfg,dataBipolar_Ret_SS{2}{4},combineddata);
dataMaint = ft_appenddata(cfg,dataBipolar_Ret_SS{3}{4},combineddata2);


%% Granger
strChannelNameList = dataEnc.label;
iSS = 4;
iPair = 1;
freq                   = [];
freq.freqcfg           = [];
freq.freqcfg.method    = 'mtmfft';
freq.freqcfg.output    = 'fourier';
freq.freqcfg.taper     = 'hanning';
freq.freqcfg.pad       = 20;
freq.freqcfg.tapsmofrq = 2;
Maintenance_freq       = ft_freqanalysis(freq.freqcfg, dataMaint);
Encoding_freq          = ft_freqanalysis(freq.freqcfg, dataEnc);

grangercfg = [];
grangercfg.method  = 'granger';
grangercfg.granger.conditional = 'no';
grangercfg.granger.sfmethod = 'bivariate';
grngChan1 = 'AHR1-AHR3';
grngChan2 = 'cortex';
grangercfg.channelcmb = {grngChan1,grngChan2};
gdata.Maint{iPair,iSS}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
gdata.Enc{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Encoding_freq);

%% Visualization

figure;

light_blue = [0.30,0.75,0.93];
light_red       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',light_red};

hold on;
freq = Maintenance_freq.freq;
semilogx(freq,gdata.Enc{iPair,iSS}.grangerspctrm(2,:),'LineWidth',3,'color',Colors{1});
semilogx(freq,gdata.Enc{iPair,iSS}.grangerspctrm(1,:),'LineWidth',3,'color',Colors{2})
semilogx(freq,gdata.Maint{iPair,iSS}.grangerspctrm(2,:),'LineWidth',3,'color',Colors{3})
semilogx(freq,gdata.Maint{iPair,iSS}.grangerspctrm(1,:),'LineWidth',3,'color',Colors{4})