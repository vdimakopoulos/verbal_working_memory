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
%% load EEG data for a patient
pID = 4;
[dataScalp,TrialInformationTable,dataScalp_Task_SS] = LoadEEGData(pID,Drive_Letter)
% dataScalp = dataScalp{4};
cd(strPaths.Main)

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

%% read electrodes
elec = ft_read_sens('standard_1020.elc');

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

%% keep only the scalp channels that exist in the dataset
chansEEG = dataScalp.label;
ind_elc = find(ismember(elec_aligned.label,chansEEG));
elec_aligned.chanpos = elec_aligned.chanpos(ind_elc,:);
elec_aligned.chantype = elec_aligned.chantype(ind_elc,:);
elec_aligned.chanunit = elec_aligned.chanunit(ind_elc,:);
elec_aligned.elecpos = elec_aligned.elecpos(ind_elc,:);
elec_aligned.label = elec_aligned.label(ind_elc,:);
elec_aligned = rmfield(elec_aligned,'tra')
%% visualize electrodes
figure;
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
cfg.resolution       = 1;   
cfg.sourcemodel.unit = 'cm';
sourcemodel = ft_prepare_leadfield(cfg);
sourcemodel = ft_convert_units(sourcemodel,'mm')

inside = sourcemodel;
outside = sourcemodel;

inside.pos = sourcemodel.pos(sourcemodel.inside, :);
outside.pos = sourcemodel.pos(~sourcemodel.inside, :);


figure
hold on
ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
ft_plot_mesh(outside, 'vertexsize', 20);
ft_plot_headmodel(vol, 'facealpha', 0.1)

view(125, 10)

%% leadfield
cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = vol;      %% how do currents spread?
cfg.elec        = elec_aligned; %% where are the sensors?

% how do sources and sensors connect?
sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);

plot_leadfield(elec_aligned,sourcemodel_and_leadfield,vol,vol.bnd(3))

%% save headmodel source model leadfield and electrodes
save(['F:/Vasileios/Task Analysis/Data/Analysis Data/SourceLocalization/',date,'_patient_',int2str(28), '_SourceLocalization_vars'],'sourcemodel_and_leadfield','sourcemodel','vol','elec_aligned','mri_realigned','mesh','-v7.3')


%% prepare the data for beamforming
if pID == 42
    cfg = [];
    cfg.channel = {'all','-T1','-T2'}
    dataScalp = ft_preprocessing(cfg,dataScalp)
    [dataScalp_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable);
    Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
    [dataScalp_SS,TrialInformationTable_SS,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataScalp_SS,TrialInformationTable);
    
    %maintenance
    nSet_Size = size(dataScalp_SS,2);
    dataScalp_Task_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [-2,-1/dataScalp_SS{iSS}.fsample];
        
        dataScalp_Task_SS{iSS} = ft_selectdata(cfg,dataScalp_SS{iSS});
    end

else
    dataScalp_Task_SS = dataScalp_Task_SS{3}
end
%% ICA
cfg = [];
cfg.method       = 'runica'
cfg.channel      = {'all'};
cfg.numcomponent = 'all';
cfg.demean       = 'yes';
cfg.feedback     =  'text'
IC_components    = ft_componentanalysis(cfg,dataScalp)%dataScalp_Task_SS{5})

%% Plot the IC components of ICA using topoplot
figure;
cfg = [];
cfg.component  = [1:size(IC_components.topo,1)];
cfg.layout      = 'elec1020.lay';
cfg.colormap    = 'jet';
cfg.colorbar    = 'no'
cfg.comment = 'no';
ft_topoplotIC(cfg,IC_components);

%%
cfg = [];
cfg.layout = 'elec1020.lay';
cfg.viewmode = 'component';

ft_databrowser(cfg, IC_components);
%% reject comp
cfg = [];
cfg.component = [1];
cfg.demean = 'yes';
dataBipolarScalp_Reconstructed = ft_rejectcomponent(cfg,IC_components);

cfg = [];
dataScalp2 = ft_selectdata(cfg,dataBipolarScalp_Reconstructed)
%% reject trials with artifacts
cfg          = [];
cfg.method   = 'summary';
cleanEEGdata        = ft_rejectvisual(cfg,dataScalp)%dataScalp_Task_SS{3}{4});
%% transform data

cfg                  = [];
cfg.covariance       = 'yes';
timelock             = ft_timelockanalysis(cfg, dataScalp2)%dataScalp_Task_SS{3}{4});



%% beamforming lcmv
% for i = 1:length(timelock.label)
cfg              = [];
cfg.method       = 'sloreta';
cfg.sourcemodel  = sourcemodel_and_leadfield;
cfg.headmodel    = vol;
cfg.elec = elec_aligned;
% cfg.keeepfilter = 'yes'
% cfg.sourcemodel.unit = 'mm';
% cfg.keepfilter       = 'yes';
% cfg.channel = timelock.label{i};
source = ft_sourceanalysis(cfg,dataScalp_Task_SS{3}{4});


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



%% construct a virtual channel in the scalp data containing only the source reconstructed
% Find the maximum source reconstructed from the beamformer
data = cleanEEGdata%dataScalp_Task_SS{3}{4};
[maxval, maxindx] = max(source.avg.pow);
maxpos = source.pos(maxindx,:);

% Find which spatial beamformer filter is in that position
x = find(source.pos(:,1) == maxpos(:,1));
y = find(source.pos(:,2) == maxpos(:,2));
z = find(source.pos(:,3) == maxpos(:,3));
filter_pos = intersect(intersect(x,y),z);

beamformer = source.avg.filter{filter_pos};

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
  virtualchanneldata.trial{i} = 10^(-5)*u(:,1)' * beamformer * data.trial{i}(chansel,:);
end


cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, data);

cfg = [];
cfg.viewmode = 'vertical'; 
ft_databrowser(cfg, combineddata);
