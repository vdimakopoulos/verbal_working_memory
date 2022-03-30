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
    vol = ft_convert_units(vol,'mm');
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
%% convert all models in SI Units
sourcemodel = ft_convert_units(sourcemodel,'m');
vol = ft_convert_units(vol,'m');
elec_aligned = ft_convert_units(elec_aligned,'m');

cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = vol;      %% how do currents spread?
cfg.elec        = elec_aligned; %% where are the sensors?

% how do sources and sensors connect?
sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);

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
cfg.sourcemodel  = sourcemodel_and_leadfield;
cfg.headmodel    = vol;
cfg.sloreta.lambda = 2;
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


%% beamforming
cfg = [];
cfg.parameter = 'trial';
cfg.operation = '(x1-x2)/x2'
dataLCMV = ft_math(cfg,dataScalp_Task_SS{2}{4},dataScalp_Task_SS{1}{4})

cfg              = [];
cfg.method       = 'lcmv';
cfg.sourcemodel  = sourcemodel_and_leadfield;
cfg.headmodel    = vol;
cfg.elec = elec_aligned;
cfg.keeepfilter = 'yes'
cfg.keepfilter       = 'yes';
% cfg.channel = timelock.label{i};
sourceLCMV = ft_sourceanalysis(cfg,dataLCMV);


cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
source_Intrp  = ft_sourceinterpolate(cfg, sourceLCMV ,vol.bnd(3));

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
cfg.funcolorlim = [0 1e-2]
ft_sourceplot(cfg, sourceLCMV);
