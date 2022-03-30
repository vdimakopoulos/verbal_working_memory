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
warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState')
%% load EEG data for a patient
pID = 45;
[dataScalp,TrialInformationTable,dataScalp_Task_SS] = LoadEEGData(pID,Drive_Letter,'incorrect trials')
% dataScalp = dataScalp{4};
cd(strPaths.Main)
%% ----------------------
%% Preprocessing

%
cfg = [];
cfg.method = 'summary';
cfg.layout = 'standard_1020.elc';
% cfg.keeptrial = 'yes';
dataClean_fix = ft_rejectvisual(cfg,dataScalp_Task_SS{1}{4});
dataClean_enc = ft_rejectvisual(cfg,dataScalp_Task_SS{2}{4});
dataClean_maint = ft_rejectvisual(cfg,dataScalp_Task_SS{3}{4});

%ICA
cfg = [];
cfg.method = 'runica';
cfg.runica.maxsteps = 50;
comp_enc = ft_componentanalysis(cfg,dataClean_enc);
comp_maint = ft_componentanalysis(cfg,dataClean_maint);


cfg = [];
cfg.channel = {comp_enc.label{1:end}};
cfg.layout = 'standard_1020.elc';
cfg.compscale = 'local';
ft_databrowser(cfg,comp_enc)
ft_databrowser(cfg,comp_maint)


cfg = [];
cfg.component = [1];
dataICA = ft_rejectcomponent(cfg,comp_enc)
dataICA2 = ft_rejectcomponent(cfg,comp_maint)

%% timelock
cfg = [];
timelock = ft_timelockanalysis(cfg,dataClean_maint);

cfg = [];
cfg.showlabels = 'yes';
cfg.layout = 'standard_1020.elc';
figure;
ft_multiplotER(cfg,timelock);

%% ----------------------
%% Source analysis
%% Load template headmodel
load('F:\Vasileios\Toolboxes\fieldtrip-20200315\template\headmodel\headmodel_dipoli.mat')
vol = headmodel_dipoli;
vol = ft_convert_units(vol,'mm')
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
set(gca,'view',[-1.659371794480957e+02,20.967737453644254]);
legend('brain','skull','scalp');

%% read electrodes
elec = ft_read_sens('standard_1020.elc');
elec_aligned = elec;
%% keep only the scalp channels that exist in the dataset
chansEEG = dataScalp.label;
ind_elc = find(ismember(elec_aligned.label,chansEEG));
elec_aligned.chanpos = elec_aligned.chanpos(ind_elc,:);
elec_aligned.chantype = elec_aligned.chantype(ind_elc,:);
elec_aligned.chanunit = elec_aligned.chanunit(ind_elc,:);
elec_aligned.elecpos = elec_aligned.elecpos(ind_elc,:);
elec_aligned.label = elec_aligned.label(ind_elc,:);
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
template_grid = ft_prepare_leadfield(cfg);
template_grid = ft_convert_units(template_grid,'mm')

inside = template_grid;
outside = template_grid;

inside.pos = template_grid.pos(template_grid.inside, :);


inside.leadfield = inside.leadfield([template_grid.inside])
outside.pos = template_grid.pos(~template_grid.inside, :);


figure
hold on
ft_plot_mesh(inside, 'vertexsize', 20, 'vertexcolor', 'red');
ft_plot_mesh(outside, 'vertexsize', 20);
ft_plot_headmodel(vol, 'facealpha', 0.1)

view(125, 10)

%% Load atlas and create a binary mask
atlas = ft_read_atlas('F:\Vasileios\Toolboxes\fieldtrip-20200315\template\atlas\aal\ROI_MNI_V4.nii');
atlas = ft_convert_units(atlas,'mm'); % assure that atlas and template_grid are expressed in the %same units
template_grid.coordsys = 'mni';
cfg = []
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel;
cfg.inputcoord = 'mni';
mask = ft_volumelookup(cfg,template_grid);

% create temporary mask according to the atlas entries
tmp                  = repmat(template_grid.inside,1,1);
tmp(tmp==1)          = 0;
tmp(mask)            = 1;

% define inside locations according to the atlas based mask
template_grid.inside = tmp;

% plot the atlas based grid
figure; ft_plot_mesh(template_grid.pos(template_grid.inside,:));

%% recompute source model based on the atlas grid
load('F:\Vasileios\Task Analysis\template_mri.mat');
cfg                  = [];
cfg.warpmni = 'yes';
cfg.template = template_grid;
cfg.nonlinear = 'yes';
cfg.mri = mri;
sourcemodel = ft_prepare_sourcemodel(cfg);
% sourcemodel = ft_convert_units(sourcemodel,'mm')
%% Plot the final source model together with the individual head model and the sensor array

close all
hdm = ft_convert_units(vol, 'm');
sourcemodel = ft_convert_units(sourcemodel, 'm');
elec_aligned = ft_convert_units(elec_aligned, 'm');

figure; hold on
ft_plot_headmodel(hdm,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; %camlight;
alpha 0.4           % make the surface transparent

ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot only locations inside the volume

ft_plot_sens(elec_aligned,'style','*r'); % plot the sensor array
view ([0 -90 0])
%% leadfield
%% convert all models in SI Units
cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = hdm;      %% how do currents spread?
cfg.elec        = elec_aligned; %% where are the sensors?

% how do sources and sensors connect?
sourcemodel_and_leadfield = ft_prepare_leadfield(cfg);
sourcemodel_and_leadfield = ft_convert_units(sourcemodel_and_leadfield,'m')
plot_leadfield(elec_aligned,sourcemodel_and_leadfield,hdm,hdm.bnd(3))

%% save headmodel source model leadfield and electrodes
save(['F:/Vasileios/Task Analysis/Data/Analysis Data/SourceLocalization/',date,'_patient_',int2str(45), '_SourceLocalization_vars'],'sourcemodel_and_leadfield','sourcemodel','template_grid','hdm','elec_aligned','-v7.3')

%% lf 2
cfg                 = [];
cfg.channel         = dataScalp.label; % ensure that rejected sensors are not present
cfg.elec            = elec_aligned;
cfg.headmodel       = hdm;
cfg.lcmv.reducerank = 3; % default for MEG is 2, for EEG is 3
cfg.grid = sourcemodel;
[grid] = ft_prepare_leadfield(cfg);
plot_leadfield(elec_aligned,grid,hdm,hdm.bnd(3))

%% data for beamforming
cfg = [];
cfg.covariance = 'yes';
avgFix = ft_timelockanalysis(cfg,dataClean_fix);
avgEnc = ft_timelockanalysis(cfg,dataClean_enc);
avgMaint = ft_timelockanalysis(cfg,dataClean_maint);

%% LCMV

cfg = [];
cfg.method = 'lcmv';
cfg.grid = sourcemodel_and_leadfield;%grid;%
cfg.headmodel = hdm;
cfg.elec = elec_aligned;
cfg.lcmv.keepfilter = 'yes';
cfg.channel = dataScalp.label;
sourceavgFix = ft_sourceanalysis(cfg,avgFix);
sourceavgEnc = ft_sourceanalysis(cfg,avgEnc);
sourceavgMaint = ft_sourceanalysis(cfg,avgMaint);

%% contrast
cfg = [];
cfg.parameter = 'avg.pow';
cfg.operation = '((x1-x2)./x2)*100';
sourceEnc_baselined = ft_math(cfg,sourceavgEnc,sourceavgFix);
sourceMaint_baselined = ft_math(cfg,sourceavgMaint,sourceavgFix);
%% plot

templatefile = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\template\anatomy\single_subj_T1.nii';
template_mri = ft_read_mri(templatefile);
sourceEnc_baselined.pos=template_grid.pos;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, sourceEnc_baselined, template_mri);

cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
cfg.location = [64 -32 8];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,source_int);

%% plot in parceled brain space
templatefile = 'F:\Vasileios\Toolboxes\fieldtrip-20200315\external\spm8\templates\T1.nii';
template_mri = ft_read_mri(templatefile);
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
source_int  = ft_sourceinterpolate(cfg, sourceEnc_baselined, template_mri);
cfg=[];
parcel = ft_sourceparcellate(cfg, source_int, atlas);

dummy=atlas;
for i=1:length(parcel.pow)
    dummy.tissue(find(dummy.tissue==i))=parcel.pow(i);
end

source_int.parcel=dummy.tissue;
source_int.coordsys = 'mni';
cfg=[];
cfg.method = 'ortho';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'jet';
cfg.renderer = 'zbuffer';
cfg.location = [-42 -20 6];
cfg.atlas = atlas;
cfg.funcolorlim = [-30 30];
ft_sourceplot(cfg,source_int);


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-40 40];
cfg.funcolormap    = 'jet';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% cfg.projthresh     = 0.8;
cfg.camlight       = 'no';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

%% plot roi activity
temp_idx = strmatch('Frontal_Inf_Oper_L', atlas.tissuelabel);
temp_idx2 = strmatch('Cingulum_Ant_L', atlas.tissuelabel);
temp_idx = sort([temp_idx; temp_idx2]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx);
mask = ft_volumelookup(cfg, source_int);


cfg = [];
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-40 40];
cfg.funcolormap    = 'jet';
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.atlas = atlas;
cfg.camlight       = 'no';
cfg.roi = atlas.tissuelabel(temp_idx);
cfg.inputcoord = 'mni';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull

%% ----------------------
%% Statistics
cfg = [];
cfg.covariance = 'yes';
cfg.keeptrials = 'yes';
avg = ft_timelockanalysis(cfg,dataScalp);
avgFix_st = ft_timelockanalysis(cfg,dataScalp_Task_SS{1}{4});
avgEnc_st = ft_timelockanalysis(cfg,dataScalp_Task_SS{2}{4});
avgMaint_st = ft_timelockanalysis(cfg,dataScalp_Task_SS{3}{4});
cfg = [];
cfg.method = 'lcmv';
cfg.grid = grid;
cfg.headmodel = hdm;
cfg.channel = dataScalp.label;
cfg.elec = elec_aligned;
cfg.lcmv.keepfilter = 'yes';
sourceavg = ft_sourceanalysis(cfg,avg)

cfg = [];
cfg.method = 'lcmv';
cfg.grid = sourcemodel_and_leadfield;%grid;%
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.headmodel = hdm;
cfg.elec = elec_aligned;
% cfg.lcmv.keepfilter = 'yes';
% cfg.channel = dataScalp.label;
cfg.rawtrial = 'yes';
sourceavgFix_st = ft_sourceanalysis(cfg,avgFix_st);
sourceavgEnc_st = ft_sourceanalysis(cfg,avgEnc_st);
sourceavgMaint_st = ft_sourceanalysis(cfg,avgMaint_st);


% %% contrast
% cfg = [];
% cfg.parameter = 'avg.pow';
% cfg.operation = '((x1-x2)./x2)*100';
% sourceEnc_baselined_st = ft_math(cfg,sourceavgEnc_st,sourceavgFix_st);
% sourceMaint_baselined_st = ft_math(cfg,sourceavgMaint_st,sourceavgFix_st);
%% permutation test
cfg = [];
cfg.parameter    = 'pow';
cfg.dim          = grid.dim;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.atlas = atlas
cfg.numrandomization = 1000;

ntrials = numel(sourceavgFix_st.trial);
design  = zeros(2,2*ntrials);
design(1,1:ntrials) = 1;
design(1,ntrials+1:2*ntrials) = 2;
design(2,1:ntrials) = [1:ntrials];
design(2,ntrials+1:2*ntrials) = [1:ntrials];

cfg.design   = design;
cfg.ivar     = 1;
cfg.uvar     = 2;
stat = ft_sourcestatistics(cfg,sourceavgMaint_st,sourceavgFix_st);
stat.pos=template_grid.pos; % keep positions for plotting later
%% interpolate the result and the binary mask containing information of significant deferences per voxel.
% stat.inside = template_grid.inside;
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'stat';
cfg.interpmethod = 'nearest';
statint  = ft_sourceinterpolate(cfg, stat, template_mri);
% cfg.parameter    = 'mask';
maskint  = ft_sourceinterpolate(cfg, stat, template_mri);
% statint.mask = maskint.mask;


%% plot
statint.coordsys = 'mni';
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'stat';
cfg.maskparameter = 'mask';
cfg.atlas         = atlas;
cfg.location = 'max';

cfg.funcolorlim   = [0 20];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,statint);

%% plot only the significant
probplot=stat;
probplot.prob1=1-probplot.prob;
lowlim=0.95;
probplot.mask=(probplot.prob1>=lowlim);
probplot.anatomy=template_mri.anatomy;
probplot.dim = template_mri.dim;
probplot.coordsys = 'mni'
cfg               = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'prob1';
cfg.maskparameter = 'mask';
cfg.atlas         = atlas;
cfg.location = 'max';

cfg.funcolorlim   = [0.95 1];
cfg.funcolormap = 'jet';
ft_sourceplot(cfg,probplot);

%% create dummy struct
clear parcel;
cfg=[];
parcel = ft_sourceparcellate(cfg, statint, atlas);
parcelmask = ft_sourceparcellate(cfg, maskint, atlas);
dummy=atlas;
dummymask = atlas;
for i=1:length(parcel.stat)
    dummy.tissue(find(dummy.tissue==i))=parcel.stat(i);
    %       dummymask.tissue(find(dummymask.tissue==i))=parcelmask.mask(i);
end
%% plot the result
statint.parcel=dummy.tissue;
statint.coordsys = 'mni';
statint.mask  = dummymask.tissue;
cfg=[];
cfg.method = 'slice';
cfg.funparameter = 'parcel';
cfg.funcolormap    = 'bluewhitered_pos(128)';
cfg.maskparameter = 'mask';
cfg.renderer = 'zbuffer';
cfg.funcolorlim   = [0 10];
cfg.atlas = atlas;
ft_sourceplot(cfg,statint);

%% ----------------------
%% Connectivity analysis
%First, we interpolate the statistical result and the atlas.
strParcels = {'Frontal','central','Temporal','Occipital'};
clear gdata Granger_random_trial_sampling;
visualFlag = 0;
for iG = 1:length(strParcels)
    sprintf('\n\n\nCalculating Granger between hippocampus and %s regions\n\n\n', strParcels{iG})
    clear indxtempSupL;
    cfg = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter = 'tissue';
    stat_atlas = ft_sourceinterpolate(cfg, atlas, stat);
    
    %Determine the index of the label of interest,
    x = find(contains(atlas.tissuelabel,strParcels{iG}));%find(ismember(atlas.tissuelabel,'Frontal_Inf_Tri_L'));
    y = find(contains(atlas.tissuelabel,'_L'));% left hemisphere only
    indxtempSupL = intersect(x,y);
    %     indxtempSupL = find(stat_atlas.tissue==x);
    %     end
    
    %Next, we normalise the individual MRI to derive parameters allowing to
    %convert the mni- coordinates of the desired parcels into individual
    %coordinates. For this we use ft_warp_apply.
    
    template_grid=ft_convert_units(template_grid,'mm'); % ensure no unit mismatch
    norm=ft_volumenormalise([],mri);
    
    posTSL=template_grid.pos(indxtempSupL,:); % xyz positions in mni coordinates
    
    
    posback=ft_warp_apply(norm.params,posTSL,'sn2individual');
    btiposTSL= ft_warp_apply(pinv(norm.initial),posback); % xyz positions in individual coordinates
    
    %Now we create a source model for these particular locations only.
    
    cfg=[];
    cfg.headmodel = hdm;
    cfg.channel = dataScalp.label;
    cfg.sourcemodel.pos = [btiposTSL]./1000; % units of m
    cfg.elec = elec_aligned;
    sourcemodel_virt = ft_prepare_leadfield(cfg);
    
    %And repeat the source analysis steps for above but now for 1 parcel
    %represented in a total of 2 locations.
    
    % keep covariance in the output
    cfg = [];
    cfg.channel=dataScalp.label;
    cfg.covariance='yes';
    cfg.covariancewindow=[0 1];
    avg = ft_timelockanalysis(cfg,dataScalp);
    
    % perform source analysis on parcel of interest
    cfg=[];
    cfg.method='lcmv';
    cfg.grid = sourcemodel_virt;
    cfg.headmodel=hdm;
    cfg.lcmv.keepfilter='yes';
    cfg.lcmv.fixedori='yes';
    cfg.lcmv.lamda='5%';
    cfg.elec = elec_aligned
    source=ft_sourceanalysis(cfg, avg);
    
    cfg=[];
    cfg.method='lcmv';
    cfg.grid = sourcemodel_virt;
    cfg.headmodel=hdm;
    cfg.lcmv.keepfilter='yes';
    cfg.lcmv.fixedori='yes';
    cfg.lcmv.lamda='5%';
    cfg.elec = elec_aligned
    cfg.sourcemodel.filter = sourceavg.avg.filter;
    sourceavgFix_st_roi = ft_sourceanalysis(cfg,avgFix_st);
    sourceavgEnc_st_roi = ft_sourceanalysis(cfg,avgEnc_st);
    sourceavgMaint_st_roi = ft_sourceanalysis(cfg,avgMaint_st);
    
    % virtual channels
    if pID == 42
        cfg = [];
        cfg.channel = {'all','-T1','-T2'};
        dataScalp_Task_SS{2}{4} = ft_selectdata(cfg,dataScalp_Task_SS{2}{4});
        dataScalp_Task_SS{3}{4} = ft_selectdata(cfg,dataScalp_Task_SS{3}{4});
    elseif pID == 38
        cfg = [];
        cfg.channel = {'all','-Subm1','-Subm2'};
        dataScalp_Task_SS{2}{4} = ft_selectdata(cfg,dataScalp_Task_SS{2}{4});
        dataScalp_Task_SS{3}{4} = ft_selectdata(cfg,dataScalp_Task_SS{3}{4});
    elseif pID == 44
        cfg = [];
        cfg.channel = {'all','-Submm','-Submp','-T1','-T2','-EKGm','-EKGp'};
        dataScalp_Task_SS{2}{4} = ft_selectdata(cfg,dataScalp_Task_SS{2}{4});
        dataScalp_Task_SS{3}{4} = ft_selectdata(cfg,dataScalp_Task_SS{3}{4});
    elseif pID == 45
        cfg = [];
        cfg.channel = {'all','-Submm','-Submp','-EKGm','-EKGp'};
        dataScalp_Task_SS{2}{4} = ft_selectdata(cfg,dataScalp_Task_SS{2}{4});
        dataScalp_Task_SS{3}{4} = ft_selectdata(cfg,dataScalp_Task_SS{3}{4});
    end
    
    spatialfilter=cat(1,sourceavgEnc_st_roi.avg.filter{:});
    spatialfilter2=cat(1,sourceavgMaint_st_roi.avg.filter{:});
    % enc
    virtsens=[];
    for i=1:length(dataScalp_Task_SS{2}{4}.trial)
        virtsens.trial{i}=spatialfilter*100*dataScalp_Task_SS{2}{4}.trial{i};%spatialfilter*1000*dataScalp_Task_SS{2}{4}.trial{i};
        
    end
    virtsens.time=dataScalp_Task_SS{2}{4}.time;
    virtsens.fsample=dataScalp_Task_SS{2}{4}.fsample;
    indx=[indxtempSupL];
    for i=1:length(virtsens.trial{1}(:,1))
        virtsens.label{i}=[num2str(i)];
    end
    
    cfg = [];
    cfg.channel = virtsens.label(:);
    % cfg.channel = virtsens.label(1);
    
    cfg.avgoverchan = 'yes';
    virtsensTSL = ft_selectdata(cfg,virtsens);
    virtsensTSL.label = {'TSL'};
    
    virtsensparcel=ft_appenddata([],virtsensTSL,dataScalp_Task_SS{2}{4});
    cfg = [];
    cfg.viewmode = 'vertical';
    % ft_databrowser(cfg,virtsensparcel)
    
    % maint
    virtsens2=[];
    for i=1:length(dataScalp_Task_SS{3}{4}.trial)
        virtsens2.trial{i}=spatialfilter2*100*dataScalp_Task_SS{3}{4}.trial{i};%spatialfilter2*500*dataScalp_Task_SS{3}{4}.trial{i};
        
        
    end
    virtsens2.time=dataScalp_Task_SS{3}{4}.time;
    virtsens2.fsample=dataScalp_Task_SS{3}{4}.fsample;
    indx=[indxtempSupL];
    for i=1:length(virtsens2.trial{1}(:,1))
        virtsens2.label{i}=[num2str(i)];
    end
    
    
    cfg = [];
    cfg.channel = virtsens2.label(:);
    % cfg.channel = virtsens2.label(1);
    
    cfg.avgoverchan = 'yes';
    virtsens2TSL = ft_selectdata(cfg,virtsens2);
    virtsens2TSL.label = {'TSL'};
    
    
    virtsensparcel2=ft_appenddata([],virtsens2TSL,dataScalp_Task_SS{3}{4});
    cfg = [];
    cfg.viewmode = 'vertical';
    % ft_databrowser(cfg,virtsensparcel2)
    
    
    % Optional check
    % cfg=[];
    % tlkvc=ft_timelockanalysis(cfg, virtsensparcel2);
    % figure;
    % for i=1
    %     cfg=[];
    %     cfg.channel = tlkvc.label{i};
    %     cfg.parameter = 'avg';
    % %     cfg.xlim    = [-.1 1];
    %
    %     subplot(2,2,i);ft_singleplotER(cfg,tlkvc);
    % end
    
    %--------------
    %% Granger
    %% Parameters for loading the hipp data
    analysis_type = 'correct trials'; %'incorrect trials'; %or 'correct trials'
    % Granger Causality calculation between a hippocampal contact of a depth
    %electrode and a contact in the cortex of a depth electrode (outer contacts)
    flag_Hipp_GC_Ctx_Depth = 1;
    %% load hippocampal data
    if pID>=10
        [dataBip,TrialInformationTable_iEEG_clean] = loadHippocampalData_grid_patients(pID,analysis_type)
    else
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
    end
    %% Downsample to 40 Hz
    cfg = [];
    cfg.resamplefs = 40;%200;
    dataBipolarResampled = ft_resampledata(cfg,dataBip);
    combineddata = ft_resampledata(cfg,virtsensparcel);
    combineddata2 = ft_resampledata(cfg,virtsensparcel2);
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
    if pID <10
        TrialInformationTableCorrect.SetSize = TrialInformationTableCorrect.Setsize;
    end
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
    freq.freqcfg.taper     = 'dpss';
    freq.freqcfg.pad       = 20;
    freq.freqcfg.tapsmofrq = 1;
    Maintenance_freq       = ft_freqanalysis(freq.freqcfg, dataMaint);
    Encoding_freq          = ft_freqanalysis(freq.freqcfg, dataEnc);
    
    grangercfg = [];
    grangercfg.method  = 'granger';
    grangercfg.granger.conditional = 'no';
    grangercfg.granger.sfmethod = 'bivariate';
    % grngChan1 = 'PHL1-PHL3';
    % grngChan2 = 'TSL';
    grangercfg.channelcmb = {grngChan1,grngChan2};
    gdata{iG}.Maint{iPair,iSS}    = ft_connectivityanalysis(grangercfg, Maintenance_freq);
    gdata{iG}.Enc{iPair,iSS}       = ft_connectivityanalysis(grangercfg, Encoding_freq);
    
    
    %% visualization
    if visualFlag
        figure;
        
        light_blue = [0.30,0.75,0.93];
        light_red       = [1 0.45 0.45];
        Colors     = {'b',light_blue,'r',light_red};
        
        hold on;
        freq = Maintenance_freq.freq;
        semilogx(freq,gdata.Enc{iPair,iSS}.grangerspctrm(1,:),'LineWidth',3,'color',Colors{1});
        semilogx(freq,gdata.Enc{iPair,iSS}.grangerspctrm(2,:),'LineWidth',3,'color',Colors{2})
        semilogx(freq,gdata.Maint{iPair,iSS}.grangerspctrm(1,:),'LineWidth',3,'color',Colors{3})
        semilogx(freq,gdata.Maint{iPair,iSS}.grangerspctrm(2,:),'LineWidth',3,'color',Colors{4})
        xlim([4 20])
    end
    % cfg =[];
    % cfg.parameter = 'grangerspctrm';
    % cfg.xlim = [4 20];
    % figure
    % ft_connectivityplot(cfg,gdata.Enc{iPair,iSS},gdata.Maint{iPair,iSS})
    
    %% trials balance Granger convergence
    channelPair = {grngChan1,grngChan2}
    % nTrials_toSelect = length(find(TrialInformationTable.Scalp.Correct == 0));
    nTrials_toSelect = 20;
    
    [Granger_random_trial_sampling{iG}] = calculateGranger_random_trial_sampling(pID,iSS,nTrials_toSelect,dataMaint,dataEnc,channelPair)
    ids = [42 38 37 40 45 44; 10 11 12 13 14 15]';
    if pID<10
        id = pID;
    else
        idx = find(ids(:,1) == pID);
        id = ids(idx,2);
    end
    strPath = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\multiParcel_Granger_analysis\';
    save([strPath 'Granger_random_trial_sampling_P' ,num2str(id)],'Granger_random_trial_sampling','-v7.3');
end
%% visualization
figure;
id = 6;
light_blue = [0.30,0.75,0.93];
light_red       = [1 0.45 0.45];
Colors     = {'b',light_blue,'r',light_red};

hold on;
freq = Maintenance_freq.freq;
semilogx(freq,Granger_random_trial_sampling{id}.Median_gdata.Enc{2},'LineWidth',3,'color',Colors{1});
semilogx(freq,Granger_random_trial_sampling{id}.Median_gdata.Enc{1},'LineWidth',3,'color',Colors{2})
semilogx(freq,Granger_random_trial_sampling{id}.Median_gdata.Maint{2},'LineWidth',3,'color',Colors{3})
semilogx(freq,Granger_random_trial_sampling{id}.Median_gdata.Maint{1},'LineWidth',3,'color',Colors{4})
xlim([4 20])
%% convergence plot calculation
nTrials_toSelect = 40;

for i = [1:19 20:5:nTrials_toSelect]
    [Granger_random_trial_sampling_convergence{i}] = calculateGranger_random_trial_sampling(1,iSS,i,dataMaint,dataEnc,channelPair);
end


%% %% convergence plot visualization
freq_start = 4;
freq_end = 8;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);

ind = [1 10:15 2:9];
k=0;
for i = [1:19 20:5:nTrials_toSelect]
    DGranger_enc(i-k) = min(Granger_random_trial_sampling_convergence{i}{1}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - Granger_random_trial_sampling_convergence{i}{1}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
    DGranger_maint(i-k) = max(Granger_random_trial_sampling_convergence{i}{1}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - Granger_random_trial_sampling_convergence{i}{1}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));
    for j = 1:size(Granger_random_trial_sampling_convergence{i}{1}.gdata.Maint{1},1)
        dgranger_enc_conv_per_permutation(j) = min(Granger_random_trial_sampling_convergence{i}{1}.gdata.Enc{1}(j,freq_wnd_a:freq_wnd_b) - Granger_random_trial_sampling_convergence{i}{1}.gdata.Enc{2}(j,freq_wnd_a:freq_wnd_b));
        dgranger_maint_conv_per_permutation(j) = max(Granger_random_trial_sampling_convergence{i}{1}.gdata.Maint{2}(j,freq_wnd_a:freq_wnd_b) - Granger_random_trial_sampling_convergence{i}{1}.gdata.Maint{1}(j,freq_wnd_a:freq_wnd_b));
    end
    error_dgranger_enc(i-k) = std(dgranger_enc_conv_per_permutation);
    error_dgranger_maint(i-k) = std(dgranger_maint_conv_per_permutation);
    %     k = k+1;
    
end



x = [1:19 20:5:40];
y_enc = DGranger_enc(x);
y_maint = DGranger_maint(x);
error_dgranger_enc = error_dgranger_enc(x);
error_dgranger_maint = error_dgranger_maint(x);
figure;
subplot(2,1,2)
hold on;
errorbar(x,y_enc,error_dgranger_enc,'color','b');
a = fill([x(2:end) fliplr(x(2:end))],[y_enc(2:end)-error_dgranger_enc(2:end) fliplr(y_enc(2:end)+error_dgranger_enc(2:end))], 'b');
set(a,'FaceAlpha',0.3,'FaceColor','b','linestyle','-.','EdgeColor','b')
xlim([1 40])


xlabel('Number of trials')
ylabel('\DeltaGranger (%)')
set(gca,'Fontsize',14,'XTick',[0:5:40],'XTickLabel',[0:5:40],'TickDir','out')
subplot(2,1,1)
hold on;
errorbar(x,y_maint,error_dgranger_maint,'color','r')
a = fill([x(2:end) fliplr(x(2:end))],[y_maint(2:end)-error_dgranger_maint(2:end) fliplr(y_maint(2:end)+error_dgranger_maint(2:end))], 'r');
set(a,'FaceAlpha',0.3,'FaceColor','r','linestyle','-.','EdgeColor','r')
xlim([1 40])
set(gca,'Fontsize',14,'XTick',[0:5:40],'XTickLabel',[0:5:40],'TickDir','out')
xlabel('Number of trials')
ylabel('\DeltaGranger (%)')


%% Time frequency Granger
channel_cmb = {'AL1-AL3', 'TSL'};
cfg = [];
cfg.channel = channel_cmb;
data_ENC_single_chans = ft_selectdata(cfg,dataEnc);
data_MAINT_single_chans = ft_selectdata(cfg,dataMaint);

iSS = 4;
cfg           = [];
cfg.method    = 'mtmconvol';
cfg.output    = 'fourier';
cfg.taper     = 'hanning';
cfg.foi       = 0:20;
cfg.t_ftimwin = 0.2.*ones(size(4./cfg.foi',1), 1);
cfg.tapsmofrq = 0.05*cfg.foi;
cfg.toi       = -5:0.22:-3;
cfg.pad = 20;
cfg.keeptrials = 'yes';
CrossFreq_enc     = ft_freqanalysis(cfg, data_ENC_single_chans);

cfg           = [];
cfg.method    = 'mtmconvol';
cfg.output    = 'fourier';
cfg.taper     = 'hanning';
cfg.foi       = 0:20;
cfg.t_ftimwin = 0.2.*ones(size(4./cfg.foi',1), 1);
cfg.tapsmofrq = 0.05*cfg.foi;
% cfg.toi       = -2:0.21:0;
cfg.toi       = -2:0.22:0;

cfg.pad = 20;
cfg.keeptrials = 'yes';
CrossFreq_maint     = ft_freqanalysis(cfg, data_MAINT_single_chans);

grangercfg = [];
grangercfg.method  = 'granger';
grangercfg.channelcmb = channel_cmb;
Granger_Rand_SS_enc{iSS}= ft_connectivityanalysis(grangercfg,CrossFreq_enc);
Granger_Rand_SS_maint{iSS}= ft_connectivityanalysis(grangercfg,CrossFreq_maint);


CortexHipp_Enc = squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(2,1,:,:));
HippCortex_Enc = (squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(1,2,:,:)));
Granger_SS_enc{iSS} = (HippCortex_Enc-CortexHipp_Enc)*100;


%
grangerTimeAxis = Granger_Rand_SS_enc{iSS}.time;
grangerFreqAxis = Granger_Rand_SS_enc{iSS}.freq;
figure;
clim = [-0.1 0.1]*100%[-0.1 0.15];%[-0.05 0.05];
contourf(grangerTimeAxis,grangerFreqAxis,Granger_SS_enc{iSS},100,'LineColor','none');

set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[4,5 10 20] ,'YTickLabel',[4,5,10,20]) %background color of grid

colormap(bluewhitered(128));
set(gca,'XTick',[-5 -3], 'XTickLabel',[-5 -3])
% ylim([4 100]);
ylim([4 20])
colorbar;


% maintenance
CortexHipp_Maint = squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(2,1,:,:));
HippCortex_Maint = (squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(1,2,:,:)));
Granger_SS_maint{iSS} = (HippCortex_Maint-CortexHipp_Maint)*100

grangerTimeAxis = Granger_Rand_SS_maint{iSS}.time;
grangerFreqAxis = Granger_Rand_SS_maint{iSS}.freq;
figure;
clim = [-0.1 0.1]*100%[-0.1 0.15];%[-0.05 0.05];
contourf(grangerTimeAxis,grangerFreqAxis,Granger_SS_maint{iSS},100,'LineColor','none');

set(gca,'clim',clim,'yscale','log');
set(gca,'ytick',[4,5 10 20] ,'YTickLabel',[4,5,10,20]) %background color of grid

colormap(bluewhitered(128));
set(gca,'XTick',[-2 0], 'XTickLabel',[-2 0])
% ylim([4 100]);
ylim([4 20])
colorbar;


%% Time frequency statistics

%% prepare the data

%enc
Granger_Rand_SS_enc{iSS}.powspctrm{1} = squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(1,2,:,:))
Granger_Rand_SS_enc{iSS}.powspctrm{2} = squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(2,1,:,:))
nFr = length(grangerFreqAxis);
nTim = length(grangerTimeAxis);
Granger_Rand_SS_enc{iSS}.powspctrm = reshape(cell2mat(Granger_Rand_SS_enc{iSS}.powspctrm),2,nFr,nTim);
Granger_Rand_SS_enc{iSS}.dimord = 'chan_freq_time';

Granger_Rand_SS_enc{iSS}.trial{1} = squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(1,2,:,:));
Granger_Rand_SS_enc{iSS}.trial{2} = squeeze(Granger_Rand_SS_enc{iSS}.grangerspctrm(2,1,:,:));

%maint
Granger_Rand_SS_maint{iSS}.powspctrm{1} = squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(1,2,:,:))
Granger_Rand_SS_maint{iSS}.powspctrm{2} = squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(2,1,:,:))
nFr = length(grangerFreqAxis);
nTim = length(grangerTimeAxis);
Granger_Rand_SS_maint{iSS}.powspctrm = reshape(cell2mat(Granger_Rand_SS_maint{iSS}.powspctrm),2,nFr,nTim);
Granger_Rand_SS_maint{iSS}.dimord = 'chan_freq_time';

Granger_Rand_SS_maint{iSS}.trial{1} = squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(1,2,:,:));
Granger_Rand_SS_maint{iSS}.trial{2} = squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(2,1,:,:));

Gr_maint_stat_chan_1 = Granger_Rand_SS_maint{iSS};
Gr_maint_stat_chan_1.label = Gr_maint_stat_chan_1.label{1}
Gr_maint_stat_chan_1.grangerspctrm =  squeeze(Gr_maint_stat_chan_1.grangerspctrm(1,2,:,:));
Gr_maint_stat_chan_1.powspctrm =  squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(1,2,:,:));

Gr_maint_stat_chan_2 = Granger_Rand_SS_maint{iSS};
Gr_maint_stat_chan_2.label = Gr_maint_stat_chan_2.label{2};
Gr_maint_stat_chan_2.grangerspctrm =  squeeze(Gr_maint_stat_chan_2.grangerspctrm(2,1,:,:));
Gr_maint_stat_chan_2.powspctrm =  squeeze(Granger_Rand_SS_maint{iSS}.grangerspctrm(2,1,:,:));

%%
cfg                  = [];
cfg.method           = 'montecarlo'; % use montecarlo to permute the data
cfg.statistic        = 'indepsamplesT'; % function to use when ...
% calculating the ...
% parametric t-values
cfg.alpha            = 0.05; % corresponds to an alpha level of 0.05, since ...
% two tests are made ...
% (negative and positive: 2*0.025=0.05)
cfg.frequency = [0 20];
cfg.latency   = [-2 0];
cfg.parameter = 'powspctrm';
cfg.correctm         = 'fdr'; % the correction to use
% cfg.clusteralpha     = 0.05; % the alpha level used to determine whether or ...
%                              % not a channel/time pair can be included in a ...
%                              % cluster
% cfg.clustertail      = 1; % two-way t-test
cfg.clusterstatistic = 'max';
cfg.numrandomization = 1000;  % number of permutations run
%
cfg.design    = zeros(1, 2)';
cfg.design(1)  = 1; % indicating which trials belong ...
cfg.design(2) = 2; % to what category
cfg.ivar             = 1; % indicating that the independent variable is found in ...
% first row of cfg.design
cfg.minnbchan        = 2; % minimum number of channels required to form a cluster

stat_t_cluster_freq = ft_freqstatistics(cfg,Gr_maint_stat_chan_1, Gr_maint_stat_chan_2);



%%
cfg=[];
cfg.channel = {'PHL2-PHL3'};
cfg.parameter = 'stat';
cfg.maskparameter = 'mask';
cfg.maskstyle     = 'outline';
% cfg.zlim = [-5 5]
% cfg.xlim    = [0 .7];
figure;
subplot(2,2,1);ft_singleplotTFR(cfg,stat_t_cluster_freq);
