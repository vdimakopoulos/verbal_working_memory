clear elecNames elecColors
GridLines = {'A','B','C','D','E','F','G','H'}
for i =1:8
    for j = 1:8
        elecNames{(i-1)*8+j} = sprintf('%s%d',GridLines{i},j)

    end
end
Granger_topo = reshape(PLV_PSD_Granger_Topologies.Granger_enc,1,64);
ind = find(Granger_topo<(min(Granger_topo) - min(Granger_topo)*0.2))
Granger_topo(setdiff([1:length(Granger_topo)],ind)) = 0;
elecColors = Granger_topo

cfg = [];
cfg.elecColorScale =[-10 0]
cfg.view = 'l'
cfg.elecCoord='PIAL'
cfg.showLabels='y';
cfg.ignoreDepthElec = 'y'
cfg.elecColors = elecColors';
cfg.elecNames = elecNames';
cfg.elecCmapName ='bluewhitered_neg(512)'%'bluewhitered(128,[-10 0])'%'parula'%
cfg.elecShape = 'sphere'
cfg.fsurfSubDir = 'F:/Vasileios/Task Analysis/Data/Freesurfer Data/'
cfg.elecSize = 4;
% cfg. elecCbar = 'y';
cfg.title = [];
cfgOut = plotPialSurf('USZ',cfg)
