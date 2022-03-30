strPath_random_shuflling = 'F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\random sampling Granger convergence\multiParcel_Granger_analysis\';
cd(strPath_random_shuflling);
files = dir('*.mat');
nSubjs = length(files);

freq_start = 4;%3;
freq_end = 8;%5;
x_enc = 1;
x_maint = 2;
freq_wnd_a = find(freq == freq_start);
freq_wnd_b = find(freq == freq_end);
strParcels = {'Frontal','Central','Temporal','Occipital'};

ind = [1 10:15 2:9];
for j = 1:length(strParcels)
    
    for  i= 1:length(files)
        dataDGranger{i} = load(files(i).name);
        dataDGranger_subj_parc{i,j} = dataDGranger{i}.Granger_random_trial_sampling{j}{ind(i)};
        DGranger_enc_prc(i,j) = min(dataDGranger_subj_parc{i,j}.Median_gdata.Enc{1}(freq_wnd_a:freq_wnd_b) - dataDGranger_subj_parc{i,j}.Median_gdata.Enc{2}(freq_wnd_a:freq_wnd_b));
        DGranger_maint_prc(i,j) = max(dataDGranger_subj_parc{i,j}.Median_gdata.Maint{2}(freq_wnd_a:freq_wnd_b) - dataDGranger_subj_parc{i,j}.Median_gdata.Maint{1}(freq_wnd_a:freq_wnd_b));
    end
end
median_DGranger_enc_parc = median(DGranger_enc_prc);
median_DGranger_maint_parc = median(DGranger_maint_prc);




%%
load('F:\Vasileios\Task Analysis\Analysis Results\Granger_source_level\variables\groupStats_vars\granger_surface.mat')
atlas = granger_surface.atlas;
mask = granger_surface.mask;
source_int = granger_surface.source_int;
source_int.inside(:) = 1;
source_int.parcel(:) = 0;
source_int_maint = source_int;
temp_idx = strmatch('Frontal_Inf_Oper_L', atlas.tissuelabel);
temp_idx3 = strmatch('Frontal_Inf_Tri_L', atlas.tissuelabel);
temp_idx4 = strmatch('Frontal_Mid_L', atlas.tissuelabel);

temp_idx = sort([temp_idx; temp_idx4;temp_idx3]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx);
mask1 = ft_volumelookup(cfg, source_int);
source_int.parcel(mask1) = median_DGranger_enc_LPFC;

temp_idx_AC = strmatch('Temporal_Sup_L', atlas.tissuelabel);
temp_idx_AC = sort([temp_idx_AC]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_AC);
mask2 = ft_volumelookup(cfg, source_int);

temp_idx_PPC = strmatch('Parietal_Sup_L', atlas.tissuelabel);
temp_idx_PPC = sort([temp_idx_PPC]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_PPC);
mask3 = ft_volumelookup(cfg, source_int);


temp_idx_V1 = strmatch('Occipital_Mid_L', atlas.tissuelabel);
temp_idx_V1 = sort([temp_idx_V1]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_V1);
mask4 = ft_volumelookup(cfg, source_int);


temp_idx_Broca = strmatch('Frontal_Inf_Tri_L', atlas.tissuelabel);
temp_idx_Broca = sort([temp_idx_Broca]);
cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_Broca);
% mask5 = ft_volumelookup(cfg, source_int);

source_int.parcel(mask1) = median_DGranger_enc_LPFC;
source_int.parcel(mask2) = median_DGranger_enc_AC;
source_int.parcel(mask3) = median_DGranger_enc_PPC;
source_int.parcel(mask4) = median_DGranger_enc_V1;
% source_int.parcel(mask5) = median_DGranger_enc_Broca;
%%
strTablePath = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\atlas_parcels_Granger_projection.xlsx';
tblParcels = readtable(strTablePath)
groups  = unique(tblParcels.Group);
groups = groups(2:end);
indx_group =[];
for i =groups(1):groups(end)
    ind = find(tblParcels.Group==i)
    
    group_parcels{i} = tblParcels.AtlasParcels(ind);
    for j = 1:length(group_parcels{i})
        indx_group{i,j} =    strmatch(group_parcels{i}{j}, atlas.tissuelabel);
       
    end
     temp_index_group{i} = [indx_group{i,:}];
     cfg = [];
     cfg.inputcoord = 'mni';
     cfg.atlas = atlas;
     cfg.roi = atlas.tissuelabel(temp_index_group{i});
     mask_grp{i} = ft_volumelookup(cfg, source_int);
     if i<=4
     source_int.parcel(mask_grp{i}) = median_DGranger_enc_parc(i)+1.5;

     else
         source_int.parcel(mask_grp{i}) = median_DGranger_enc_PPC;
     end

end

figure;
cfg = [];
cfg.figure = gcf;
% cfg.axis = axes(ha(5));
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [-5 0];
% cfg.colorbar = 'no'
cfg.funcolormap    = colormap(bluewhitered_neg(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.atlas = atlas;
cfg.camlight       = 'no';
cfg.roi = atlas.tissuelabel(sort([temp_idx_AC; temp_idx;temp_idx_PPC;temp_idx_V1;cell2mat(temp_index_group)']));
cfg.inputcoord = 'mni';
ft_sourceplot(cfg, source_int);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(gca,bluewhitered_neg(128))
set(gca,'FontSize',12)

%%


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx);
mask1 = ft_volumelookup(cfg, source_int_maint);
source_int_maint.parcel(mask1) = median_DGranger_maint_LPFC;


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_AC);
mask2 = ft_volumelookup(cfg, source_int_maint);

cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_PPC);
mask3 = ft_volumelookup(cfg, source_int_maint);


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_V1);
mask4 = ft_volumelookup(cfg, source_int_maint);


cfg = [];
cfg.inputcoord = 'mni';
cfg.atlas = atlas;
cfg.roi = atlas.tissuelabel(temp_idx_Broca);
% mask5 = ft_volumelookup(cfg, source_int_maint);


source_int_maint.parcel(mask1) = median_DGranger_maint_LPFC;
source_int_maint.parcel(mask2) = median_DGranger_maint_AC;
source_int_maint.parcel(mask3) = median_DGranger_maint_PPC;
source_int_maint.parcel(mask4) = median_DGranger_maint_V1;
% source_int_maint.parcel(mask5) = median_DGranger_maint_Broca;


strTablePath = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\atlas_parcels_Granger_projection.xlsx';
tblParcels = readtable(strTablePath)
groups  = unique(tblParcels.Group);
groups = groups(2:end);
indx_group =[];
for i =groups(1):groups(end)
    ind = find(tblParcels.Group==i)
    
    group_parcels{i} = tblParcels.AtlasParcels(ind);
    for j = 1:length(group_parcels{i})
        indx_group{i,j} =    strmatch(group_parcels{i}{j}, atlas.tissuelabel);
       
    end
     temp_index_group{i} = [indx_group{i,:}];
     cfg = [];
     cfg.inputcoord = 'mni';
     cfg.atlas = atlas;
     cfg.roi = atlas.tissuelabel(temp_index_group{i});
     mask_grp{i} = ft_volumelookup(cfg, source_int);
     if i<=4
     source_int_maint.parcel(mask_grp{i}) = median_DGranger_maint_parc(i)-2;

     else
         source_int.parcel(mask_grp{i}) = median_DGranger_enc_PPC;
     end

end
figure;cfg = [];
cfg.figure = gcf;
% cfg.axis = axes5;
% cfg.colorbar = 'no'
cfg.method         = 'surface';
cfg.funparameter   = 'parcel';
cfg.funcolorlim    = [0 5];
cfg.funcolormap    = colormap(gca,bluewhitered_pos(128));
cfg.projmethod     = 'nearest';
cfg.surfinflated   = 'surface_inflated_both_caret.mat';
cfg.atlas = atlas;
cfg.camlight       = 'no';
cfg.roi = atlas.tissuelabel(sort([temp_idx_AC; temp_idx;temp_idx_PPC;temp_idx_V1;cell2mat(temp_index_group)']));
cfg.inputcoord = 'mni';
ft_sourceplot(cfg, source_int_maint);
view ([-70 20 50])
light ('Position',[-70 20 50])
material dull
colormap(gca,bluewhitered_pos(128));
set(gca,'FontSize',12)
cb = get(gca,'colorbar');
