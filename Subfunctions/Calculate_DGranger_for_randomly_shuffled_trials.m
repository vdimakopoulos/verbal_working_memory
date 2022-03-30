load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger convergence\Granger_random_trial_sampling_scalp.mat')
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger convergence\Granger_random_trial_sampling_grid.mat')

nScalp = length(Granger_random_trial_sampling_scalp);
nGrid = length(Granger_random_trial_sampling);
maxFreqMaint =[];
maxFreqEnc =[];
valMaint_CH = [];
valEnc_HC = [];
for i = 1:nGrid
   maxValMaint = max(Granger_random_trial_sampling{i}.Median_gdata.Maint{1}(2:27)); 
   maxFreqMaint(i) = find(Granger_random_trial_sampling{i}.Median_gdata.Maint{1}(2:27)==maxValMaint);
   maxValEnc = max(Granger_random_trial_sampling{i}.Median_gdata.Enc{1}(2:27));
   maxFreqEnc(i) = find(Granger_random_trial_sampling{i}.Median_gdata.Enc{1}(2:27)==maxValEnc);
   valMaint_CH(i) = Granger_random_trial_sampling{i}.Median_gdata.Maint{2}(maxFreqMaint(i)+1);
   valEnc_HC(i) = Granger_random_trial_sampling{i}.Median_gdata.Enc{2}(maxFreqEnc(i)+1);
   DGranger_MaintGrid(i) = maxValMaint - valMaint_CH(i);
   DGranger_EncGrid(i) = valEnc_HC(i) - maxValEnc;

end

maxFreqMaint =[];
maxFreqEnc =[];
valMaint_CH = [];
valEnc_HC = [];
for i = 1:nScalp
   maxValMaint = max(Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{1}(2:27)); 
   maxFreqMaint(i) = find(Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{1}(2:27)==maxValMaint);
   maxValEnc = max(Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{1}(2:27));
   maxFreqEnc(i) = find(Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{1}(2:27)==maxValEnc);
   valMaint_CH(i) = Granger_random_trial_sampling_scalp{i}.Median_gdata.Maint{2}(maxFreqMaint(i)+1);
   valEnc_HC(i) = Granger_random_trial_sampling_scalp{i}.Median_gdata.Enc{2}(maxFreqEnc(i)+1);
   DGranger_Maint(i) = maxValMaint - valMaint_CH(i);
   DGranger_Enc(i) = valEnc_HC(i) - maxValEnc;

end

x_enc = 5;
x_maint =46;
nSubjs = nScalp
[EncDiff_sorted,Ienc] =  sort(DGranger_Enc,'descend');
[MaintDiff_sorted,Imaint] = sort(DGranger_Maint);
xValues = [repmat(x_enc,nSubjs,1)'; repmat(x_maint,nSubjs,1)']'%[1:15]+3*15]';
fig = figure;
for i=1:nSubjs
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_sorted(i); MaintDiff_sorted(j)]','k-','LineWidth',1)
    hold on
end

enc = scatter(repmat(x_enc,nSubjs,1)',EncDiff_sorted,20,'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint = scatter(repmat(x_maint,nSubjs,1)',MaintDiff_sorted,20,'Marker','^','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',1);
yline(0,'k-','LineWidth',1);


hold on;

[EncDiff_ECoG_sorted,Ienc] =  sort(DGranger_EncGrid,'descend');
[MaintDiff_ECoG_sorted,Imaint] = sort(DGranger_MaintGrid);

enc_grid = scatter(repmat(x_enc,length(EncDiff_ECoG_sorted),1)',EncDiff_ECoG_sorted,20,'Marker','o','MarkerEdgeColor','b','LineWidth',1);
hold on;
maint_grid = scatter(repmat(x_maint,length(MaintDiff_ECoG_sorted),1)',MaintDiff_ECoG_sorted,20,'Marker','o','MarkerEdgeColor','r','LineWidth',1);

for i=1:length(EncDiff_ECoG_sorted)
    j= find(Imaint == Ienc(i));
    plot(xValues(i,:),[EncDiff_ECoG_sorted(i); MaintDiff_ECoG_sorted(j)]','k-','LineWidth',1)
    hold on
end

%% Median for Maintenance
MaintDiff_sorted = sort([MaintDiff_sorted,MaintDiff_ECoG_sorted])
median_Maint = median(MaintDiff_sorted);
prcntile75 = prctile(MaintDiff_sorted,75);
prcntile25 = prctile(MaintDiff_sorted,25);
interquartile = iqr(MaintDiff_sorted);
Upper_adjacent = prcntile75 -1.5*interquartile;
Lower_adjacent = prcntile25 +1.5*interquartile;

line(repmat(x_maint+2,2,1)',[prcntile75 prcntile25],'Color','k')
plot([47:49],repmat(median_Maint,3,1),'r','LineWidth',1)
plot([47:49],repmat(median_Maint,3,1)-0.4,'r','LineWidth',1)
plot([47:49],repmat(prcntile75,3,1),'k');
plot([47:49],repmat(prcntile25,3,1),'k');

%% Median for Encoding
EncDiff_sorted = sort([EncDiff_sorted,EncDiff_ECoG_sorted]);
median_Enc   = median(EncDiff_sorted);
prcntile75 = prctile(EncDiff_sorted,75);
prcntile25 = prctile(EncDiff_sorted,25);
interquartile = iqr(EncDiff_sorted);
Upper_adjacent = prcntile75 -1.5*interquartile;
Lower_adjacent = prcntile25 +1.5*interquartile;

line(repmat(x_enc-2,2,1)',[prcntile75 prcntile25],'Color','k')
plot([2:4],repmat(median_Enc,3,1),'r','LineWidth',1)
plot([2:4],repmat(median_Enc,3,1)-0.4,'r','LineWidth',1)
plot([2:4],repmat(prcntile75,3,1),'k');
plot([2:4],repmat(prcntile25,3,1),'k');

%% Axes Properties
set(gca,'FontSize',16,'box','off');
ylabel('\DeltaGranger (%)');
ylim([-20 15])
set(gca,'XTick',[5 46],'XTickLabel',[{'Encoding','Maintenance'}],'TickDir','out');
