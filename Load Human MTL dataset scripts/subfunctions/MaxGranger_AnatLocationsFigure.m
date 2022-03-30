%% Max Granger values - Anterior Posterior Figure
% load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Final_SubjGranger_spectra.mat');
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\MaxGranger_AnatLocations.mat')
load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Final_SubjGranger_spectra_recalculatedGranger.mat');


for i =1:numel(Granger_spectra)
    maxEnc(i) = max(Granger_spectra{i}.Enc.grangerspctrm(1,:));
    maxMaint(i) = max(Granger_spectra{i}.Maint.grangerspctrm(1,:));
        
end

figure;
set(gcf,'Color','white');
for i =1:numel(Granger_spectra)
    if contains(AnatLocations_Granger.Labels(i),'L')
        markerType = 'o';
    else
        markerType = 'square';
    end
   scatter(AnatLocations_Granger.y(i), maxEnc(i)*100,100,'Marker',markerType,'LineWidth',2,'MarkerFaceColor','flat','MarkerEdgeColor','none' );%,'filled');
    hold on;
end
xlim([-36 2])
ylabel('Max Granger (%)');
xlabel('Posterior      y(mm)       Anterior');
set(gca,'XTick',[- 30 -20 -10 0],'XTickLabel',[- 30 -20 -10 0],'box','off','FontSize',20);
% c = polyfit(AnatLocations_Granger.y(:),maxEnc'*100,1);
% y_est = polyval(c,AnatLocations_Granger.y(:));
hold on;
% plot(AnatLocations_Granger.y(:),y_est,'k')
% [b, bint,r,rint,stats] = regress(maxEnc'*100,AnatLocations_Granger.y(:))
model =fitlm(AnatLocations_Granger.y(:),maxEnc'*100)
plot(model)
% 
% [maxEnc_sorted,ind] = sort(maxEnc);
% NoOutlier_AnatLocations = AnatLocations_Granger.y(ind(2:14))
% maxEnc_no_outliers = maxEnc_sorted(2:14);
% model =fitlm(NoOutlier_AnatLocations,maxEnc_no_outliers'*100)
% c_no_outliers = polyfit(NoOutlier_AnatLocations,maxEnc_no_outliers'*100,1);
% y_est_no_outliers = polyval(c_no_outliers,NoOutlier_AnatLocations);
% hold on;
% plot(NoOutlier_AnatLocations,y_est_no_outliers,'k')
% plot(AnatLocations_Granger.y,y_est,'k')


figure;
set(gcf,'Color','white');
for i =1:numel(Granger_spectra)
    if contains(AnatLocations_Granger.Labels(i),'L')
        markerType = 'o';
    else
        markerType = 'square';
    end
   scatter(AnatLocations_Granger.y(i), maxMaint(i)*100,100,'Marker',markerType,'LineWidth',2,'MarkerFaceColor','flat','MarkerEdgeColor','none' );
    hold on;
end
xlim([-36 2])
ylabel('Max Granger (%)')
xlabel('Posterior      y(mm)       Anterior')
set(gca,'XTick',[- 30 -20 -10 0],'XTickLabel',[- 30 -20 -10 0],'box','off','FontSize',20);
model =fitlm(AnatLocations_Granger.y(:),maxMaint'*100,'RobustOpts','on')
plot(model)
% c = polyfit(AnatLocations_Granger.y(:),maxMaint'*100,1);
% y_est = polyval(c,AnatLocations_Granger.y(:));
% hold on;
% plot(AnatLocations_Granger.y(:),y_est,'k')

% 
% 
% [maxMaint_sorted,ind] = sort(maxMaint);
% NoOutlier_AnatLocations = AnatLocations_Granger.y(ind(2:14))
% maxMaint_no_outliers = maxMaint_sorted(2:14);
% model =fitlm(NoOutlier_AnatLocations,maxMaint_no_outliers'*100)
% c_no_outliers = polyfit(NoOutlier_AnatLocations,maxMaint_no_outliers'*100,1);
% y_est_no_outliers = polyval(c_no_outliers,NoOutlier_AnatLocations);
% hold on;
% plot(NoOutlier_AnatLocations,y_est_no_outliers,'k')
