function [modelCorr modelIncorr] = Plot_PowerCorrelation(Power_MaintCorrect,Power_MaintIncorrect,band,channelPair)

nCorrectTrials = length(Power_MaintCorrect);
nIncorrectTrials = size(Power_MaintIncorrect,2);

% fig = figure;
for i = 1:nCorrectTrials
    
    LogPower_x(i) =  Power_MaintCorrect{1,i}.meanBaselinedPower;
    LogPower_y(i) = Power_MaintCorrect{2,i}.meanBaselinedPower;
    
end
for i = 1:nIncorrectTrials
    
    LogPower_x_Incorrect(i) = Power_MaintIncorrect{1,i}.meanBaselinedPower;
    LogPower_y_Incorrect(i) = Power_MaintIncorrect{2,i}.meanBaselinedPower;
    
end
%% Reject Outliers
% Correct trials
indOutlier_1 = find(LogPower_x>5);
indOutlier_2 = find(LogPower_x<-2);
indOutlier_total = union(indOutlier_1,indOutlier_2);
totalPoints =length(LogPower_x);
% ind_keepTrials = setdiff([1:totalPoints],indOutlier_total);

indOutlier_1_y = find(LogPower_y>5);
indOutlier_2_y = find(LogPower_y<-5);
indOutlier_total_y = union(indOutlier_1_y,indOutlier_2_y);
indOutlier_total_xy = union(indOutlier_total,indOutlier_total_y)

ind_keepTrials_xy = setdiff([1:totalPoints],indOutlier_total_xy);

LogPower_x = LogPower_x(ind_keepTrials_xy);
LogPower_y = LogPower_y(ind_keepTrials_xy);



indOutlier_1 = find(LogPower_x_Incorrect>5);
indOutlier_2 = find(LogPower_x_Incorrect<-2);
indOutlier_total = union(indOutlier_1,indOutlier_2);
totalPoints =length(LogPower_x_Incorrect);
% ind_keepTrials = setdiff([1:totalPoints],indOutlier_total);

indOutlier_1_y = find(LogPower_y_Incorrect>5);
indOutlier_2_y = find(LogPower_y_Incorrect<-5);
indOutlier_total_y = union(indOutlier_1_y,indOutlier_2_y);
indOutlier_total_xy = union(indOutlier_total,indOutlier_total_y)

ind_keepTrials_xy = setdiff([1:totalPoints],indOutlier_total_xy);

LogPower_x_Incorrect = LogPower_x_Incorrect(ind_keepTrials_xy);
LogPower_y_Incorrect = LogPower_y_Incorrect(ind_keepTrials_xy);




%% Plot the power covariation
modelCorr = fitlm(LogPower_x,LogPower_y,'RobustOpts','on');
modelIncorr = fitlm(LogPower_x_Incorrect,LogPower_y_Incorrect,'RobustOpts','on');
if modelCorr.Coefficients{2,4}<0.05 && modelIncorr.Coefficients{2,4}>0.05
    figure;
    scatter(LogPower_x,LogPower_y,30,'b','Marker','o','LineWidth',2,'MarkerFaceColor','flat','MarkerEdgeColor','none' );
    hold on;
    scatter(LogPower_x_Incorrect,LogPower_y_Incorrect,30,'r','Marker','o','LineWidth',2,'MarkerFaceColor','flat','MarkerEdgeColor','none' );
    
    
    
    h = plot(modelCorr);
    set(h,'Color','b')
    ht = plot(modelIncorr);
    hLegend = findobj(gcf, 'Type', 'Legend');
    set(hLegend,'Visible','off')
    ttl = title(sprintf('R_corr=%.03f, p_corr=%.03f,R_incorr=%.03f, p_incorr=%.03f ',modelCorr.Rsquared.Ordinary,modelCorr.Coefficients{2,4},...
        modelIncorr.Rsquared.Ordinary,modelIncorr.Coefficients{2,4}));
    xlabel(sprintf('Channel: %s in %s band',strrep(channelPair{1,2},'_',''),band))
    ylabel(sprintf('Channel: %s in %s band',strrep(channelPair{1,1},'_',''),band))
    
end
% modelIncorr = [];