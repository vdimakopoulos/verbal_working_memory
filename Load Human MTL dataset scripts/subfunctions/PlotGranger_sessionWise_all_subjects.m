load('F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Granger_SessionWise.mat');

nSubjs = length(Granger_SessionWise);
figure;

% Markers = {'o','+', '*', '.', 'x', '_', '|', 's', 'd'}
maxGranger_MaintSession = [];
CorrectReponseSession = [];
for i =1:nSubjs
    
    scatter(Granger_SessionWise{i}.maxGranger_MaintSession*100,Granger_SessionWise{i}.CorrectReponseSession*100,30);
    hold on
    maxGranger_MaintSession = [maxGranger_MaintSession Granger_SessionWise{i}.maxGranger_MaintSession]
    CorrectReponseSession = [CorrectReponseSession Granger_SessionWise{i}.CorrectReponseSession]
    
end
set(gca,'box','off','FontSize',16);

model = fitlm(maxGranger_MaintSession*100,CorrectReponseSession*100,'RobustOpts','on');
plot(model)
hLegend = findobj(gcf, 'Type', 'Legend');
title('');

set(hLegend,'Visible','off')
xlabel('?Granger (%)');
ylabel('Correct Response Accuracy (%)');
