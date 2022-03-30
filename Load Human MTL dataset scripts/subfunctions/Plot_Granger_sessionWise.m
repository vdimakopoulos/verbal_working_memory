function [maxGranger_Maint ] = ...
    Plot_Granger_sessionWise(gdata_session,DeltaGranger_Maint,CorrectResponseRate_Session,...
    SignificantBars_GrangerScalp_Maint,selected_subj,DeltaGranger_currentSubj,SubjectAccuracy,nTrial_per_Session)

GrangerSessions_maint = gdata_session.Maint;
nSessions = length(GrangerSessions_maint);
indMaxFreq = SignificantBars_GrangerScalp_Maint{selected_subj};

for iSes = 1:nSessions
   maxGranger_Maint(iSes) = max(DeltaGranger_Maint{iSes});
end

figure;
plot(maxGranger_Maint*100,CorrectResponseRate_Session*100,'o','MarkerEdgeColor','red',...
        'MarkerFaceColor','red');

for iSes = 1:nSessions

text(maxGranger_Maint(iSes)*100,CorrectResponseRate_Session(iSes)*100,['\leftarrow' sprintf('Session %d - Trials: %d',iSes,nTrial_per_Session(iSes))]);

end
hold on;
plot(DeltaGranger_currentSubj,SubjectAccuracy,'o','MarkerEdgeColor','k',...
    'MarkerFaceColor','k');

ylabel('Correct Response Rate (%)');
xlabel('\DeltaGranger (%)');
set(gca,'box','off','FontSize',16);

end