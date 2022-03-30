%% Close all figures, clear variables and command window
close all
clear
clc

%% Paths
% Folder that contains the folder 'Klaver Task'
strPaths.Main = 'D:\MATLAB Codes\';
strPaths.GeneralFunctions = [strPaths.Main,'General Functions\'];

% Main folder of task and subfolders
strPaths.Klaver                                     = [strPaths.Main,'Klaver Task\'];

strPaths.Variables = [strPaths.Klaver,'MATLAB Codes\Publication Figures\MAINT\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'D:\MATLAB Codes\Toolboxes\fieldtrip-20170925\';
% Matlab Import Import NLX files
strPaths.Toolboxes.MatlabImportExport   = 'D:\MATLAB Codes\Toolboxes\MatlabImportExport_v6.0.0\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(genpath(strPaths.Klaver))
addpath(genpath(strPaths.GeneralFunctions))

% Add FieldTrip and toolbox for importing NLX to path
addpath(strPaths.Toolboxes.FieldTrip)
ft_defaults

% Format for saving variables and images
strFormat.VariableNames.BehaviorAnalysis.ValidTrials = 'Behavior_Analysis_Results_Valid_Trials_Patient_%.2d_Session_%.2d_Part_%.2d.mat';

%% Load variables for the figure
Variables_SPIKES_FigureMaint = load([strPaths.Variables,'SPIKES_FigureMaint']);
Variables_MaintUnits_Numbers = load([strPaths.Variables,'MaintUnits_Numbers']);
Variables_LoadDep_Numbers = load([strPaths.Variables,'LoadDep_Numbers']);

%% Significance of firing rate
flagComputeSignificance = 0;
if(flagComputeSignificance)
    % fig = Main_Plot_Single_Spike_Shape_Raster_Rate_All_181212(SPIKES_FigureMaint);
    [Variables_SPIKES_FigureMaint.stat_Maint_HL,Variables_SPIKES_FigureMaint.t_axis_Maint_HL] = ...
        Get_FT_Significance_Analysis(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint,1);
    save([strPaths.Variables,'stat_Maint_HL'],'stat_Maint_HL','t_axis_Maint_HL')
else
    Variables_SPIKES_FigureMaint.stat_Maint_HL = load([strPaths.Variables,'stat_Maint_HL'],'stat_Maint_HL');
    Variables_SPIKES_FigureMaint.stat_Maint_HL = Variables_SPIKES_FigureMaint.stat_Maint_HL.stat_Maint_HL;
    Variables_SPIKES_FigureMaint.t_axis_Maint_HL = load([strPaths.Variables,'stat_Maint_HL'],'t_axis_Maint_HL');
    Variables_SPIKES_FigureMaint.t_axis_Maint_HL = Variables_SPIKES_FigureMaint.t_axis_Maint_HL.t_axis_Maint_HL;
end
fig = Main_Plot_Single_Spike_Shape_Raster_Rate_All_181212(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint);
% Linear fit of firing rates
% Significance bar
figure
hax = axes('Position',[0.1300,0.73,0.4942,0.01]);
ind = find(Variables_SPIKES_FigureMaint.stat_Maint_HL.mask);
conn = bwconncomp(Variables_SPIKES_FigureMaint.stat_Maint_HL.mask);
for ii = 1:length(conn)
    ind2 = conn.PixelIdxList{ii};
    rr = rectangle('Position',[Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(1)),-1,Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(end))-Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(1)),1]);
    set(rr,'FaceColor','k')
end
set(gca,'XTick',[])
set(gca,'YTick',[])
xlim([-1.7000 2.9000])
ttSig = text(-1.85,0,'Sig','Rotation',90,'HorizontalAlignment','center','FontSize',6);
box on

%% Spike shapes
Variables_SPIKES_FigureMaint.meanWave = mean(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.all_spike_wave);
Variables_SPIKES_FigureMaint.stdWave = std(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.all_spike_wave);
Variables_SPIKES_FigureMaint.fs = 32000;

%% Plot the main figure
if 1
    close all
    %%
    nListOfSS = [1,2,4,6];
    strListOfColors_SS = {'b','g','r','k'};
    
    t_min = -1.7-0.05;
    t_max = 2.9+0.05;
    
    tTrialPeriods = [-0.8,0,0.9];
    strTrialPeriods = {'Fixation','Encoding','Maintenance','Test'};
    
    strListOfGroups = {'Hipp','Ent','Amg'};
    
    %%
    clear fig
    clear ax
    
    clear plRates
    clear plTrialPhaseLines
    clear ttTrialPhase
    clear ttSetSizes
    clear plShape
    clear tt_stars
    
    %%
    strFontName = 'Arial';
    nFontSize = 10;
    
    %%
    fig = figure;
    fig.Units = 'centimeters';
    
    fig.Position = [-26,-4,18*1,24];
    
    ax(1) = axes('Units','centimeters');
    ax(2) = axes('Units','centimeters');
    ax(3) = axes('Units','centimeters');
    ax(4) = axes('Units','centimeters');
    ax(5) = axes('Units','centimeters');
    ax(6) = axes('Units','centimeters');
    ax(7) = axes('Units','centimeters');
    
    ax(1).Position = [2,fig.Position(4)-1.5-5.5,14,5.5];
    ax(2).Position = [2,fig.Position(4)-1.5-5.5-5.5-1,14,5.5];
    ax(3).Position = [2    fig.Position(4)-26+4.9    5.9761    5.4593];
    ax(4).Position = [10   fig.Position(4)-26+4.9    5.9761    5.4593];
    ax(5).Position = [13.75,fig.Position(4)-26+23.5 1.4 1.4];
    ax(6).Position = [ax(1).Position(1) fig.Position(4)-26+18.3 ax(1).Position(3) 0.4];
    ax(7).Position = [ax(4).Position(1)+1 ax(4).Position(2)+1 1,1];
        
    %% Panel a / Spike rates
    axes(ax(1))
    plot(NaN)
    hold on
    SpikeRates = Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.SpikeRates_Trials_All;
    for iSS = 1:length(nListOfSS)
        nSS = nListOfSS(iSS);
        indTrials = Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.TrialInformationTable.SetSize==nSS;
        SpikeRates_SS = SpikeRates(indTrials,:);
        meanSpikeRates   = nanmean(SpikeRates_SS);
        stdSpikeRates    = nanstd(SpikeRates_SS)/sqrt(size(SpikeRates_SS,1));
        if(isempty(SpikeRates_SS))
            meanSpikeRates = zeros(size(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.t_bin_center));
            stdSpikeRates = zeros(size(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.t_bin_center));
        end
        plRates{iSS} = shadedErrorBar(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.t_bin_center,meanSpikeRates,stdSpikeRates,{strListOfColors_SS{iSS},'linewidth',2});
        plRates{iSS}.edge.delete;
    end
    xlim([t_min,t_max])
    ylim([0,11])
    ylabel('Firing rate (Hz)')
    
    box off
    
    ax(1).XAxis.Visible = 'off';
    
    % Trial phases
    ylims = get(gca,'YLim');
    ylims(1) = 0;
    set(gca,'YLim',ylims)
    hold on
    for iTrialPhase = 1:length(tTrialPeriods)
        plTrialPhaseLines(iTrialPhase,1) = plot(tTrialPeriods(iTrialPhase)*[1,1],ylims,'Color',0.4*[1,1,1]);
    end
    % Names of trial phases
    for iTrialPhase = 1:length(strTrialPeriods)
        switch iTrialPhase
            case 1
                ttTrialPhase(iTrialPhase) = text(mean([t_min,tTrialPeriods(iTrialPhase)]),ylims(1)+(ylims(2)-ylims(1))*1.1,...
                    strTrialPeriods{iTrialPhase});
            case {2,3}
                ttTrialPhase(iTrialPhase) = text(mean([tTrialPeriods(iTrialPhase-1),tTrialPeriods(iTrialPhase)]),ylims(1)+(ylims(2)-ylims(1))*1.1,...
                    strTrialPeriods{iTrialPhase});
            case 4
                ttTrialPhase(iTrialPhase) = text(tTrialPeriods(iTrialPhase-1)+0.4,ylims(1)+(ylims(2)-ylims(1))*1.1,...
                    strTrialPeriods{iTrialPhase});
        end
        ttTrialPhase(iTrialPhase).HorizontalAlignment = 'center';
    end
    for ii = 1:4
        ttTrialPhase(ii).FontSize = nFontSize;
        ttTrialPhase(ii).FontName = strFontName;
    end
    for ii = 1:4
        ttTrialPhase(ii).Units = 'centimeters';
        ttTrialPhase(ii).Position(2) = ax(1).Position(4)+0.3;
    end
    for ii = 1:4
        ttSetSizes(ii) = text(-1.5,9.5-(ii-1)*0.7,['Set size ',num2str(nListOfSS(ii))],'Color',strListOfColors_SS{ii});
    end
    for ii = 1:4
        ttSetSizes(ii).FontSize = nFontSize;
        ttSetSizes(ii).FontName = strFontName;
        ttSetSizes(ii).VerticalAlignment = 'middle';
    end
    
    ax(1).Color = 'none';

    %% Panel a / Significance bar
    axes(ax(6))
    ind = find(Variables_SPIKES_FigureMaint.stat_Maint_HL.mask);
    conn = bwconncomp(Variables_SPIKES_FigureMaint.stat_Maint_HL.mask);
    for ii = 1:length(conn)
        ind2 = conn.PixelIdxList{ii};
        rr = rectangle('Position',[Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(1)),-1,Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(end))-Variables_SPIKES_FigureMaint.stat_Maint_HL.time(ind2(1)),1]);
        set(rr,'FaceColor','k')
    end
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    xlim([t_min,t_max])
    box on
    
    %% Panel a / Spike raster
    axes(ax(2))
    y_max = 193;
    spikes = Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.spikes(~isnan(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.spikes(:,2)),:);
    % Number of trials for each set size
    nNumberOfTrials_SS = zeros(length(nListOfSS),1);
    for iSS = 1:length(nListOfSS)
        nSS = nListOfSS(iSS);
        nNumberOfTrials_SS(iSS) = length(find(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.TrialInformationTable.SetSize==nSS));
    end
    % Sort spikes using set sizes
    for iSS = 1:length(nListOfSS)
        nSS = nListOfSS(iSS);
        nTrials_SS = find(Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.TrialInformationTable.SetSize==nSS);
        for nTrial = 1:length(nTrials_SS)
            spikes(spikes(:,2)==nTrials_SS(nTrial),4) = nTrial+sum(nNumberOfTrials_SS(1:iSS-1));
        end
    end
    % Plot raster
    for iSS = 1:length(nListOfSS)
        nSS = nListOfSS(iSS);
        indTrials = Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.TrialInformationTable.SetSize==nSS;
        spikes_SS = spikes(ismember(spikes(:,2),find(indTrials)),:);
        scatter(spikes_SS(:,3),spikes_SS(:,4),25,[strListOfColors_SS{iSS},'.'])
        hold on
    end
    xlabel('Time (s)')
    ylabel('Trial no. (reordered)')
    axis([t_min, t_max 0 y_max])
    
    box off
    % Trial phases
    ylims = get(gca,'YLim');
    for iTrialPhase = 1:length(tTrialPeriods)
        plTrialPhaseLines(iTrialPhase,2) = plot(tTrialPeriods(iTrialPhase)*[1,1],ylims,'Color',0.4*[1,1,1]);
    end
    
    %% Panel a / Spike shape
    axes(ax(5))
    
    t_shift = 0.5938; % ms
    plShape = shadedErrorBar(((1:64)-1)/Variables_SPIKES_FigureMaint.fs*1000-t_shift,Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.spike_wave(:,1),Variables_SPIKES_FigureMaint.SPIKES_FigureMaint.spike_wave(:,2));
    plShape.edge.delete;
    xlim([-t_shift 64/Variables_SPIKES_FigureMaint.fs*1000-t_shift])
    ylim([-17,40])
    set(gca,'XTick',[0,1])
    set(gca,'YTick',[0,20])
    axis square
    box off
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    
    %% Panel b
    axes(ax(3))
    strListOfGroups = {'Hipp','Ent','Amg'};
    
    bar(Variables_MaintUnits_Numbers.nInc./Variables_MaintUnits_Numbers.nAll*100,'FaceColor',0.4*[1,1,1])
    set(gca,'XTickLabel','')
    ylabel('Percentage of all units')
    xlim([-0,4])
    ylim([0,100])
    set(gca,'XTick',[])
    set(gca,'YTick',0:20:100)
    box off
    
    level1 = -5;
    tt_labels(1,1) = text(1,level1,strListOfGroups{1},'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(2,1) = text(2,level1,strListOfGroups{2},'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(3,1) = text(3,level1,strListOfGroups{3},'HorizontalAlignment','center','FontSize',nFontSize);
    
    level2 = -12;
    tt_labels(1,2) = text(1,level2,sprintf('(%d/%d)',Variables_MaintUnits_Numbers.nInc(1),Variables_MaintUnits_Numbers.nAll(1)),'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(2,2) = text(2,level2,sprintf('(%d/%d)',Variables_MaintUnits_Numbers.nInc(2),Variables_MaintUnits_Numbers.nAll(2)),'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(3,2) = text(3,level2,sprintf('(%d/%d)',Variables_MaintUnits_Numbers.nInc(3),Variables_MaintUnits_Numbers.nAll(3)),'HorizontalAlignment','center','FontSize',nFontSize);
    
    level3 = Variables_MaintUnits_Numbers.nInc./Variables_MaintUnits_Numbers.nAll*100+3;
%     tt_labels(1,3) = text(1,level3(1),[num2str(Variables_MaintUnits_Numbers.nInc(1)/Variables_MaintUnits_Numbers.nAll(1)*100,'%.1f'),'%'],'HorizontalAlignment','center','FontSize',nFontSize,'VerticalAlignment','bottom');
%     tt_labels(2,3) = text(2,level3(2),[num2str(Variables_MaintUnits_Numbers.nInc(2)/Variables_MaintUnits_Numbers.nAll(2)*100,'%.1f'),'%'],'HorizontalAlignment','center','FontSize',nFontSize,'VerticalAlignment','bottom');
%     tt_labels(3,3) = text(3,level3(3),[num2str(Variables_MaintUnits_Numbers.nInc(3)/Variables_MaintUnits_Numbers.nAll(3)*100,'%.1f'),'%'],'HorizontalAlignment','center','FontSize',nFontSize,'VerticalAlignment','bottom');
    
    for ii = 1:6
        tt_labels(ii).FontName = strFontName;
    end

    % Put stars on top of bars
    tt_stars(1,1) = text(1,level3(1),'***','HorizontalAlignment','center','FontSize',nFontSize+2);
    tt_stars(2,1) = text(2,level3(2),'***','HorizontalAlignment','center','FontSize',nFontSize+2);
    tt_stars(3,1) = text(3,level3(3),'***','HorizontalAlignment','center','FontSize',nFontSize+2);
    
    for ii = 1:3
        tt_stars(ii,1).FontName = strFontName;
    end
    
    %% Panel c
    axes(ax(4))
    bar(Variables_LoadDep_Numbers.nAll./Variables_LoadDep_Numbers.nAll*100,'FaceColor',0.5*[1,1,1])
    hold on
    bar((Variables_LoadDep_Numbers.nInc+Variables_LoadDep_Numbers.nDec)./Variables_LoadDep_Numbers.nAll*100,'FaceColor','r')
    bar(Variables_LoadDep_Numbers.nDec./Variables_LoadDep_Numbers.nAll*100,'FaceColor','b')
    set(gca,'XTickLabel','')
    ylabel('Percentage of maintenance units')
    xlim([-0,4])
    ylim([0,100])
    set(gca,'XTick',[])
    set(gca,'YTick',0:20:100)
    box off
    
    level1 = -5;
    tt_labels(1,1) = text(1,level1,strListOfGroups{1},'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(2,1) = text(2,level1,strListOfGroups{2},'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(3,1) = text(3,level1,strListOfGroups{3},'HorizontalAlignment','center','FontSize',nFontSize);
    
    level2 = -12;
    tt_labels(1,2) = text(1,level2,['(',num2str(Variables_LoadDep_Numbers.nInc(1)),' - ',num2str(Variables_LoadDep_Numbers.nDec(1)),')'],'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(2,2) = text(2,level2,['(',num2str(Variables_LoadDep_Numbers.nInc(2)),' - ',num2str(Variables_LoadDep_Numbers.nDec(2)),')'],'HorizontalAlignment','center','FontSize',nFontSize);
    tt_labels(3,2) = text(3,level2,['(',num2str(Variables_LoadDep_Numbers.nInc(3)),' - ',num2str(Variables_LoadDep_Numbers.nDec(3)),')'],'HorizontalAlignment','center','FontSize',nFontSize);
    
    for ii = 1:6
        tt_labels(ii).FontName = strFontName;
    end
    
    tt_stars(1,2) = text(1,(Variables_LoadDep_Numbers.nInc(1)+Variables_LoadDep_Numbers.nDec(1))./Variables_LoadDep_Numbers.nAll(1)*100+3,'***','HorizontalAlignment','center','FontSize',nFontSize+2);
    tt_stars(2,2) = text(2,(Variables_LoadDep_Numbers.nInc(2)+Variables_LoadDep_Numbers.nDec(2))./Variables_LoadDep_Numbers.nAll(2)*100+3,'***','HorizontalAlignment','center','FontSize',nFontSize+2);
    %     tt_stars(3,1) = text(3,nSelected(3)./nTotal(3)*100+3,'***','HorizontalAlignment','center','FontSize',14);
    
    for ii = 1:2
        tt_stars(ii,2).FontName = strFontName;
    end

    %% Legends
    axes(ax(7))
    strPlotColors_Change = {0.5*[1,1,1],'r','b'};
    ranClass_ToPlot = 1;
    ranRegions_ToPlot = [10,1,2,3];
    
    hold on
    strLegendLabels = {'No difference with workload','Increase with workload','Decrease with workload'};
    for ii = 1:length(strPlotColors_Change)
        pp = plot(NaN,'s','Color',strPlotColors_Change{ii},'Marker','s','MarkerFaceColor',strPlotColors_Change{ii},'MarkerSize',10);
    end
    ll = legend(strLegendLabels);
    
    %%
    ll.Box = 'off';
    ll.FontName = strFontName;
    ll.FontSize = nFontSize;
    ll.Units = 'centimeters';
    ll.Position(1) = 9.7;
    ll.Position(2) = ax(4).Position(2)-2.5;
    
    %% Tick directions
    for ii = 1:length(ax)
        ax(ii).TickDir = 'out';
    end
    ax(3).TickLength(1) = 0.02;
    ax(4).TickLength(1) = 0.02;
    ax(5).TickLength(1) = 0.02;
    
    %% Add subplot labels
    axes(ax(1))
    ttLetter(1) = text(2,2,'a');
    axes(ax(3))
    ttLetter(2) = text(2,2,'b');
    axes(ax(4))
    ttLetter(3) = text(2,2,'c');
    
    for ii = 1:3
        ttLetter(ii).Units = 'centimeters';
        ttLetter(ii).Position(1) = -0.5;
        switch ii
            case 1
                ttLetter(ii).Position(2) = ax(1).Position(4)+0.4;
            case {2,3}
                ttLetter(ii).Position(2) = ax(ii+1).Position(4)+0.4;                
        end
    end
    
    for ii = 1:3
        ttLetter(ii).VerticalAlignment = 'bottom';
        ttLetter(ii).HorizontalAlignment = 'center';
        ttLetter(ii).FontName = strFontName;
        ttLetter(ii).FontSize = nFontSize+1;
        ttLetter(ii).FontWeight = 'bold';
    end
    
    %%
    for ii = 1:length(ax)
        ax(ii).FontSize = nFontSize;
        ax(ii).FontName = strFontName;
    end
    
    %%
    ax(5).FontSize = nFontSize-3;
    ax(5).YLabel.Units = 'centimeters';
    ax(5).YLabel.Position(1) = -0.42;
    
    %%
    for ii = 1:4
    ax(ii).YLabel.VerticalAlignment = 'middle';
    ax(ii).YLabel.Units = 'centimeters';
    ax(ii).YLabel.Position(1) = -1.2;
    end
    for ii = 1:3
        ttLetter(ii).Position(1) = ax(ii).YLabel.Position(1);
    end
    
    %%
    for ii = 1:length(ax)
        ax(ii).XColor = 'k';
        ax(ii).YColor = 'k';
    end
end

%%
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\MAINT\Files\Figure_MAINT_RGB.fig')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\MAINT\Files\Figure_MAINT_RGB.png')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\MAINT\Files\Figure_MAINT_RGB.pdf')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\MAINT\Files\Figure_MAINT_RGB.eps')