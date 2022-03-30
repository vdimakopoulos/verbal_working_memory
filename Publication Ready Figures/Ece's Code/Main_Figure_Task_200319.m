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

strPaths.Variables = [strPaths.Klaver,'MATLAB Codes\Publication Figures\TASK\'];

% FieldTrip toolbox
strPaths.Toolboxes.FieldTrip            = 'D:\MATLAB Codes\Toolboxes\fieldtrip-20170925\';
% Matlab Import Import NLX files
strPaths.Toolboxes.MatlabImportExport   = 'D:\MATLAB Codes\Toolboxes\MatlabImportExport_v6.0.0\';

% Change main directory
cd(strPaths.Main)

% Add all subfolders to path
addpath(genpath(strPaths.Klaver))
addpath(genpath(strPaths.GeneralFunctions))

% Format for saving variables and images
strFormat.VariableNames.BehaviorAnalysis.ValidTrials = 'Behavior_Analysis_Results_Valid_Trials_Patient_%.2d_Session_%.2d_Part_%.2d.mat';

%% Load variables for the figure
Variables_Behavioral = load([strPaths.Variables,'Behavioral_Results']);
Variables_Anatomical = load([strPaths.Variables,'Anatomical_Locations']);

ImgTask = imread('D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\TASK\TaskIm.png');

%% Plot the main figure
if 1
    close all
    
    %%
    clear fig
    clear ax
    
    clear pos34_Phases
    clear strPhaseFigureColors
    clear plFixationMarker
    clear posAnatomicalY
    clear ttTaskPhases
    
    %%
    strFontName = 'Arial';
    nFontSize = 10;
    
    %%
    nMarkerSize_Behavior = 2;
    nLineWidth_Behavior = 0.5;
    
    %%
    aspectRatio_Axis_Screen = 1.4;
    pos34_Phases(1) = 3.488372093023256; % 5.2;
    pos34_Phases(2) = pos34_Phases(1)/aspectRatio_Axis_Screen;
    
    TaskAxOffsetX = 1.5;
    TaskInterAxDistance = pos34_Phases(1)/10;
    TaskAxOffsetY = 1.5;
    
    strPhaseFigureColors.Background = [0.9,0.9,0.9];
    strPhaseFigureColors.Objects(1,1:4) = {[2,165,85]/256,[251,224,1]/256,[241,26,39]/256,[4,142,205]/256}; % {'g','y','r','b'};
    strPhaseFigureColors.Objects(2,1:4) = {[2,165,85]/256,[251,224,1]/256,[241,26,39]/256,[0,0,0]}; % {'g','y','r','k'};
    
    posAnatomicalY = 3.3;
    
    nMarkerSize_Anatomical = 2;
    
    strTaskPhases = {'Fixation','Memory array','Retention interval','Probe array';'2-5 s','0.8 s','0.9 s','2 s'};
    
    %%
    fig = figure;
    fig.Units = 'centimeters';
    
    fig.Position = [-40,0,18*1,12];
    
    ax(1) = axes('Units','centimeters');
    ax(2) = axes('Units','centimeters');
    ax(3) = axes('Units','centimeters');
    ax(4) = axes('Units','centimeters');
    ax(5) = axes('Units','centimeters');
    ax(6) = axes('Units','centimeters');
    ax(7) = axes('Units','centimeters');
    
    ax(1).Position = [TaskAxOffsetX                                             fig.Position(4)-TaskAxOffsetY-pos34_Phases(2)   pos34_Phases(1),pos34_Phases(2)];
    ax(2).Position = [TaskAxOffsetX+1*pos34_Phases(1)+1*TaskInterAxDistance     fig.Position(4)-TaskAxOffsetY-pos34_Phases(2)   pos34_Phases(1),pos34_Phases(2)];
    ax(3).Position = [TaskAxOffsetX+2*pos34_Phases(1)+2*TaskInterAxDistance     fig.Position(4)-TaskAxOffsetY-pos34_Phases(2)   pos34_Phases(1),pos34_Phases(2)];
    ax(4).Position = [TaskAxOffsetX+3*pos34_Phases(1)+3*TaskInterAxDistance     fig.Position(4)-TaskAxOffsetY-pos34_Phases(2)   pos34_Phases(1),pos34_Phases(2)];
    
    ax(5).Position = [1.7,fig.Position(4)-16+6,3.9,3.9];
    ax(6).Position = [7.2,fig.Position(4)-16+6,3.9,3.9];
    
    ax(7).Position = [13 fig.Position(4)-16+6 length(Variables_Anatomical.pXAxSagittal)/length(Variables_Anatomical.pYAxSagittal)*posAnatomicalY posAnatomicalY];
    
    length(Variables_Anatomical.pXAxSagittal)
    length(Variables_Anatomical.pYAxSagittal)
    
    %% Panel a / Task phases
    % Plot the fixation point
    for ii = 1:4
        axes(ax(ii))
        plFixationMarker(ii) = plot(mean([0,pos34_Phases(1)]),mean([0,pos34_Phases(2)]),'Marker','+','Color','k');
        plFixationMarker(ii).LineWidth = 1;
    end
    
    %%
    if 0
        axes(ax(2))
        cla
        imagesc(linspace(0,pos34_Phases(1),size(ImgTask,1)),linspace(0,pos34_Phases(2),size(ImgTask,2)),flipud(ImgTask))
        set(gca,'YDir','normal')
    end
    
    %% Plot the squares
    axes(ax(2))
    hold on
    plSquares_1(1,1) = plot(1.8/5.2*pos34_Phases(1),2.64/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{1,1},'MarkerFaceColor',strPhaseFigureColors.Objects{1,1});
    plSquares_1(1,1).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_1(1,2) = plot(3.86/5.2*pos34_Phases(1),2.44/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{1,2},'MarkerFaceColor',strPhaseFigureColors.Objects{1,2});
    plSquares_1(1,2).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_1(2,1) = plot(1.31/5.2*pos34_Phases(1),1.5/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{1,3},'MarkerFaceColor',strPhaseFigureColors.Objects{1,3});
    plSquares_1(2,1).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_1(2,2) = plot(3.38/5.2*pos34_Phases(1),1.05/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{1,4},'MarkerFaceColor',strPhaseFigureColors.Objects{1,4});
    plSquares_1(2,2).MarkerSize = 19/5.2*pos34_Phases(1);
    
    axes(ax(4))
    hold on
    plSquares_2(1,1) = plot(1.8/5.2*pos34_Phases(1),2.64/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{2,1},'MarkerFaceColor',strPhaseFigureColors.Objects{2,1});
    plSquares_2(1,1).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_2(1,2) = plot(3.86/5.2*pos34_Phases(1),2.44/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{2,2},'MarkerFaceColor',strPhaseFigureColors.Objects{2,2});
    plSquares_2(1,2).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_2(2,1) = plot(1.31/5.2*pos34_Phases(1),1.5/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{2,3},'MarkerFaceColor',strPhaseFigureColors.Objects{2,3});
    plSquares_2(2,1).MarkerSize = 19/5.2*pos34_Phases(1);
    
    plSquares_2(2,2) = plot(3.38/5.2*pos34_Phases(1),1.05/3.7143*pos34_Phases(2),'Marker','s','MarkerEdgeColor',strPhaseFigureColors.Objects{2,4},'MarkerFaceColor',strPhaseFigureColors.Objects{2,4});
    plSquares_2(2,2).MarkerSize = 19/5.2*pos34_Phases(1);
    
    %%
    % Axis properties
    for ii = 1:4
        ax(ii).XLim = [0,pos34_Phases(1)];
        ax(ii).YLim = [0,pos34_Phases(2)];
        ax(ii).Box = 'on';
        ax(ii).XColor = 'k';
        ax(ii).YColor = 'k';
        ax(ii).XTick = [];
        ax(ii).YTick = [];
        ax(ii).XAxis.LineWidth = ceil(2/5.2*pos34_Phases(1)*2)/2;
        ax(ii).YAxis.LineWidth = ax(ii).XAxis.LineWidth;
        ax(ii).Color = strPhaseFigureColors.Background;
    end
    
    %% Add text
    for ii = 1:4
        axes(ax(ii))
        for jj = 1:2
            ttTaskPhases(ii,jj) = text(0,-1,strTaskPhases{jj,ii});
            ttTaskPhases(ii,jj).Units = 'normalized';
        end
    end
    for ii = 1:4
        for jj = 1:2
            ttTaskPhases(ii,jj).FontName = strFontName;
            ttTaskPhases(ii,jj).FontSize = nFontSize;
        end
    end
    for ii = 1:4
        for jj = 1:2
            ttTaskPhases(ii,jj).Position(1) = 0.02;
            switch jj
                case 1
                    ttTaskPhases(ii,jj).Position(2) = -0.12;
                case 2
                    ttTaskPhases(ii,jj).Position(2) = -0.30;
            end
        end
    end
    
    %% Plot K for all sessions / Publication
    axes(ax(5))
 
    [~,indSort] = sort(Variables_Behavioral.K_SS_Sessions(:,4));
    [~,indSort] = sort(indSort);
    
    c2 = 0.01;
    xoff = -(length(indSort)+1)*c2/2;
     
    hold on
    for iPatientSession = 1:size(Variables_Behavioral.nPatientSessionList,1)
        tempval = 0; % randn(1)/50;
        if(ismember(Variables_Behavioral.nPatientSessionList(iPatientSession,1),Variables_Behavioral.nPatientList_Capacity{1}))
            plot(Variables_Behavioral.SetSizeValues+c2*indSort(iPatientSession,:)+xoff+tempval,Variables_Behavioral.K_SS_Sessions(iPatientSession,:),'k:','Color',Variables_Behavioral.strPlotColors_MC{1},'LineWidth',nLineWidth_Behavior)
        else
            plot(Variables_Behavioral.SetSizeValues+c2*indSort(iPatientSession,:)+xoff+tempval,Variables_Behavioral.K_SS_Sessions(iPatientSession,:),'b:','Color',Variables_Behavioral.strPlotColors_MC{2},'LineWidth',nLineWidth_Behavior)
        end
        
        if(ismember(Variables_Behavioral.nPatientSessionList(iPatientSession,1),Variables_Behavioral.nPatientList_Capacity{1}))
            strMarker = 'v';
        else
            strMarker = 'o';
        end
        
        for iSS = 1:Variables_Behavioral.nNumberOfSetSizes
            pp = plot(Variables_Behavioral.SetSizeValues(iSS)+c2*indSort(iPatientSession,:)+xoff+tempval,Variables_Behavioral.K_SS_Sessions(iPatientSession,iSS),...
                Variables_Behavioral.strPlotColor_SS{iSS},'Marker',strMarker,'MarkerFaceColor',Variables_Behavioral.strPlotColor_SS{iSS}); % strPlotColor_SS{iSS});
            pp.MarkerSize = nMarkerSize_Behavior;
        end
        %     [valmax,indmax] = max(BehaviorResults{iPatientSession}.K_SS);
        %     plot(BehaviorResults{iPatientSession}.SetSizeValues(indmax)+tempval,valmax,['w',strMarker],'MarkerSize',4,'MarkerFaceColor','w')
    end
    
    for iSS = 1:Variables_Behavioral.nNumberOfSetSizes
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.K_SS_Sessions(:,iSS)),median(Variables_Behavioral.K_SS_Sessions(:,iSS))],...
            'b','LineWidth',nLineWidth_Behavior+0.5);
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.K_SS_Sessions(:,iSS))+std(Variables_Behavioral.K_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1),...
            median(Variables_Behavioral.K_SS_Sessions(:,iSS))+std(Variables_Behavioral.K_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1)],...
            'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.K_SS_Sessions(:,iSS))-std(Variables_Behavioral.K_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1),...
            median(Variables_Behavioral.K_SS_Sessions(:,iSS))-std(Variables_Behavioral.K_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1)],...
            'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    end
    
    set(gca,'XTick',Variables_Behavioral.SetSizeValues)
    set(gca,'YTick',1:6)
    set(gca,'TickDir','out')
    ylabel('Cowan''s K')
    xlabel('Set size')
    xlim([0.5,6.5])
    ylim([0,5])
    
    %% Plot response time for all sessions / Publication
    axes(ax(6))
%      figure
    c2 = 0.01;
    xoff = -(length(indSort)+1)*c2/2;
    
    vals_ToOrder = Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,end);
    [~,indSort] = sort(vals_ToOrder);
    [~,indSort] = sort(indSort);
    
    hold on
    for iPatientSession = 1:size(Variables_Behavioral.nPatientSessionList,1)
        tempval = randn(1)/20000;
        if(ismember(Variables_Behavioral.nPatientSessionList(iPatientSession,1),Variables_Behavioral.nPatientList_Capacity{1}))
            plot(Variables_Behavioral.SetSizeValues+tempval+c2*indSort(iPatientSession,:)+xoff,Variables_Behavioral.ResponseTime_Median_SS_Sessions(iPatientSession,:),'k:','Color',Variables_Behavioral.strPlotColors_MC{1},'LineWidth',nLineWidth_Behavior)
        else
            plot(Variables_Behavioral.SetSizeValues+tempval+c2*indSort(iPatientSession,:)+xoff,Variables_Behavioral.ResponseTime_Median_SS_Sessions(iPatientSession,:),'b:','Color',Variables_Behavioral.strPlotColors_MC{2},'LineWidth',nLineWidth_Behavior)
        end
        
        if(ismember(Variables_Behavioral.nPatientSessionList(iPatientSession,1),Variables_Behavioral.nPatientList_Capacity{1}))
            strMarker = 'o';
        else
            strMarker = 'o';
        end
        
        for iSS = 1:Variables_Behavioral.nNumberOfSetSizes
            pp = plot(Variables_Behavioral.SetSizeValues(iSS)+tempval+c2*indSort(iPatientSession,:)+xoff,Variables_Behavioral.ResponseTime_Median_SS_Sessions(iPatientSession,iSS),...
                Variables_Behavioral.strPlotColor_SS{iSS},'Marker',strMarker,'MarkerFaceColor',Variables_Behavioral.strPlotColor_SS{iSS}); % strPlotColor_SS{iSS});
            pp.MarkerSize = nMarkerSize_Behavior;
        end
        %     [valmax,indmax] = max(BehaviorResults{iPatientSession}.K_SS);
        %     plot(BehaviorResults{iPatientSession}.SetSizeValues(indmax)+tempval,valmax,['w',strMarker],'MarkerSize',4,'MarkerFaceColor','w')
    end
    
    for iSS = 1:Variables_Behavioral.nNumberOfSetSizes
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS)),median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))],...
            'b','LineWidth',nLineWidth_Behavior+0.5);
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))+std(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1),...
            median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))+std(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1)],...
            'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
        plot([Variables_Behavioral.SetSizeValues(iSS)-xoff*2,Variables_Behavioral.SetSizeValues(iSS)+xoff*2],...
            [median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))-std(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1),...
            median(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))-std(Variables_Behavioral.ResponseTime_Median_SS_Sessions(:,iSS))/size(Variables_Behavioral.nPatientSessionList,1)],...
            'b','LineWidth',nLineWidth_Behavior,'Color',[0.4,0.4,1]);
    end
    
    set(gca,'XTick',[1,2,4,6])
    set(gca,'YTick',0:0.5:3)
    set(gca,'TickDir','out')
    ylabel('Response time (s)')
    xlabel('Set size')
    xlim([0.5,6.5])
    ylim([0.4,1.6])
    
    %% Load new images for plotting anatomical locations
    axes(ax(7))
 
    Locations_ToPlot = Variables_Anatomical.AnatomicalLocation_Table_FirstContact_Regions;
    
    imagesc(Variables_Anatomical.pXAxSagittal,Variables_Anatomical.pYAxSagittal,Variables_Anatomical.ImgBrain)
    set(gca,'YDir','Normal')
    
    set(gca,'XLim',[Variables_Anatomical.pXAxSagittal(1),Variables_Anatomical.pXAxSagittal(end)])
    set(gca,'YLim',[Variables_Anatomical.pYAxSagittal(end),Variables_Anatomical.pYAxSagittal(1)])
    
    hold on
    s = [];
    for nRegion = 1:3
        Temp = Locations_ToPlot{nRegion};
        s = scatter((Temp.y),...
            (Temp.z),'b','Marker','o','MarkerFaceColor',Variables_Anatomical.strPlotColors_Regions{nRegion},'MarkerEdgeColor',Variables_Anatomical.strPlotColors_Regions{nRegion});
        s.SizeData = nMarkerSize_Anatomical;
    end
   
%     ll = legend('Hipp','Ent','Amg');
    
    % tt = [];
    % tt(1) = text(-40,-16,'Hipp','Color',strPlotColors_Regions{1});
    % tt(2) = text(2,-38,'Ent','Color',strPlotColors_Regions{2});
    % tt(3) = text(-1.8,-12,'Amg','Color',strPlotColors_Regions{3});
    % tt(4) = text(-40,-25,'Other','Color',strPlotColors_Regions{4});
    
    xlabel('') % xlabel('MNI y-coordinate (mm)')
    ylabel('') % ylabel('MNI z-coordinate (mm)')
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    
    %% Add subplot labels
    axes(ax(1))
    ttLetter(1) = text(2,2,'a');
    axes(ax(5))
    ttLetter(2) = text(2,2,'b');
    axes(ax(6))
    ttLetter(3) = text(2,2,'c');
    axes(ax(7))
    ttLetter(4) = text(2,2,'d');
    
    for ii = 1:4
        ttLetter(ii).Units = 'centimeters';
        ttLetter(ii).Position(1) = -0.5;
        switch ii
            case 1
                ttLetter(ii).Position(2) = ax(1).Position(4)+0.4;
            case {2,3,4}
                ttLetter(ii).Position(2) = ax(ii+3).Position(4)+0.4;
        end
    end
    
    for ii = 1:4
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
    for ii = 1:length(ax)
        ax(ii).XColor = 'k';
        ax(ii).YColor = 'k';
    end
    
    %% Lengthen the tick marks
    ax(5).TickLength(1) = 0.02;
    ax(6).TickLength(1) = 0.02;
    
    %% Re-adjust axis y-labels
    ax(5).YLabel.Units = 'centimeter';
    ax(5).YLabel.VerticalAlignment = 'middle';
    ax(5).YLabel.Position(1) = -1.05;
    
    ax(6).YLabel.Units = 'centimeter';    
    ax(6).YLabel.VerticalAlignment = 'middle';
    ax(6).YLabel.Position(1) = -1.15;
    
    %% Re-adjust letters  
    ttLetter(2).Position(1) = ax(5).YLabel.Position(1);
    ax(5).Position(1) = ax(1).Position(1)+ttLetter(1).Position(1)-ttLetter(2).Position(1);
    
    ax(6).Position(1) = ax(5).Position(1)+ax(5).Position(3)+1.8;
    ax(7).Position(1) = ax(6).Position(1)+ax(6).Position(3)+0.9;
    
    %%
    ax(7).Position(2) = ax(6).Position(2)+ax(6).Position(4)-ax(7).Position(4);
    
end

%%
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\TASK\Files\Figure_TASK_RGB.fig')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\TASK\Files\Figure_TASK_RGB.png')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\TASK\Files\Figure_TASK_RGB.pdf')
saveas(fig,'D:\MATLAB Codes\Klaver Task\MATLAB Codes\Publication Figures\TASK\Files\Figure_TASK_RGB.eps')