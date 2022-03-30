function Scalp_topoplot(toolbox_path,freqBand,freqAxis,data,dataBipolar_Sclp,iSS_to_Plot,nChannelPairs,DataFlag)

%% Read electrodes from excel
Coords_excel_path = 'F:\Vasileios\Task Analysis\Data\EEGLAB_Data\Cap_coords_all.xls';
strSheetName = '64-chan added channels';
ElectrodeCoordinatesTable = Get_BioSemi_Coordinates_From_Excel( Coords_excel_path, strSheetName);
Project_Cartesian_3D_On_2D(ElectrodeCoordinatesTable.X, ElectrodeCoordinatesTable.Y, ElectrodeCoordinatesTable.Z);

%% Keep only the channels needed
strChannelLabels = dataBipolar_Sclp.label;
strChannelLabels = strrep(strChannelLabels,' ',''); 

ScalpChansPositions = ElectrodeCoordinatesTable(ismember(ElectrodeCoordinatesTable.label,strChannelLabels),:);

%% Get Chan Locations in EEGLAB format
rmpath(toolbox_path.FieldTrip);
addpath(genpath(toolbox_path.EEGLAB));
chanlocs  = Get_Standard_BESA_Coordinates_From_EEGLAB(toolbox_path.EEGLAB)
for i = 1:size(ScalpChansPositions,1)
    if ~(strcmp(ScalpChansPositions.label{i},'T1') ||  strcmp(ScalpChansPositions.label{i},'T2'))
        ind(i) = find(strcmp({chanlocs.labels},ScalpChansPositions.label{i}));
%         ind(i) = find(strcmp({chanlocs.labels},dataBipolar_Sclp.label{i})==1);

    end
end
T1_loc = find(contains(dataBipolar_Sclp.label, 'T1')==1);
if T1_loc~=0
    val = setdiff([1:size(ind,2)],[size(ind,2)-4+1:size(ind,2)]) %% here we only exclude A1,A2 because T1,T2 were excluded from the preprocessing
else
    val = setdiff([1:size(ind,2)],[size(ind,2)-2+1:size(ind,2)]) %% here we only exclude A1,A2 because T1,T2 were excluded from the preprocessing
end
ScalpChansLocs = chanlocs(ind(val)) %Exclude T1,T2,A1,A2


%% Prepare the EEGLAB datavector
[~,indFreq1] = min(abs(freqAxis-freqBand(1)));
[~,indFreq2] = min(abs(freqAxis-freqBand(2)));
% topoorder = [10 6 8 3 20 15 22 13 17 11 12 7 9 5 4 21 16 23 14 18 19 1 2]; %Map of chanlocs to dataBipolar_Sclp
% topoorder = [10 6 8 3 20 15 13 17 11 12 7 9 5 4 16 14 18 19 1 2];
A = dataBipolar_Sclp.label;
A= strtrim(A);
for i =1:numel(ScalpChansLocs)
    topoorder(i) = find(strcmp(ScalpChansLocs(i).labels,A));
end
for i =1:length(topoorder)
    datavector(i) = data(topoorder(i));%data{topoorder(i)};
end


topoplot(datavector,ScalpChansLocs,'electrodes','labels', .....
    'colormap','jet','maplimits',[0 0.02],'conv','on','efontsize',14); % For PLV
%             'colormap','jet','maplimits',[0 40],'conv','on');% For PSD

                   

                
 
        h = colorbar;
        set(h,'FontSize',13);
% 
% if DataFlag == 1
%     
% %     PLV_pairs_ss = data{iSS_to_Plot};
%     PLV_In_Band_SS = [];
%     for iPair = 1:size(data,1)
%         for iSS = 1:size(data,2)
%             PLV_temp = data{iPair,iSS}(indFreq1:indFreq2);
%             PLV_temp = abs(PLV_temp);
%             PLV_In_Band_SS(iPair,iSS) = mean(PLV_temp);
%         end
%     end
%     
%     %% Use EEGLAB's topoplot
%     for i=1:size(nChannelPairs,1)
%         if i <= size(dataBipolar_Sclp.label,1)
%             datavector(i,iSS_to_Plot) = PLV_In_Band_SS(topoorder(i),iSS_to_Plot); %%tbc
%         else
%             j = i-(floor((i-1)/23)*23); 
%             disp(num2str([topoorder(j)*round((i-1)/23) j]))
%             datavector(j,iSS_to_Plot) = PLV_In_Band_SS(topoorder(j)+round((i-1)/23)*23,iSS_to_Plot);
% %             figure;
% %             topoplot(datavector,ScalpChansLocs,'electrodes','labels', ...
% %                 'colormap','jet','maplimits','maxmin','conv','on');
% %             h = colorbar;
% %             set(h,'FontSize',18);
%         end
%         if mod(i,size(dataBipolar_Sclp.label,1)) == 0
%                DataVector{i/23} = datavector(1:19,iSS_to_Plot)'; % Exclude T1,T2,A1,A2
%                clear datavector;
%         end
%     end
%     
%     for i = 1:size(DataVector,2)
%         figure;
%         topoplot(DataVector{i},ScalpChansLocs,'electrodes','labels', ...
%             'colormap','jet','maplimits','maxmin','conv','on');
%         h = colorbar;
%         set(h,'FontSize',18);
%         
%     end
% else
%     for i = 1:length(order)
%         datavector(i) = mean(10*log10(data.powspctrm(topoorder(i),indFreq1:indFreq2)));
%     end
%     datavector = datavector(1:19); % Exclude T1,T2,A1,A2
%     figure;
%     topoplot(datavector,ScalpChansLocs,'electrodes','labels', ...
%                 'colormap','jet','maplimits','maxmin','conv','on');
%      h = colorbar;
%      set(h,'FontSize',18);
% end
%% Reset Paths
rmpath(genpath(toolbox_path.EEGLAB));
addpath(toolbox_path.FieldTrip);

% % 
% % ScalpChansLocs = table2struct(ScalpChansPositions)';
% % % Fix field names
% % [ScalpChansLocs.labels] = ScalpChansLocs.label;
% % ScalpChansLocs = rmfield(ScalpChansLocs,'label');
% % r = num2cell(sqrt([ScalpChansLocs.x].^2+[ScalpChansLocs.y].^2)/85);
% % [ScalpChansLocs.radius] = r{:}
% % ScalpChansLocs = rmfield(ScalpChansLocs,{'theta','phi','r','x','y'});
% % 
% % 
% % str = strsplit(repmat('EEG ', 1, size(ScalpChansLocs,2)),' ')';
% % 
% % t = num2cell(1:size(ScalpChansLocs,2));
% % [ScalpChansLocs.type] = str{1:size(t,2)};
% % 
% % rad = num2cell(ones(size(ScalpChansLocs,2),1)+84);
% % [ScalpChansLocs.sph_radius] = rad{:};
% % [ScalpChansLocs.radius] = ScalpChansLocs.r;
% % [ScalpChansLocs.sph_theta_besa] = ScalpChansLocs.sph_theta;
% % [ScalpChansLocs.sph_phi_besa] = ScalpChansLocs.sph_phi;
% % ScalpChansLocs = orderfields(ScalpChansLocs,[12:14 7:8 2 15 4:6 16:17 1 3 9:11]);
% % ScalpChansLocs = rmfield(ScalpChansLocs,{'label','phi','r','x','y'});



end