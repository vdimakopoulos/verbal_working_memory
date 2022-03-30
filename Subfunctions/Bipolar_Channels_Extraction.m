function Bipolar_Matrix_Sternberg = Bipolar_Channels_Extraction(Path)

%%
%Electrode Files in the Path
Files=dir([Path,'/*.ncs']);
File_Cell=struct2cell(Files);
for i=1:length(Files)
   Elec{i}=File_Cell{1,i};
   %Remove the first '_'
   [Labels{i},remain{i}]=strtok(Elec{i},'_');
   Elec_label{i}=strcat(Labels{i},remain{i});
end

%% To be changed
Bipolar_Matrix_Sternberg=Elec_label;





end