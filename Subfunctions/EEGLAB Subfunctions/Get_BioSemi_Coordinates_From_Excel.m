function [ ElectrodeCoordinatesTable ] = Get_BioSemi_Coordinates_From_Excel( strFilePath, strSheetName)

% Downloaded from
% http://www.biosemi.com/download/Cap_coords_all.xls
% Sheet '64-chan added channels' modified
% Channels added

ElectrodeCoordinatesTable = readtable(strFilePath,'Sheet',strSheetName);
ElectrodeCoordinatesTable.Properties.VariableNames{strcmpi(ElectrodeCoordinatesTable.Properties.VariableNames,'x')} = 'X';
ElectrodeCoordinatesTable.Properties.VariableNames{strcmpi(ElectrodeCoordinatesTable.Properties.VariableNames,'y')} = 'Y';
ElectrodeCoordinatesTable.Properties.VariableNames{strcmpi(ElectrodeCoordinatesTable.Properties.VariableNames,'z')} = 'Z';

X = ElectrodeCoordinatesTable.X;
Y = ElectrodeCoordinatesTable.Y;
Z = ElectrodeCoordinatesTable.Z;

[x,y] = Project_Cartesian_3D_On_2D(X,Y,Z);

ElectrodeCoordinatesTable.x = x;
ElectrodeCoordinatesTable.y = y;

end