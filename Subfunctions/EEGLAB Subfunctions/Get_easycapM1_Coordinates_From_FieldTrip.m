function [ lay ] = Get_easycapM1_Coordinates_From_FieldTrip( strFieldTripToolboxPath )

%% Channel locations template / FieldTrip
lay = load([strFieldTripToolboxPath,'/template/layout/easycapM1.mat']);
lay = lay.lay;

%% Replace channel names Add channels to layout
strExisting = {'T7','T8','P7','P8'};
strAdded = {'T3','T4','T5','T6'};
for nChannel = 1:length(strExisting)
    ind = strcmpi(lay.label,strExisting{nChannel});
    lay.label = [lay.label;strAdded{nChannel}];
    lay.pos = [lay.pos;lay.pos(ind,:)];
    lay.width = [lay.width;lay.width(ind)];
    lay.height = [lay.height;lay.height(ind)];
end

end