function [ data ] = Get_Added_Channels_Dataset_From_Channel_Names_FieldTrip( data, strAddedChannelLabelList )

for iChannel = 1:length(strAddedChannelLabelList)
    strAddedChannelLabel = strAddedChannelLabelList{iChannel};
    if(~isempty(strfind(strAddedChannelLabel,'2-4')))
        nChannelCombination = 3;
        strElectrodeName = strrep(strAddedChannelLabel,'2-4','');
    elseif(~isempty(strfind(strAddedChannelLabel,'1-3')))
        nChannelCombination = 1;
        strElectrodeName = strrep(strAddedChannelLabel,'1-3','');
    elseif(~isempty(strfind(strAddedChannelLabel,'1-4')))
        nChannelCombination = 2;
        strElectrodeName = strrep(strAddedChannelLabel,'1-4','');
    end
    
    switch nChannelCombination
        case 1 % 1-3
            nChannel_1 = find(strcmpi(data.label,[strElectrodeName,'1-2']));
            nChannel_2 = find(strcmpi(data.label,[strElectrodeName,'2-3']));
            
            for nTrial = 1:length(data.trial)
                data.trial{nTrial} = [data.trial{nTrial};...
                    [data.trial{nTrial}(nChannel_1,:)+data.trial{nTrial}(nChannel_2,:)]];
            end
            data.label = [data.label,[strElectrodeName,'1-3']];
            
        case 2 % 1-4
            nChannel_1 = find(strcmpi(data.label,[strElectrodeName,'1-2']));
            nChannel_2 = find(strcmpi(data.label,[strElectrodeName,'2-3']));
            nChannel_3 = find(strcmpi(data.label,[strElectrodeName,'3-4']));
            
            for nTrial = 1:length(data.trial)
                data.trial{nTrial} = [data.trial{nTrial};...
                    [data.trial{nTrial}(nChannel_1,:)+data.trial{nTrial}(nChannel_2,:)+data.trial{nTrial}(nChannel_3,:)]];
            end
            data.label = [data.label,[strElectrodeName,'1-4']];
        case 3 % 2-4
            nChannel_1 = find(strcmpi(data.label,[strElectrodeName,'2-3']));
            nChannel_2 = find(strcmpi(data.label,[strElectrodeName,'3-4']));
            
            for nTrial = 1:length(data.trial)
                data.trial{nTrial} = [data.trial{nTrial};...
                    [data.trial{nTrial}(nChannel_1,:)+data.trial{nTrial}(nChannel_2,:)]];
            end
            data.label = [data.label,[strElectrodeName,'2-4']];
    end
end

end