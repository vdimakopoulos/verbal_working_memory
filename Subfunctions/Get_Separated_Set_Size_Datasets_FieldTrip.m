function [ data_SS, TrialInformationTable_SS, nTrialList_TT_SS ] = Get_Separated_Set_Size_Datasets_FieldTrip( data_All, TrialInformationTable )

data_SS = cell(1,5);
TrialInformationTable_SS = cell(1,5);
nTrialList_TT_SS = cell(1,4);
for iSS = 1:5
    if(ismember(iSS,1:3))
        nSetSize = 2*iSS+2;
    elseif(iSS==4)
        nSetSize = [6,8];
    elseif(iSS==5)
        nSetSize = [4,6,8];
    end
    [data_SS{iSS},TrialInformationTable_SS{iSS},nTrialList_TT_SS{iSS}] = Get_Selected_Set_Size_Data_FieldTrip(data_All,nSetSize,TrialInformationTable);
end

end