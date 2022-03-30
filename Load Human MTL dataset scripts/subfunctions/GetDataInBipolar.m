function [montage,chans_to_be_bipolar,dataBip,Bip_chans] = GetDataInBipolar(data,nSubj,flagGC_hipp_iEEG_Ctx)
clear montage;
montage.labelold        = data.label;
num_bipolar_chans       = round(length(montage.labelold)/8*2);
num_reference_chans     = size(montage.labelold,2);

montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);

switch nSubj
    case 1
        chans_to_be_bipolar = [1 3 2 3 9 11 10 11 17 19 18 19 25 27 26 27 33 34 33 35 41 42 41 43];
    case 2
        chans_to_be_bipolar = [1 3 2 3 9 11 10 11 17 18 17 19 25 27 26 27 33 34 33 35 41 43 42 43 49 50 49 51 57 59 58 59];
    case 3
        chans_to_be_bipolar = [1 4 2 4 9 13 10 13 17 19 18 19 25 28 26 28 33 36 34 36];
    case 4
        chans_to_be_bipolar = [1 4 2 4 9 11 10 11 17 19 18 19 25 27 26 27 33 35 34 35 41 42 41 43 49 50 50 51 57 59 58 59];
    case 5
        chans_to_be_bipolar = [1 3 2 3 9 12 10 12 17 20 18 20 25 27 26 27];
    case 6
        chans_to_be_bipolar = [1 5 2 5 9 12 10 12 17 19 18 19 25 28 26 28 33 36 34 36 41 42 41 43 49 51 50 51 57 59 58 59];
    case 7
        chans_to_be_bipolar = [1 4 2 4 9 12 10 12 17 20 18 20 25 28 26 28 33 36 34 36 41 43 42 43 49 52 50 52 57 59 58 59];
    case 8
        chans_to_be_bipolar = [1 4 2 4 9 12 10 12 17 19 18 19 25 27 26 27 33 35 34 35 41 45 42 45 49 50 49 51 57 59 58 59];
        if flagGC_hipp_iEEG_Ctx
            chans_to_be_bipolar = [chans_to_be_bipolar 5 6 6 7 13 14 14 15 21 22 22 23 29 30 30 31 37 38 38 39 45 46 46 47 53 54 54 55 61 62];
            clear montage;
            montage.labelold        = data.label;
            num_bipolar_chans       = length(chans_to_be_bipolar)/2;
            num_reference_chans     = size(montage.labelold,2);
            
            montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        end
    case 9 
        chans_to_be_bipolar = [1 4 2 4 9 11 10 11 17 20 18 20 25 26 25 27 33 35 34 35 41 44 42 44 49 52 50 52 57 60 58 60];
end

if flagGC_hipp_iEEG_Ctx
    if nSubj ~= 8
        for i = 6:8:length(data.label)
            chans_to_be_bipolar = [chans_to_be_bipolar i i+1 i+1 i+2];
        end
   
        clear montage;
        montage.labelold        = data.label;
        num_bipolar_chans       = round(length(montage.labelold)/8*2)+round(length(montage.labelold)/8*2);
        num_reference_chans     = size(montage.labelold,2);
        
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
     end
end
sign = 1;
for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end

for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(data.label(chans_to_be_bipolar(round(i))),'-',data.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew        = {data.label{:},montage.labelnew{:}};
montage.tra             = montage_matrix;

cfg = [];
% cfg.reref = 'yes'
cfg.refmethod  = 'bipolar'
cfg.refchannel = [num_reference_chans+1 num_reference_chans+num_bipolar_chans];
cfg.montage = montage;
dataBip = ft_preprocessing(cfg,data);

Bip_chans = (find(contains(dataBip.label, '-')==1));