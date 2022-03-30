function [ MontageMatrix ] = Get_Bipolar_Montage_Matrix( nPatient )

Temp = [1,-1,zeros(1,6)];
MontageMatrix1 = [];
for ii = 1:7
    MontageMatrix1 = [MontageMatrix1;circshift(Temp,ii-1,2)];
end
MontageMatrix1 = [MontageMatrix1,zeros(7,56)];

MontageMatrix_8_8 = [];
for ii = 1:8
    MontageMatrix_8_8 = [MontageMatrix_8_8;circshift(MontageMatrix1,(ii-1)*8,2)];
end

switch nPatient
    case 1 % Special case
        MontageMatrix = Get_Patient_01_GR_Montage_Matrix();
    case 2
        MontageMatrix = MontageMatrix_8_8(1:49,1:56);
    case 4
        MontageMatrix = MontageMatrix_8_8(1:35,1:40);
    case 10
        MontageMatrix = MontageMatrix_8_8(1:35,1:40);
    case 12
        nChanOrder = [1:11:56 62:64 2:11 13:22 24:33 35:44 46:55 57:61 65 67:74 66];
        chans_to_be_bipolar = [nChanOrder(1:end-1);[nChanOrder(2:end)]]';
        chans_to_be_bipolar = chans_to_be_bipolar(setdiff([1:length(chans_to_be_bipolar)],[8:8:64]),:);
        MontageMatrix = [];
        for ii = 1:size(chans_to_be_bipolar,1)
            MontageMatrix(ii,chans_to_be_bipolar(ii,1)) = 1;
            MontageMatrix(ii,chans_to_be_bipolar(ii,2)) = -1;
        end
    case 13
        MontageMatrix = MontageMatrix_8_8;
    case 16
        MontageMatrix = MontageMatrix_8_8;
    case 19
        MontageMatrix = MontageMatrix_8_8(1:35,1:40);
    case 20
        nChanOrder = [1:11:23 27:32 2:11 13:22 24:26 33 41:48 34:40];
        chans_to_be_bipolar = [nChanOrder(1:end-1);[nChanOrder(2:end)]]';
        chans_to_be_bipolar = chans_to_be_bipolar(setdiff([1:length(chans_to_be_bipolar)],[8:8:40]),:);
        MontageMatrix = [];
        for ii = 1:size(chans_to_be_bipolar,1)
            MontageMatrix(ii,chans_to_be_bipolar(ii,1)) = 1;
            MontageMatrix(ii,chans_to_be_bipolar(ii,2)) = -1;
        end
    case 22
        MontageMatrix = MontageMatrix_8_8;
    case 23
        MontageMatrix = MontageMatrix_8_8;
    case 26
        nChanOrder = [1:8:9 10:2*8 2:8 17 19:26 18 27:8:35 36:42 28:34]
        chans_to_be_bipolar = [nChanOrder(1:end-1);nChanOrder(2:end)]';
        chans_to_be_bipolar = chans_to_be_bipolar(setdiff([1:length(chans_to_be_bipolar)],[8:8:16 26:8:34]),:);
        MontageMatrix  = [];
         for ii = 1:size(chans_to_be_bipolar,1)
            MontageMatrix(ii,chans_to_be_bipolar(ii,1)) = 1;
            MontageMatrix(ii,chans_to_be_bipolar(ii,2)) = -1;
        end
    case 28
        MontageMatrix = MontageMatrix_8_8(1:42,1:48);
    case 29
        MontageMatrix = MontageMatrix_8_8;
        MontageMatrix = MontageMatrix(1:55,[1,3:63]);
        MontageMatrix(1:end-1,2:end) = MontageMatrix(2:end,2:end);
        MontageMatrix = MontageMatrix(1:end-1,:);
    case 30
        MontageMatrix = MontageMatrix_8_8;
    case 32
        MontageMatrix = MontageMatrix_8_8;
    case 33
        MontageMatrix = MontageMatrix_8_8(1:28,1:32);
    case 34
        MontageMatrix_8_8 = MontageMatrix_8_8(1:42,1:48);
        MontageMatrix_8_8(6:7,:) = [];
        MontageMatrix_8_8(:,7:8) = [];
        MontageMatrix_8_8(11:12,:) = [];
        MontageMatrix_8_8(:,13:14) = [];
        MontageMatrix = MontageMatrix_8_8;
    case 35
        MontageMatrix = MontageMatrix_8_8;
    case 36
        MontageMatrix = MontageMatrix_8_8(1:14,1:16);
    case 37
        MontageMatrix = MontageMatrix_8_8(1:40,1:46);
    case 38
        MontageMatrix = MontageMatrix_8_8(1:30,1:36);
        MontageMatrix(25:end,:) = circshift(MontageMatrix(25:end,:),[0,1]);
        MontageMatrix(28:end,:) = circshift(MontageMatrix(28:end,:),[0,1]);
        MontageMatrix(29:end,:) = circshift(MontageMatrix(29:end,:),[0,-1]);
    case 40
        MontageMatrix = MontageMatrix_8_8;
    case 41
        MontageMatrix = MontageMatrix_8_8(1:35,:);
    case 42
        chans_to_be_bipolar = [[1:80-1];[2:80]]';
        chans_to_be_bipolar = chans_to_be_bipolar(setdiff([1:length(chans_to_be_bipolar)],[8:8:72]),:);
        MontageMatrix  = [];
        for ii = 1:size(chans_to_be_bipolar,1)
            MontageMatrix(ii,chans_to_be_bipolar(ii,1)) = 1;
            MontageMatrix(ii,chans_to_be_bipolar(ii,2)) = -1;
        end
    case 44
        MontageMatrix = MontageMatrix_8_8;
end

end

%% For patient 1 GR
function [ MontageMatrix ] = Get_Patient_01_GR_Montage_Matrix( )

MontageMatrix = zeros(56,64);

% The order of channels is different in data.label
OrderedChannelNumbers = zeros(64,1);
OrderedChannelNumbers(8*0+1:8*1) = [1;12;23;34;45;56;62;63]; % mAR
OrderedChannelNumbers(8*1+1:8*2) = [64;(2:8)']; % mAL
OrderedChannelNumbers(8*2+1:8*3) = [(9:11)';(13:17)']; % mAHR
OrderedChannelNumbers(8*3+1:8*4) = [(18:22)';(24:26)']; % mAHL
OrderedChannelNumbers(8*4+1:8*5) = [(27:33)';35]; % mECR
OrderedChannelNumbers(8*5+1:8*6) = (36:43)'; % mECL
OrderedChannelNumbers(8*6+1:8*7) = [44;(46:52)']; % mPHR
OrderedChannelNumbers(8*7+1:8*8) = [(53:55)';(57:61)']; % mPHL

nNumberOfGroups = 8;
nNumberOfElectrodesPerGroup = 8;
nNumberOfBipolarChannelsPerGroup = nNumberOfElectrodesPerGroup-1;

for iGroup = 1:nNumberOfGroups
    for iElectrode = 1:nNumberOfBipolarChannelsPerGroup
        nChannel1 = OrderedChannelNumbers((iGroup-1)*nNumberOfElectrodesPerGroup+iElectrode);
        nChannel2 = OrderedChannelNumbers((iGroup-1)*nNumberOfElectrodesPerGroup+iElectrode+1);
        MontageMatrix((iGroup-1)*nNumberOfBipolarChannelsPerGroup+iElectrode,nChannel1) = 1;
        MontageMatrix((iGroup-1)*nNumberOfBipolarChannelsPerGroup+iElectrode,nChannel2) = -1;
    end
end

end

