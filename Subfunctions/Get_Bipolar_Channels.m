montage.labelold        = 	dataBipolar.label;
num_bipolar_chans 		= 	length(dataBipolar.label)-1;
num_reference_chans     =  size(montage.labelold,1);
montage_matrix          =  zeros(num_bipolar_chans,num_reference_chans);%eye(num_bipolar_chans+num_reference_chans,num_reference_chans)
chans_to_be_bipolar     = repelem([1:num_reference_chans],2);
chans_to_be_bipolar     = chans_to_be_bipolar(2:end-1);
sign = 1;

% chans_to_be_bipolar = [[1:num_reference_chans-1];[2:num_reference_chans]]';
% for ii = 1:size(chans_to_be_bipolar,1)
%     montage_matrix(ii,chans_to_be_bipolar(ii,1)) = 1;
%     montage_matrix(ii,chans_to_be_bipolar(ii,2)) = -1;
% end

for i = 1:size(chans_to_be_bipolar,2)
    montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
    sign = sign*(-1);
end

%%
dataBipolar.label(chans_to_be_bipolar)


%%
%Append the bipolar channel labels to the reference channels labels
for i = 1:2:size(chans_to_be_bipolar,2)
    montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
end
montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
montage.tra             = montage_matrix;

ind_BipolarChan_start = length(dataBipolar.label)+1;
ind_BipolarChan_end = length(dataBipolar.label)+ num_bipolar_chans;
cfg = [];
cfg.refmethod  = 'bipolar'
cfg.refchannel = [ind_BipolarChan_start ind_BipolarChan_end]
cfg.montage = montage;
dataBipolar_referenced = ft_preprocessing(cfg,dataBipolar);




