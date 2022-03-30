function [data_arr, data_arr_enc] = ft_prepare_GrangerData_permutation_test(strPath,flag_incorrect)
load(strPath)
if ~flag_incorrect
%% maintenance
data_arr{1}.powspctrm = Granger_Spectra_Fig.spectra_recalculated{1}.Maint.grangerspctrm(1,:);
data_arr{2}.powspctrm = Granger_Spectra_Fig.spectra_recalculated{1}.Maint.grangerspctrm(2,:);

for i = 2:15
    data_arr{1}.powspctrm = [data_arr{1}.powspctrm; Granger_Spectra_Fig.spectra_recalculated{i}.Maint.grangerspctrm(1,:)];
    data_arr{1}.label{1} = 'Hipp - cortex maint';
    data_arr{1}.dimord = 'subj_freq';
    data_arr{1}.freq = Granger_Spectra_Fig.freqAxis; 
    data_arr{2}.powspctrm = [data_arr{2}.powspctrm;Granger_Spectra_Fig.spectra_recalculated{i}.Maint.grangerspctrm(2,:)];
    data_arr{2}.label{1} = 'Cortex hipp maint';
    data_arr{2}.freq = Granger_Spectra_Fig.freqAxis;
    data_arr{2}.dimord = 'subj_freq';

end

%% encoding
data_arr_enc{1}.powspctrm = Granger_Spectra_Fig.spectra_recalculated{1}.Enc.grangerspctrm(1,:);
data_arr_enc{2}.powspctrm = Granger_Spectra_Fig.spectra_recalculated{1}.Enc.grangerspctrm(2,:);

for i = 2:15
    data_arr_enc{1}.powspctrm = [data_arr_enc{1}.powspctrm; Granger_Spectra_Fig.spectra_recalculated{i}.Enc.grangerspctrm(1,:)];
    data_arr_enc{1}.label{1} = 'Hipp - cortex enc';
    data_arr_enc{1}.dimord = 'subj_freq';
    data_arr_enc{1}.freq = Granger_Spectra_Fig.freqAxis; 
    data_arr_enc{2}.powspctrm = [data_arr_enc{2}.powspctrm;Granger_Spectra_Fig.spectra_recalculated{i}.Enc.grangerspctrm(2,:)];
    data_arr_enc{2}.label{1} = 'Cortex hipp enc';
    data_arr_enc{2}.freq = Granger_Spectra_Fig.freqAxis;
    data_arr_enc{2}.dimord = 'subj_freq';

end
else
%% incorrect trials maint

strPath_incor = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Granger for human MTL dataset\Granger_spectra_incorrectTrials.mat';
load(strPath_incor)
data_arr_inc{1}.powspctrm = Granger_spectra_incorrectTrials{1}.Maint.grangerspctrm(1,4:100);
data_arr_inc{2}.powspctrm = Granger_spectra_incorrectTrials{1}.Maint.grangerspctrm(2,4:100);

for i = 2:15
    if size(Granger_spectra_incorrectTrials{i}.Maint.grangerspctrm(1,:),2) > 97
        test = Granger_spectra_incorrectTrials{i}.Maint.grangerspctrm(1,4:100);
    end
    
      if size(Granger_spectra_incorrectTrials{i}.Maint.grangerspctrm(2,:),2) > 97
        test2 = Granger_spectra_incorrectTrials{i}.Maint.grangerspctrm(2,4:100);
    end
    data_arr_inc{1}.powspctrm = [data_arr_inc{1}.powspctrm; test];
    data_arr_inc{1}.label{1} = 'Hipp - cortex maint';
    data_arr_inc{1}.dimord = 'subj_freq';
    data_arr_inc{1}.freq = [4:100]; 
    data_arr_inc{2}.powspctrm = [data_arr_inc{2}.powspctrm;test2];
    data_arr_inc{2}.label{1} = 'Cortex hipp maint';
    data_arr_inc{2}.freq = [4:100];%Granger_Spectra_Fig.freqAxis;
    data_arr_inc{2}.dimord = 'subj_freq';

end


%% Incorrect Encoding

data_arr_inc_enc{1}.powspctrm = Granger_spectra_incorrectTrials{1}.Enc.grangerspctrm(1,4:100);
data_arr_inc_enc{2}.powspctrm = Granger_spectra_incorrectTrials{1}.Enc.grangerspctrm(2,4:100);

for i = 2:15
     if size(Granger_spectra_incorrectTrials{i}.Enc.grangerspctrm(1,:),2) > 97
        test = Granger_spectra_incorrectTrials{i}.Enc.grangerspctrm(1,4:100);
    end
    
      if size(Granger_spectra_incorrectTrials{i}.Enc.grangerspctrm(2,:),2) > 97
        test2 = Granger_spectra_incorrectTrials{i}.Enc.grangerspctrm(2,4:100);
    end
    data_arr_inc_enc{1}.powspctrm = [data_arr_inc_enc{1}.powspctrm; test];
    data_arr_inc_enc{1}.label{1} = 'Hipp - cortex enc';
    data_arr_inc_enc{1}.dimord = 'subj_freq';
    data_arr_inc_enc{1}.freq = [4:100];%Granger_Spectra_Fig.freqAxis; 
    data_arr_inc_enc{2}.powspctrm = [data_arr_inc_enc{2}.powspctrm;test2];
    data_arr_inc_enc{2}.label{1} = 'Cortex hipp enc';
    data_arr_inc_enc{2}.freq = [4:100];%Granger_Spectra_Fig.freqAxis;
    data_arr_inc_enc{2}.dimord = 'subj_freq';

end

data_arr = data_arr_inc;
data_arr_enc = data_arr_inc_enc;
end
