function [dataBip,TrialInformationTable_iEEG_clean] = loadHippocampalData_grid_patients(pID,analysis_type)

Drive_Letter ='F:\';

switch pID
    
    case 42
        %% Load data for 2 sessions
        Session1                   =      [Drive_Letter,'Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_01_Part_01'];
        Session2                   =        [Drive_Letter,'Vasileios\Task Analysis\Data\Macro Data\Macro_Data_Sessions_Patient_42_Session_02_Part_01'];
        data_s1                    =      load(Session1);
        data_s2                    =      load(Session2);
        TrialInformationTable      =      [data_s1.TrialInformationTable;data_s2.TrialInformationTable] ;
        data_s1                    =      data_s1.dataMacro;
        data_s2                    =      data_s2.dataMacro;
        
        %% Merge sessions
        cfg = [];
        dataBipolar = ft_appenddata(cfg,data_s1,data_s2);
        
        %% Reref based on the average of the signals
        % [dataBipRef] = ft_preproc_rereference(dataBipolar.trial, 'all', 'avg',1)
        cfg               =   [];
        cfg.reref         =  'yes'
        cfg.refchannel    =  [3 4] %white matter contacts referencing
        cfg.refmethod     =  'avg'%'avg'
        dataBipolar       = ft_preprocessing(cfg,dataBipolar);
        
        %% Apply montage
        clear montage;
        montage.labelold = dataBipolar.label;
        num_bipolar_chans=6;
        num_reference_chans=size(montage.labelold,1);
        %prepare the montage matrix
        montage_matrix = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar = [1 2 1 3 2 3 73 74 73 75 74 75 ] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [82 85]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        
        %             [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        
        dataBip = dataBipolar_SS
        TrialInformationTable_iEEG_clean = TrialInformationTable
        
        
    case 38
        %% Load data for 6 sessions - Macro Data
        strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\38 NP\'];
        cd(strMacroDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable = []
        for i=1:length(data_ses)
            TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
        end
        %% Merge sessions
        for nSes = 1:length(data_ses)
            if nSes<=2
                cfg = [];
                cfg.channel = {'all','-mTLIL*','-mTLSL5','-mTLSL6'}
                data_ses(nSes).data = ft_selectdata(cfg,data_ses(nSes).data)
            end
            if nSes == 1
                dataBipolar = data_ses(nSes).data;
                
            else
                cfg = [];
                dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
            end
        end
        
        %% Reref to different sources
        strChannelNameList = dataBipolar.label;
        AhippChans = find(contains(strChannelNameList,'AH'));
        PhippChans = find(contains(strChannelNameList,'PH'));
        GridChans = find(contains(strChannelNameList,'T'));
        dataRerefAhipp = reref_toSeparate_Chan(dataBipolar, 3, 'all','avg');
        dataRerefGrid = reref_toSeparate_Chan(dataBipolar, 19, 'all','avg');
        %selection of channels
        cfg = [];
        cfg.channel = [AhippChans; PhippChans];
        dataRerefAhipp = ft_selectdata(cfg,dataRerefAhipp);
        cfg = [];
        cfg.channel = GridChans;
        dataRerefGrid = ft_selectdata(cfg,dataRerefGrid);
        
        
        % append them again
        cfg = [];
        dataBipolar = ft_appenddata(cfg,dataRerefAhipp,dataRerefGrid);
        %% Apply montage %%
        clear montage;
        montage.labelold        = dataBipolar.label;
        num_bipolar_chans       = 9;
        num_reference_chans     = size(montage.labelold,1);
        
        %prepare the montage matrix
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra             = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [37:45]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        macro_data = dataBipolar;
        
        
        %             [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        dataBip = dataBipolar_SS;
        TrialInformationTable_iEEG_clean = TrialInformationTable;
    case 37
        %% Load data for 6 sessions - Macro Data
        strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\37 PN\'];
        cd (strMacroDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable = []
        for i=1:length(data_ses)
            TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
        end
        %% Merge sessions
        dataBipolar = data_ses(1).data;
        for nSes = 2:length(data_ses)
            cfg = [];
            dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
            
        end
        dataBipolar.label = strrep(dataBipolar.label,'m','');
        
        %% Reref to different sources
        strChannelNameList = dataBipolar.label;
        AhippChans = find(contains(strChannelNameList,'AH'));
        PhippChans = find(contains(strChannelNameList,'PH'));
        LLChans = find(contains(strChannelNameList,'LL'));
        
        GridChans = find(contains(strChannelNameList,'PL'));
        dataRerefAhipp = reref_toSeparate_Chan(dataBipolar, 4, 'all','avg');
        dataRerefGrid = reref_toSeparate_Chan(dataBipolar, 11, 'all','avg');
        %selection of channels
        cfg = [];
        cfg.channel = [AhippChans;LLChans; PhippChans];
        dataRerefAhipp = ft_selectdata(cfg,dataRerefAhipp);
        cfg = [];
        cfg.channel = GridChans;
        dataRerefGrid = ft_selectdata(cfg,dataRerefGrid);
        
        
        % append them again
        cfg = [];
        dataBipolar = ft_appenddata(cfg,dataRerefAhipp,dataRerefGrid);
        
        
        
        %% Apply montage %%
        clear montage;
        montage.labelold        = dataBipolar.label;
        num_bipolar_chans       = 15;
        num_reference_chans     = size(montage.labelold,1);
        
        %prepare the montage matrix
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 33 34 33 35 34 35] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra             = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [47 58]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        macro_data = dataBipolar;
        
        
        %             [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        dataBip = dataBipolar_SS;
        TrialInformationTable_iEEG_clean = TrialInformationTable;
    case 40
        %% Load data for 8 sessions - Macro Data
        strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\40 DG\'];
        cd (strMacroDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable = []
        for i=1:length(data_ses)
            TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
        end
        %% Merge sessions
        dataBipolar = data_ses(1).data;
        for nSes = 2:length(data_ses)
            cfg = [];
            dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
        end
        %% Rereferencing
        cfg               =   [];
        cfg.reref         =  'yes'
        cfg.refchannel    =  [61 62]%[21 22][26 27][27 28][51 52] %white matter contacts referencing
        cfg.refmethod     =  'avg'%'avg'
        dataBipolar       = ft_preprocessing(cfg,dataBipolar);
        %% Apply montage %%
        clear montage;
        montage.labelold        = dataBipolar.label;
        num_bipolar_chans       = 24;
        num_reference_chans     = size(montage.labelold,1);
        
        %prepare the montage matrix
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 ...
            33 34 33 35 34 35 41 42 41 43 42 43 49 50 49 51 50 51 57 58 57 59 58 59] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra             = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [65:88]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        macro_data = dataBipolar;
        
        %             [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        dataBip = dataBipolar_SS;
        TrialInformationTable_iEEG_clean = TrialInformationTable;
    case 44
        strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\44 MJ\'];
        cd (strMacroDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable = []
        for i=1:length(data_ses)
            TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
        end
        
        %% Merge sessions
        dataBipolar = data_ses(1).data;
        for nSes = 2:length(data_ses)
            cfg = [];
            dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
            
        end
        
        %% Rereferencing
        cfg               =   [];
        cfg.reref         =  'yes'
        cfg.refchannel    =  [3] %white matter contacts referencing
        cfg.refmethod     =  'avg'%'avg'
        dataBipolar       = ft_preprocessing(cfg,dataBipolar);
        %% Apply montage %%
        clear montage;
        montage.labelold        = dataBipolar.label;
        num_bipolar_chans       = 24;
        num_reference_chans     = size(montage.labelold,1);
        
        %prepare the montage matrix
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar     = [1 2 1 3 2 3 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 33 34 33 35 34 35 41 42 41 43 42 43 ...
            49 50 49 51 50 51 57 58 57 59 58 59] %Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra             = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [66:1:88]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        macro_data = dataBipolar;
        
        %             [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        dataBip = dataBipolar_SS;
        TrialInformationTable_iEEG_clean = TrialInformationTable;
        
    case 45
        
        %% Load data for 6 sessions - Macro Data
        strMacroDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Macro Data\45 SS\'];
        cd (strMacroDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable = []
        for i=1:length(data_ses)
            TrialInformationTable = [TrialInformationTable;data_ses(i).TrialInformationTable];
        end
        %% Merge sessions
        dataBipolar = data_ses(1).data;
        for nSes = 2:length(data_ses)
            cfg = [];
            dataBipolar = ft_appenddata(cfg,dataBipolar,data_ses(nSes).data);
            
        end
        %% Rereferencing
        cfg               =   [];
        cfg.reref         =  'yes'
        cfg.refchannel    =  [64 65] %white matter contacts referencing
        cfg.refmethod     =  'avg'%'avg'
        dataBipolar       = ft_preprocessing(cfg,dataBipolar);
        
        %% Apply montage %%
        clear montage;
        montage.labelold        = dataBipolar.label;
        num_bipolar_chans       = 21;
        num_reference_chans     = size(montage.labelold,1);
        
        %prepare the montage matrix
        montage_matrix          = eye(num_bipolar_chans+num_reference_chans,num_reference_chans);
        chans_to_be_bipolar     = [1 4 2 4 3 4 9 10 9 11 10 11 17 18 17 19 18 19 25 26 25 27 26 27 ...
            45 46 45 47 46 47 53 54 53 55 54 55 61 62 61 63 62 63]%Bipolar Channels would be: the pair combinations like Chan1-2 or Chan1-3 etc.
        sign = 1;
        for i = 1:size(chans_to_be_bipolar,2)
            montage_matrix(num_reference_chans+round(i/2),chans_to_be_bipolar(i)) = sign;
            sign = sign*(-1);
        end
        %Append the bipolar channel labels to the reference channels labels
        for i = 1:2:size(chans_to_be_bipolar,2)
            montage.labelnew(round(i/2)) = strcat(dataBipolar.label(chans_to_be_bipolar(round(i))),'-',dataBipolar.label(chans_to_be_bipolar(round(i+1))))
        end
        montage.labelnew        = {dataBipolar.label{:},montage.labelnew{:}};
        montage.tra             = montage_matrix;
        % dataBipolar = ft_apply_montage(dataBipolar,montage);
        
        cfg = [];
        % cfg.reref = 'yes'
        cfg.refmethod  = 'bipolar'
        cfg.refchannel = [69 89]
        cfg.montage = montage;
        dataBipolar = ft_preprocessing(cfg,dataBipolar);
        
        Bip_chans = (find(contains(dataBipolar.label, '-')==1));
        macro_data = dataBipolar;
        
        %              [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable] = Get_Only_Correct_Trials_FieldTrip(dataBipolar,TrialInformationTable);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable] = Get_Only_Incorrect_Trials_FieldTrip(dataBipolar, TrialInformationTable )
        end
        dataBip = dataBipolar_SS;
        TrialInformationTable_iEEG_clean = TrialInformationTable;
        
end
end