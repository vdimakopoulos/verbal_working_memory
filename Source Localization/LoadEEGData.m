function [dataScalp,TrialInformationTable_Scalp, data_Task_Scalp] = LoadEEGData(pID,Drive_Letter,analysis_type)

switch pID
    
    case 42
        %% Load scalp data for 2 sessions
        Sclp_Session1                    =      [Drive_Letter,'Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_01_Part_01'];
        Sclp_Session2                    =      [Drive_Letter,'Vasileios\Task Analysis\Data\Scalp_Data\Scalp_Data_Sessions_Patient_42_Session_02_Part_01'];
        data_s1_scalp                    =      load(Sclp_Session1);
        data_s2_scalp                    =      load(Sclp_Session2);
        TrialInformationTable_Scalp      =      [data_s1_scalp.TrialInformationTable;data_s2_scalp.TrialInformationTable] ;
        data_scalp_s1                    =      data_s1_scalp.dataScalp;
        data_scalp_s2                    =      data_s2_scalp.dataScalp;
        
        %% merge
        cfg = [];
        dataScalp = ft_appenddata(cfg,data_scalp_s1,data_scalp_s2);
        
        %% Re-reference
        cfg               =   [];
        cfg.reref         =  'yes'
        cfg.refchannel    =  [1 2] %white matter contacts referencing
        cfg.refmethod     =  'avg'%'avg'
        dataScalp       = ft_preprocessing(cfg,dataScalp);
        
        %% reject channels
        cfg = [];
        cfg.channel = {'all', '-_Submp','-_Submm'};
        dataScalp = ft_preprocessing(cfg,dataScalp);
        
        %% downsample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataScalp);
        dataScalp.label = strrep(dataScalp.label,'_','');
        %         data_Task_Scalp = [];
        %         IC_components = [];
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
    case 37
        %% Load data for 6 sessions - Scalp Data
        strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\37 PN\'];
        cd (strScalpDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_Scalp_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable_Scalp = [];
        for i=1:length(data_Scalp_ses)
            TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
        end
        
        %% merge sessions
        
        dataBipolarScalp = data_Scalp_ses(1).data;
        for nSes = 2:length(data_Scalp_ses)
            if nSes == 6 %% remove the two redundant channels which are present only on this session
                cfg =[];
                cfg.channel = {'all','-Subm1','-Subm2'};
                data_Scalp_ses(nSes).data = ft_selectdata(cfg, data_Scalp_ses(nSes).data);
            end
            cfg = [];
            dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
            
            
        end
        
        %% re-ref EEG
        cfg               = []
        cfg.channel       = {'all'}
        cfg.reref         =  'yes'
        cfg.refchannel    =   [1, 2]
        cfg.refmethod     =  'avg'
        dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp);
        
        %% resample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataBipolarScalp);
        
        %% correct trials
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
        
    case 38
        
        %% Load data for 6 sessions - Scalp Data
        strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\38 NP\'];
        cd (strScalpDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_Scalp_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable_Scalp = [];
        for i=1:length(data_Scalp_ses)
            TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
        end
        
        %% append
        
        dataBipolarScalp = data_Scalp_ses(1).data;
        for nSes = 2:length(data_Scalp_ses)
            cfg = [];
            dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
            
        end
        %         cd(strPaths.Project)
        %% rereference
        cfg               = []
        cfg.channel       = {'all'}
        cfg.reref         =  'yes'
        cfg.refchannel    =   [1, 2]
        cfg.refmethod     =  'avg'
        dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp);
        %% resample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataBipolarScalp);
        
        %% correct trials
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
        
    case 40
        %% Load data for 8 sessions - Scalp Data
        strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\40 DG\'];
        cd (strScalpDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_Scalp_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable_Scalp = [];
        for i=1:length(data_Scalp_ses)
            TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
        end
        
        dataBipolarScalp = data_Scalp_ses(1).data;
        for nSes = 2:length(data_Scalp_ses)
            cfg = [];
            dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
        end
        
        %% re-ref EEG
        cfg               = []
        cfg.channel       = {'all'}
        cfg.reref         =  'yes'
        cfg.refchannel    =   [1, 2]
        cfg.refmethod     =  'avg'
        dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp);
        %% resample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataBipolarScalp);
        
        %% correct trials
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
        
    case 44
        
        %% Load data for 4 sessions - Scalp Data
        strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\44 MJ\'];
        cd (strScalpDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_Scalp_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable_Scalp = [];
        for i=1:length(data_Scalp_ses)
            TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
        end
        %% Merge sessions
        dataBipolarScalp = data_Scalp_ses(1).data;
        for nSes = 2:length(data_Scalp_ses)
            cfg = [];
            dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
            
        end
        
        %% re-ref EEG
        cfg               = []
        cfg.channel       = {'all'}
        cfg.reref         =  'yes'
        cfg.refchannel    =   [1, 2]
        cfg.refmethod     =  'avg'
        dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp);
        
        %% resample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataBipolarScalp);
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
    case 45
        %% Load data for 6 sessions - Scalp Data
        strScalpDataDir = [Drive_Letter,'Vasileios\Task Analysis\Data\Sternberg Task\Sessions\Scalp Data\45 SS\'];
        cd (strScalpDataDir);
        files = dir('*.mat'); %try to look for all the .mat files under the folder
        for i=1:length(files)
            data_Scalp_ses(i) = load(files(i).name); %load the files
        end
        
        TrialInformationTable_Scalp = [];
        for i=1:length(data_Scalp_ses)
            TrialInformationTable_Scalp = [TrialInformationTable_Scalp;data_Scalp_ses(i).TrialInformationTable];
        end
        dataBipolarScalp = data_Scalp_ses(1).data;
        for nSes = 2:length(data_Scalp_ses)
            cfg = [];
            dataBipolarScalp = ft_appenddata(cfg,dataBipolarScalp,data_Scalp_ses(nSes).data);
            
        end
        %% re-ref EEG
        cfg               = []
        cfg.channel       = {'all'}
        cfg.reref         =  'yes'
        cfg.refchannel    =   [1, 2]
        cfg.refmethod     =  'avg'
        dataBipolarScalp       = ft_preprocessing(cfg,dataBipolarScalp);
        
         %% resample
        cfg = [];
        cfg.resamplefs = 80;
        dataScalp = ft_resampledata(cfg,dataBipolarScalp);
        %% correct/incorrect trials
        switch analysis_type
            case 'correct trials'
                [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
            case 'incorrect trials'
                [dataBipolar_SS, TrialInformationTable_Scalp] = Get_Only_Incorrect_Trials_FieldTrip(dataScalp, TrialInformationTable_Scalp )
        end
        
%         %% correct trials
%         [dataBipolar_SS,TrialInformationTable_Scalp] = Get_Only_Correct_Trials_FieldTrip(dataScalp,TrialInformationTable_Scalp);
        %% set size
        Set_Sizes = {4;6;8;[6 8]; [4 6 8]};
        %     TrialInformationTable_Scalp.SetSize = TrialInformationTable_Scalp.Setsize;
        [dataBipolar_SS,TrialInformationTable_Scalp,nTrialList_TT_SS] = Get_Separated_Set_Size_Datasets_FieldTrip(dataBipolar_SS,TrialInformationTable_Scalp);
        %% latencies for every task period
        %% Select the latencies for every task period for each task period
        [data_Task_Scalp{1}] = Extract_task_period_Data('fix',dataBipolar_SS); % fixation
        [data_Task_Scalp{2}] = Extract_task_period_Data('encod',dataBipolar_SS); % encoding
        [data_Task_Scalp{3}] = Extract_task_period_Data('maint',dataBipolar_SS); % maintenance
    otherwise
        [dataScalp, data_Task_Scalp,TrialInformationTable_Scalp] = LoadHumanMTL_EEG(pID,analysis_type)
        TrialInformationTable_Scalp = TrialInformationTable_Scalp.Scalp;
end
