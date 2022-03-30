function [dataBipolar_Ret_SS] = Extract_task_period_Data(TaskPeriod,dataBipolar_SS)
nSet_Size = size(dataBipolar_SS,2);
if strcmp(TaskPeriod,'fix') %Fixation
    %% Extract 1 seconds of fixation latency
    dataBipolar_Ret_SS = {};
    fsample = dataBipolar_SS{1}.fsample;
    for iSS = 1:nSet_Size
        if ~isempty(dataBipolar_SS{iSS}.time)
            
            cfg = [];
            cfg.latency = [dataBipolar_SS{iSS}.time{1}(1) dataBipolar_SS{iSS}.time{1}(fsample)]%
            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
            
            cfg =[];
            cfg.resamplefs = 1000;
            dataBipolar_Ret_SS{iSS} = ft_resampledata(cfg,dataBipolar_Ret_SS{iSS})
        else
            dataBipolar_Ret_SS{iSS} = dataBipolar_SS{iSS};
        end
        
    end
    
elseif strcmp(TaskPeriod,'encod') % Encoding
    %% Extract 2 seconds of encoding latency
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        if ~isempty(dataBipolar_SS{iSS}.time)
            cfg.latency = [-5,-3-1/dataBipolar_SS{iSS}.fsample];
            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        else
            dataBipolar_Ret_SS{iSS} = dataBipolar_SS{iSS};
        end
    end
    
elseif strcmp(TaskPeriod,'maint') %Maintenance
    %% Extract 2 seconds of maintenance latency
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        if ~isempty(dataBipolar_SS{iSS}.time)
            cfg.latency = [-2,-1/dataBipolar_SS{iSS}.fsample];
%             cfg.latency = [-3,-1/dataBipolar_SS{iSS}.fsample];

            dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        else
            dataBipolar_Ret_SS{iSS} = dataBipolar_SS{iSS};
        end
    end
    
else %Retrieval
    %% Extract 2 seconds of retrieval latency
    dataBipolar_Ret_SS = {};
    for iSS = 1:nSet_Size
        cfg = [];
        cfg.latency = [0,2-1/dataBipolar.fsample];
        dataBipolar_Ret_SS{iSS} = ft_selectdata(cfg,dataBipolar_SS{iSS});
        
    end
    
end

end