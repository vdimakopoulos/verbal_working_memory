function [Power_Maint_Hipp Fr_Fix_Hipp] = getPower_trialWise(dataBipolar_SS,dataBipolar_baseline, iSS,grangerChannelPair,foi_1,foi_2,band)


nTrials = length(dataBipolar_SS{iSS}.trial);
chanPairs = size(grangerChannelPair);
for jCh = 1:chanPairs(1)
    for iCh = 1:chanPairs(2)
        if iCh == 1
            fr = foi_1;
        else
            fr = foi_2;
        end
        for iTrial =1:nTrials
            cfg                         = [];
            cfg.method                  = 'mtmfft';
            cfg.taper                   = 'dpss';
            cfg.output                  = 'pow';
            cfg.foi                     = fr;
            cfg.tapsmofrq               = 2;
            cfg.trials                  = iTrial;
            cfg.channel                 = grangerChannelPair{jCh,iCh}
            Fr_Maint_Hipp{iCh,iTrial}   = ft_freqanalysis(cfg,dataBipolar_SS{iSS});
            Fr_Fix_Hipp{iCh,iTrial}     = ft_freqanalysis(cfg,dataBipolar_baseline{iSS});
            FixPower = 10*log10(Fr_Fix_Hipp{iCh,iTrial}.powspctrm(1:length(band)));
            MaintPower = 10*log10(Fr_Maint_Hipp{iCh,iTrial}.powspctrm(1:length(band)));
            Fr_Maint_Hipp{iCh,iTrial}.meanBaselinedPower = (MaintPower-FixPower)./FixPower;
            Fr_Maint_Hipp{iCh,iTrial}.meanBaselinedPower = mean(Fr_Maint_Hipp{iCh,iTrial}.meanBaselinedPower);
        end 
    end
    Power_Maint_Hipp{jCh} = Fr_Maint_Hipp;
    
end