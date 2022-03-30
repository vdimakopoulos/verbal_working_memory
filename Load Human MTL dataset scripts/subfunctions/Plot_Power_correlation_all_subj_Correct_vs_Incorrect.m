%% High Workload
figure;
subplot(331)
[modelCorr{1},modelIncor{1}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{1}.theta,Fr_Maint_Hipp_Incorrect{1}.theta,'theta')

subplot(332)
[modelCorr{2},modelIncor{2}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{1}.alpha,Fr_Maint_Hipp_Incorrect{1}.alpha,'alpha')

subplot(333)
[modelCorr{3},modelIncor{3}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{1}.beta,Fr_Maint_Hipp_Incorrect{1}.beta,'beta')

subplot(334)
[modelCorr{4},modelIncor{4}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{2}.theta,Fr_Maint_Hipp_Incorrect{2}.theta,'theta')


subplot(335)
[modelCorr{5},modelIncor{5}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{2}.alpha,Fr_Maint_Hipp_Incorrect{2}.alpha,'alpha')


subplot(336)
[modelCorr{6},modelIncor{6}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{2}.beta,Fr_Maint_Hipp_Incorrect{2}.beta,'beta')


subplot(337)
[modelCorr{7},modelIncor{7}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{3}.theta,Fr_Maint_Hipp_Incorrect{3}.theta,'theta')

subplot(338)
[modelCorr{8},modelIncor{8}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{3}.alpha,Fr_Maint_Hipp_Incorrect{3}.alpha,'alpha')

subplot(339)
[modelCorr{9},modelIncor{9}]  = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct{3}.beta,Fr_Maint_Hipp_Incorrect{3}.beta,'beta')



%% Scalp 

figure;
for i=1:9
   subplot(3,3,i)
   [mdl_Scalp_Corr mdl_Scalp_Incorr] = Plot_PowerCorrelation(Fr_Maint_Hipp_Correct_scalp{i}.theta,Fr_Maint_Hipp_Incorrect_Scalp{i}.theta,'theta');
   
    
end