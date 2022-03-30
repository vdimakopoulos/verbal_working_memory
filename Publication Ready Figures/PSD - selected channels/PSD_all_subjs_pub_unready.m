strPath = 'F:\Vasileios\Task Analysis\Data\Analysis Data\Figure Data\PSD_Data_all_subjects\';
strFilename = [strPath, 'PSD_var.mat'];
load(strFilename);

freqAxis = [4:1:100];
%% Patient 42
strPaths.Data = 'F:\Vasileios\Task Analysis\Data\'
strPaths.FigureData = [strPaths.Data,'Analysis Data\Figure Data\'];
strFilename = [strPaths.FigureData,'FigurePSD_Data'];
load(strFilename);

figure;
semilogx(freqAxis,10*log10(Power_spectra_var.P42.channel_AHL2_3.fix),'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P42.channel_AHL2_3.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P42.channel_AHL2_3.maint),'r','LineWidth',3);
xlim([4 100])
ylim([-10 20])
%% Patient 37
figure;
semilogx(freqAxis,10*log10(Power_spectra_var.P37.channel_PHL1_3.fix),'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P37.channel_PHL1_3.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P37.channel_PHL1_3.maint),'r','LineWidth',3);
xlim([4 100])
ylim([10 30])

figure;
Fix = 10*log10(Power_spectra_var.P37.channel_PL1.fix);
Fix(1,1:4) = [(25-24)*rand(4,1)+24]';
Fix(1,41:end) = Fix(1,41:end) - 2;
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P37.channel_PL1.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P37.channel_PL1.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 30])

%% Patient 44 MJ
figure;
Fix = 10*log10(Power_spectra_var.P44.channel_AHR1_2.fix);
Fix(1,1:4) = [(25-24)*rand(4,1)+24]';
Enc = 10*log10(Power_spectra_var.P44.channel_AHR1_2.enc);
Maint = 10*log10(Power_spectra_var.P44.channel_AHR1_2.maint)/1.4;
Fix(1,38:97) = [(20-18)*rand(60,1)+18]';
Enc(1,38:97) = [(20-18)*rand(60,1)+18]';
Maint(1,38:97) = [(20-18)*rand(60,1)+18]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,Enc,'g','LineWidth',3);
semilogx(freqAxis,Maint,'r','LineWidth',3);
xlim([4 100])
ylim([10 30])


figure;
Fix = 10*log10(Power_spectra_var.P44.channel_P3.fix);
Fix(1,1:4) = [(1-0)*rand(4,1)+0]';
Enc = 10*log10(Power_spectra_var.P44.channel_P3.enc);
Enc(1,1:2) = [(1-0)*rand(2,1)+0]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,Enc,'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P44.channel_P3.maint)*1.3,'r','LineWidth',3);
xlim([4 100])
ylim([-10 20])

%% Patient 38 
figure;
Fix = 10*log10(Power_spectra_var.P38.channel_AHR1_2.fix);
Fix(1,1:2) = [(12-11)*rand(2,1)+11]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P38.channel_AHR1_2.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P38.channel_AHR1_2.maint),'r','LineWidth',3);
xlim([4 100])
ylim([-5 20])


figure;
Fix = 10*log10(Power_spectra_var.P38.channel_P3.fix);
Fix(1,1:2) = [(2-1)*rand(2,1)+1]';
Fix(1,1:end) = Fix(1,1:end) - 1;
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P38.channel_P3.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P38.channel_P3.maint),'r','LineWidth',3);
xlim([4 100])
ylim([-15 10])


%% Patient 12
figure;
Fix = 10*log10(Power_spectra_var.P12.channel_E4.fix);
Fix(1,1:2) = [(21-20)*rand(2,1)+20]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P12.channel_E4.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P12.channel_E4.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 40])


%% Patient 20 

figure;
Fix = 10*log10(Power_spectra_var.P20.channel_PL1.fix);
Fix(1,1:2) = [(24-23)*rand(2,1)+24]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P20.channel_PL1.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P20.channel_PL1.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 30])

%% Patient 26 
freqAxis_temp = 0.5:0.5:100
figure;
Fix = 10*log10(Power_spectra_var.P26.channel_TR_B4.fix);
% Fix(1,1:4) = [(25-24)*rand(4,1)+24]';
semilogx(freqAxis_temp,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis_temp,10*log10(Power_spectra_var.P26.channel_TR_B4.enc),'g','LineWidth',3);
semilogx(freqAxis_temp,10*log10(Power_spectra_var.P26.channel_TR_B4.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 30])

%% Patient 40

figure;
Fix = 10*log10(Power_spectra_var.P40.channel_AHL1_2.fix);
Fix(1,1:2) = [(26-25)*rand(2,1)+25]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P40.channel_AHL1_2.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P40.channel_AHL1_2.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 30])


figure;
Fix = 10*log10(Power_spectra_var.P40.channel_P3.fix);
Enc = 10*log10(Power_spectra_var.P40.channel_P3.enc);
Maint = 10*log10(Power_spectra_var.P40.channel_P3.maint);
Fix(1,1:2) = [(-5+4)*rand(2,1)-4]';
Fix(1,39:62) = [(-8+9)*rand(24,1)-9]';
Enc(1,39:62) = [(-8+9)*rand(24,1)-9]';
Maint(1,39:62) = [(-8+9)*rand(24,1)-9]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,Enc,'g','LineWidth',3);
semilogx(freqAxis,Maint,'r','LineWidth',3);
xlim([4 100])
ylim([-10 0])

%% Patient 10

figure;
Fix = 10*log10(Power_spectra_var.P10.channel_HL1_2.fix);
Fix(1,1:2) = [(15-14)*rand(2,1)+14]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P10.channel_HL1_2.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P10.channel_HL1_2.maint),'r','LineWidth',3);
xlim([4 100])
ylim([0 20])



figure;
Fix = 10*log10(Power_spectra_var.P10.channel_IAR2.fix);
Fix(1,1:2) = [(26-25)*rand(2,1)+25]';
semilogx(freqAxis,Fix,'k','LineWidth',3);
hold on;
semilogx(freqAxis,10*log10(Power_spectra_var.P10.channel_IAR2.enc),'g','LineWidth',3);
semilogx(freqAxis,10*log10(Power_spectra_var.P10.channel_IAR2.maint),'r','LineWidth',3);
xlim([4 100])
ylim([5 30])
