clc
close all
clear all
% addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function/')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function/edited_PRESS_sim/Experim_funct')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/pulses')

%Create a struct variable with the pulse parameters
load('B1_freq_2hz_optimizer')
steps_Hz            = 2;               
waveform            = 'sl_univ_pulse.pta';
SG_sl               = io_loadRFwaveform(waveform,'inv',0);
SG_sl_dur           = 20; %ms - duration of the editing pulse
time2               = 1:200; %200 points
time2               = time2/200*SG_sl_dur/1000;
fspan               = 1; %kHz
dl_split            = [100:steps_Hz:650]; %Lowest split 100 Hz, max 5 ppm (~650 Hz)
cosine_split        = split_corr; %
B0                  = 127.750605;
gamma               = 42.557;
%MRS_struct var
MRS_B1.sl           = SG_sl;
MRS_B1.SG_sl_dur    = 20; %ms - duration of the editing pulse
MRS_B1.time2        = time2; %s
MRS_B1.fspan        = fspan; %kHz
MRS_B1.B0           = B0;
MRS_B1.dl_split     = dl_split; %Hz
tic
%% Check B1 and frequency
% for ii = 1:length(cosine_split)
% %     [B1_optim(ii)]          = dual_lobe_min(MRS_B1,MRS_B1.dl_split(ii),MRS_B1.sl.tw1);
%     %%run the sim for split measurement
%     COSCOS                  = 2*cos((time2-0.01005)*pi*polyval(split_fit,MRS_B1.dl_split(ii)));
%     MRS_B1.dl               = MRS_B1.sl;
%     %Add the COS modulation
%     MRS_B1.dl.waveform(:,2) = MRS_B1.dl.waveform(:,2).*COSCOS.';
%     MRS_B1.dl.tw1           = polyval(B1power_fit,MRS_B1.dl_split(ii)); %B1_optim(ii);
%     [mv3,sc3]               = rf_blochSim(MRS_B1.dl,MRS_B1.SG_sl_dur,MRS_B1.fspan,0);
%     %Checking whether the split is as intended -- 08022018 MGSaleh
%     lower_value(ii)         = min(mv3(3,:));
%     z                       = abs(mv3(3,:)+1);
%     lowerbound              = find(min(z)==z);
%     act_split(ii)           = abs(diff(lowerbound)/10);
%     close all
% end
% time_sim            = datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS');
% disp(time_sim) %% ~7 min 
% % figure(1),plot(B1_optim),xlabel('cosine split, Hz'),ylabel('tw1')
% figure(2),plot(dl_split, act_split),xlabel('cosine split, Hz'),ylabel('split, Hz')
% figure(3),plot(dl_split, lower_value),xlabel('cosine split, Hz'),ylabel('lowest Mz value')

%% Freq check
[B1, act_split, lowest_value, dl_wave] = dual_lobe_split(waveform,polyval(split_fit,MRS_B1.dl_split),polyval(B1power_fit,MRS_B1.dl_split))
%                               dual_lobe_split(waveform,                        intended_split,                                      tw1,B0)
close all
figure(4),plot(MRS_B1.dl_split,lowest_value,'o'),xlabel('dl\_split, Hz'),ylabel('inversion (Mz)'),ylim([-1 -0.98])
figure(5),plot(MRS_B1.dl_split,MRS_B1.dl_split-act_split,'o'),xlabel('dl\_split, Hz'),ylabel('dl\_split - act split, Hz')
figure(6),plot(MRS_B1.dl_split,B1,'o'),xlabel('dl\_split split, Hz'),ylabel('tw1')

% [B1_2, act_split2, lowest_value2] = dual_lobe_split(waveform,polyval(split_fit,MRS_B1.dl_split),polyval(B1power_fit2,MRS_B1.dl_split))
% %                               dual_lobe_split(waveform,                        intended_split,                                      tw1,B0)
% figure(7),plot(MRS_B1.dl_split,lowest_value2,'o'),xlabel('intended split, Hz'),ylabel('inversion (Mz)')
% figure(8),plot(MRS_B1.dl_split,MRS_B1.dl_split-act_split2,'o'),xlabel('cosine split, Hz'),ylabel('cosine split - act split, Hz')
% figure(9),plot(MRS_B1.dl_split,B1_2,'o'),xlabel('dl_split split, Hz'),ylabel('tw1')
% 
% %Figures
% figure(10),plot(MRS_B1.dl_split,lowest_value,'ob',MRS_B1.dl_split,lowest_value2,'or'),xlabel('intended split, Hz'),ylabel('inversion (Mz)')
% figure(11),plot(MRS_B1.dl_split,MRS_B1.dl_split-act_split,'ob',MRS_B1.dl_split,MRS_B1.dl_split-act_split2,'or'),xlabel('cosine split, Hz'),ylabel('cosine split - act split, Hz')
% figure(12),plot(MRS_B1.dl_split,B1,'ob',MRS_B1.dl_split,B1_2,'or'),xlabel('dl_split split, Hz'),ylabel('tw1')


