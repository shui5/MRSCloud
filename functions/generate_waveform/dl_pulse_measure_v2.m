clc
% close all
clear all
% addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function/')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/Function/edited_PRESS_sim/Experim_funct')
addpath('/Users/muhammad/Documents/PostDoc_Hopkins/Work/Projects/HERMES_HERCULES/HERCULES/HERC_Optimization/Simulation/pulses')

%% Freq check
% waveform     = 'sl_univ_pulse.pta';
% intend_split = [100:5:650]; %Lowest split 100 hz
% [B1, act_split] = dual_lobe_split(waveform,intend_split)
% close all
% figure(1),plot(intend_split,act_split,'o'),xlabel('intended split, Hz'),ylabel('actual split, Hz')
% figure(2),plot(intend_split,B1,'o'),xlabel('intended split, Hz'),ylabel('tw1')

%% B1 optimizer
waveform                = 'sl_univ_pulse.pta';
SG_sl                   = io_loadRFwaveform(waveform,'inv',0);
SG_sl_dur               = 20; %ms - duration of the editing pulse
time2                   = 1:200; %200 points
time2                   = time2/200*SG_sl_dur/1000;
fspan                   = 1; %kHz
cosine_split            = [56:2:60]; %[100:2:650];
corr_split              = cosine_split; %Lowest split 100 Hz, max 5 ppm (~650 Hz)
B0                      = 127.750605;
gamma                   = 42.557;
split_corr              = zeros(size(cosine_split));
%MRS_struct var
MRS_B1.sl               = SG_sl;
MRS_B1.SG_sl_dur        = SG_sl_dur; %ms - duration of the editing pulse
MRS_B1.time2            = time2; %s
MRS_B1.fspan            = fspan; %kHz
MRS_B1.B0               = B0;
MRS_B1.cosine_split     = cosine_split; %Hz
MRS_B1.corr_split       = cosine_split; %Hz
tic
for ii = 1:length(MRS_B1.cosine_split)
    %%run the initial sim for split measurement fro each ii
    COSCOS                  = 2*cos((time2-(SG_sl_dur/2000))*pi*MRS_B1.cosine_split(ii));
    MRS_B1.dl               = MRS_B1.sl;
    %Add the COS modulation
    MRS_B1.dl.waveform(:,2) = MRS_B1.dl.waveform(:,2).*COSCOS.';
    MRS_B1.dl.tw1           = 2*MRS_B1.sl.tw1;
    [mv3,sc3]               = rf_blochSim(MRS_B1.dl,MRS_B1.SG_sl_dur,MRS_B1.fspan,0);
    close all
    %Checking whether the split is as intended -- 08022018 MGSaleh
    z                       = abs(mv3(3,:)+1);
    lowerbound              = find(min(z)==z);
    diff_freq               = abs(diff(lowerbound)/10);
    dif_in_range            = abs(diff_freq-MRS_B1.cosine_split(ii))/MRS_B1.cosine_split(ii);
    temp_fact               = 1;
    B1_corr                 = 2*MRS_B1.sl.tw1;
    split_freq              = cosine_split(ii);
    if cosine_split(ii) >= 540 % I have to change the limit at larger split values
        temp_lim = 0.0005;
    else
        temp_lim = 0.005;
    end
    while ~(dif_in_range < temp_lim)  %Run this as long as the split is not within 0.5% of the actual split
        [B1_corr]               = dual_lobe_min(MRS_B1,MRS_B1.corr_split(ii),MRS_B1.sl.tw1);
        split_freq              = cosine_split(ii);
        split_freq              = (split_freq/(diff_freq*temp_fact))*split_freq; % This is deteermined from visual analysis of the plot in line 126 -- 06252018 MGSaleh
        COSCOS                  = 2*cos((time2-(SG_sl_dur/2000))*pi*split_freq);
        MRS_B1.dl               = MRS_B1.sl;
        MRS_B1.dl.waveform(:,2) = MRS_B1.sl.waveform(:,2).*COSCOS.';
        MRS_B1.dl.tw1           = B1_corr; %2*MRS_B1.sl.tw1;
        close all
        [mv3,sc3]               = rf_blochSim(MRS_B1.dl,MRS_B1.SG_sl_dur,MRS_B1.fspan,0);
        z                       = abs(mv3(3,:)+1);
        lowerbound              = find(min(z)==z);
        diff_freq2              = abs(diff(lowerbound)/10);
        dif_in_range            = abs(diff_freq2-MRS_B1.cosine_split(ii))/MRS_B1.cosine_split(ii); % NO abs
        temp_fact               = temp_fact + 0.005;
        MRS_B1.corr_split(ii)   = split_freq;
    end
    split_corr(ii)          = split_freq;
    final_split(ii)         = diff_freq2
    B1_optim(ii)            = B1_corr;
    clear split_freq diff_freq2 B1_corr
end
time_sim                 = datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS');
disp(time_sim) %% ~7 min
figure(10),plot(MRS_B1.cosine_split, B1_optim/MRS_B1.SG_sl_dur*1000/gamma),xlabel('cosine split, Hz'),ylabel('B1, uT'),xlim([MRS_B1.cosine_split(1) MRS_B1.cosine_split(end)])
figure(11),plot(MRS_B1.cosine_split, final_split),xlabel('cosine split, Hz'),ylabel('split, Hz'),xlim([MRS_B1.cosine_split(1) MRS_B1.cosine_split(end)])
% save('B1_freq_5hz_optimizer')
% save('B1_freq_2hz_optimizer')

% load('B1_freq_2hz_optimizer')
%Fitting a polynomial degree n for cosine split
split_fit      = polyfit(cosine_split,split_corr,17);
x1             = 100:1:650;
freq_Hz_fit    = polyval(split_fit,x1); 
table_freq     = table(x1.',freq_Hz_fit.','VariableNames',{'cosine_split','act_split'})
figure(12)
plot(split_corr,cosine_split,'o'),xlabel('cosine split J'),ylabel('act split')
hold on
plot(freq_Hz_fit,x1,'r--')
hold off

%Fitting a polynomial degree n for tw1 - Using the sinusoidal one
B1power_fit    = polyfit(cosine_split,B1_optim,19);
x1             = 100:1:650;
tw1_fit        = polyval(B1power_fit,x1); 
table_B1       = table(x1.',tw1_fit.','VariableNames',{'cosine_split','B1'})
figure(13)
plot(cosine_split,B1_optim,'o'),xlabel('cosine split J'),ylabel('tw1')
hold on
plot(x1,tw1_fit,'r--')
% axis([0  5  0  2])
hold off
% 
%Fitting a polynomial degree n for tw1 - flatten the curve
B1_temp          = B1_optim;
B1_temp(18:end)  = B1_optim(end)
B1power_fit2     = polyfit(cosine_split,B1_temp,19);
x1               = 100:1:650;
tw1_fit          = polyval(B1power_fit2,x1); 
table_B1         = table(x1.',tw1_fit.','VariableNames',{'cosine_split','B1'})
figure(4)
plot(cosine_split,B1_temp,'o'),xlabel('cosine split J'),ylabel('tw1')
hold on
plot(x1,tw1_fit,'r--')
% axis([0  5  0  2])
hold off


