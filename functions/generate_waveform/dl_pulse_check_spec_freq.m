clc
close all
clear 
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
fspan               = 1; %kHz
B0_MHz              = 127.760501; %MHz
B0                  = 3;     % Philips magnetic field strength [Tesla]
centreFreq          = 3; %ppm
gamma               = 42577000; % gyromagnetic ratio
[mv_sl1,sc_sl1]       = rf_blochSim(SG_sl,SG_sl_dur,fspan,0);
close all
split_fit           = [-4.72454812433615e-39,3.49746593911910e-35,-1.17972529114192e-31,2.41257501068876e-28,-3.35334486561442e-25,3.36289930490681e-22,-2.51925567024359e-19,1.43926948019703e-16,-6.34362195286665e-14,2.16723851396576e-11,-5.73276757923450e-09,1.16606262695043e-06,-0.000179940454083095,0.0206037322492284,-1.68939288729085,93.4579199805247,-3115.57335510316,47271.5463178346];
B1power_fit         = [3.04275712634937e-44,-2.20316511739376e-40,7.44769386993364e-37,-1.56125144847875e-33,2.27359544767397e-30,-2.44184493594410e-27,2.00405005381190e-24,-1.28500696403995e-21,6.52644480995629e-19,-2.64612779316882e-16,8.59140538534298e-14,-2.23150800017893e-11,4.61306484949738e-09,-7.51638521375630e-07,9.50378264188365e-05,-0.00910554220880586,0.637054360250509,-30.6191846583449,901.589099210139,-12233.0077131620];
GABA_freq           = 1.9;
GSH_freq            = [4.58, 4.68];
editOnFreq          = [mean([GSH_freq(1),GABA_freq]), mean([GSH_freq(2),GABA_freq])];
sc_ppm              = (sc_sl1*1000/B0_MHz) + centreFreq;% - (centreFreq-editOnFreq(1));

%% Freq check
SG_sl_shf              = rf_freqshift(SG_sl,SG_sl_dur,(centreFreq-editOnFreq(1))*B0*gamma/1e6);
[mv_sl2,sc_sl2]        = rf_blochSim(SG_sl_shf,SG_sl_dur,fspan,0);


dl_split1              = abs(diff([GSH_freq(1),GABA_freq]))*B0*gamma/1e6;
edit_wave_456_19       = dual_lobe_split_fn(SG_sl,polyval(split_fit,dl_split1),polyval(B1power_fit,dl_split1),SG_sl_dur);
edit_wave_456_19_shf   = rf_freqshift(edit_wave_456_19,SG_sl_dur,-(centreFreq-editOnFreq(1))*B0*gamma/1e6);
[mv_dl1,sc_dl1]        = rf_blochSim(edit_wave_456_19_shf,SG_sl_dur,fspan,0);

dl_split2           = abs(diff([GSH_freq(2),GABA_freq]))*B0*gamma/1e6;
edit_wave_468_19    = dual_lobe_split_fn(SG_sl,polyval(split_fit,dl_split2),polyval(B1power_fit,dl_split2),SG_sl_dur);
[mv_dl2,sc_dl2]     = rf_blochSim(edit_wave_468_19,SG_sl_dur,fspan,0);

%plots
figure(5)
plot(sc_ppm,mv_dl1(3,:))
xl            = xline(4.1,'k');
xl.LineWidth  = 1;
xl2           = xline(4.58,'k--');
xl2.LineWidth = 1;
%,sc_ppm,mv_sl2(3,:),'r')%, ylim([0.95 1])




%% Nested Function #1, create a dual-lobe pulse with a fixed 20 ms pulse
function [dl_struct] = dual_lobe_split_fn(sl_struct, dl_split, tw1, dur)
%Load sl editing RF waveform
SG_sl               = sl_struct;
SG_sl_dur           = dur; %ms - duration of the editing pulse
time2               = 1:200; %200 points
time2               = time2/200*SG_sl_dur/1000;
%CALCULATE THE SPLITTING
split               = dl_split;
COSCOS              = 2*cos((time2-(SG_sl_dur/2000))*pi*split);
%Add the COS modulation
SG_dl               = SG_sl;
SG_dl.waveform(:,2) = SG_sl.waveform(:,2).*COSCOS.';
SG_dl.tw1           = tw1;
dl_struct           = SG_dl;
end
