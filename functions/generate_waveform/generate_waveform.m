%Nested Function #1, create a dual-lobe pulse with a fixed 20 ms pulse
Bfield              = 3; % 3; 2.89;        % Philips magnetic field strength [Tesla]
gamma               = 42577000; % gyromagnetic ratio
F0                  = Bfield*gamma/1e6;
tp                  = 14; % duration of the editing pulse [ms]
fspan               = 1; %[kHz]
centreFreq          = 4.68;
ac_avg              = mean([4.56,1.9]);
split_fit   = [-4.72454812433615e-39,3.49746593911910e-35,-1.17972529114192e-31,2.41257501068876e-28,-3.35334486561442e-25,3.36289930490681e-22,-2.51925567024359e-19,1.43926948019703e-16,-6.34362195286665e-14,2.16723851396576e-11,-5.73276757923450e-09,1.16606262695043e-06,-0.000179940454083095,0.0206037322492284,-1.68939288729085,93.4579199805247,-3115.57335510316,47271.5463178346];
B1power_fit = [3.04275712634937e-44,-2.20316511739376e-40,7.44769386993364e-37,-1.56125144847875e-33,2.27359544767397e-30,-2.44184493594410e-27,2.00405005381190e-24,-1.28500696403995e-21,6.52644480995629e-19,-2.64612779316882e-16,8.59140538534298e-14,-2.23150800017893e-11,4.61306484949738e-09,-7.51638521375630e-07,9.50378264188365e-05,-0.00910554220880586,0.637054360250509,-30.6191846583449,901.589099210139,-12233.0077131620];

editWaveform1       = 'sg100_200pts.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
dl_split1   = abs(diff([4.56,1.9]))*Bfield*gamma/1e6;
% use parametrized look up to get the correct intended splitting
%editRF11   = dl_generate(editRF1,polyval(split_fit,dl_split1),polyval(B1power_fit,dl_split1));
% use the hard coded splitting as on Philips
editRF11   = dl_generate(editRF1,dl_split1,editRF1.tw1*2);

% rf_bloch simulation
[mv1,sc1]         = rf_blochSim(editRF11,tp,fspan,0);   %dual-lobe pulse
sc_ppm = (sc1*1000/F0) + centreFreq  - (centreFreq-ac_avg); % convert Hz to ppm
figure(2),plot(sc_ppm,mv1(3,:),'b');
io_writepta(editRF11,'dl_Philips_4_56_1_90_14ms.pta');

function [dl_struct] = dl_generate(sl_struct, dl_split, tw1)
%Load sl editing RF waveform
SG_sl               = sl_struct;
SG_sl_dur           = 14; %ms - duration of the editing pulse
time2               = 1:200; %200 points
time2               = time2/200*SG_sl_dur/1000;
time2_mid           = (time2(1)+time2(end))/2;
%CALCULATE THE SPLITTING
split               = dl_split;
COSCOS              = 2*cos((time2-time2_mid)*pi*split);
%Add the COS modulation
SG_dl               = SG_sl;
SG_dl.waveform(:,2) = SG_sl.waveform(:,2).*COSCOS.';
SG_dl.tw1           = tw1;
dl_struct           = SG_dl;
end