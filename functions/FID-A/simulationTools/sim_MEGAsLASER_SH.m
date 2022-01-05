% sim_MEGAsLASER_SH.m

% Written by Steve CN Hui and Muhammad G Saleh, Johns Hopkins University, 2019
% adapted from sim_MP_Siemens.m
% Georg Oeltzschner, Johns Hopkins University, 2018
% adapted from run_simSLaserShaped.m
% Dana Goerzen, Jamie Near, McGill University, 2019
%
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment with semi-LASER localization using Philips product sequence with fully shaped editing
% and refocusing pulses. Phase cycling of both the editing and refocusing
% pulses is performed. Furthermore, simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal. Summation across phase cycles and spatial positions is
% performed. As a result of the phase cycling and spatially resolved simulations,
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations. To enable the use
% of the MATLAB parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes. If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
%
% Parameters description:
%
% refocWaveform     = name of refocusing pulse waveform.
% editWaveform1      = name of editing pulse waveform.
% editOnFreq        = freqeucny of edit on pulse[ppm]
% editOffFreq       = frequency of edit off pulse[ppm]
% refTp             = duration of refocusing pulses[ms]
% editTp            = duration of editing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% thkZ              = slice thickness of z excitation pulse [cm]
% x                 = vector of x positions to simulate [cm]
% y                 = vector of y positions to simulate [cm]
% z                 = vector of z positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]

function [outA,outB] = sim_MEGAsLASER_SH(metabolites,freq_ppms)

for kk=1:numel(metabolites)
    metabolite = metabolites{kk};
    freq_ppm = freq_ppms{kk};
    % ********************PARAMETERS**********************************
    
    % Spin system to simulate
    spinSys     = metabolite;
    out_name    = ['Philips_MEGA-sLASER_' spinSys '.mat']; 
    
    % General properties of the simulations
    Npts    = 8192;     % number of spectral points
    sw      = 4000;     % spectral width [Hz]
    Bfield  = 3;        % Philips magnetic field strength [Tesla]
    lw      = 1;        % linewidth of the output spectrum [Hz]
    gamma   = 42577000; % gyromagnetic ratio
    
    % Define the pulse waveforms here
    %exciteWaveform      = 'sg100_200pts.pta';          % name of excitation pulse waveform.
    refWaveform          = 'sampleAFPpulse_HS2_R15.RF'; % adiabatic RF pulse shaped waveform 
    editWaveform         = 'sampleEditPulse.pta';       % name of 1st single editing pulse waveform. %ON/OFF 
    
    % Define frequency parameters for editing targets
    centreFreq  = 3.0;                 % Center frequency of MR spectrum [ppm]
    % editOnFreq = freq_ppm;                % Center frequency of 1st MEGA-PRESS experiment [ppm]

    % Define pulse durations and flip angles specific for every vendor
    refTp       = 3.5;              % duration of refocusing pulses [ms], set to zero for hard pulse
    editTp1     = 20;               % duration of 1st editing pulse [ms]
    editTp2     = 20;               % duration of 1st editing pulse [ms]
    TE          = 80;               % Echo time [ms]
    TE1         = 5.0688*2          % TE1 [ms]; 
    TE2         = TE - TE1;         % TE2 [ms]
    flips       = [180,180];        % flip angles of first an second refocusing pulses [degrees]
    
    % Define phase cycling pattern - stays the same, there are two editing
    % pulses in MEGA-PRESS with the same frequency, different time
    editPhCyc1=[0 90]; %phase cycling steps for 1st editing pulse [degrees]
    editPhCyc2=[0 90 180 270]; %phase cycling steps for 2nd editing pulse [degrees] From Shaped edit
    ph1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse
    ph2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse
    ph3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse
    ph4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse
  
    % Define spatial resolution of simulation grid
    fovX    = 4.5;  % size of the full simulation Field of View in the x-direction [cm]
    fovY    = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
    fovZ    = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
    thkX    = 3;    % slice thickness of x refocusing pulse [cm]
    thkY    = 3;    % slice thickness of y refocusing pulse [cm]
    thkZ    = 3;    % slice thickness of z excitation pulse [cm]
    nX      = 1;   % number of spatial points to simulate in x direction
    nY      = 1;   % number of spatial points to simulate in y direction
    nZ      = 1;   % number of spatial points to simulate in z direction
    if nX>1
        x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
    else
        x=0;
    end
    if nY>1
        y=linspace(-fovY/2,fovY/2,nY); %X positions to simulate [cm]
    else
        y=0;
    end
    if nZ>1
        z=linspace(-fovZ/2,fovZ/2,nZ); %X positions to simulate [cm]
    else
        z=0;
    end

% set up exact timing - based on 'normal' pulse set on Philips 3T - SH
taus = [5.0688, abs(5.0688 - 24.3619), (38.3882-24.3619), (43.0007-38.3882), (49.6813-43.0007), (64.3619-49.6813), (80.0-64.3619)];
    
    % ********************SET UP SIMULATION**********************************
    
    % Load RF waveforms for excitation and refocusing pulses
    
    %excRF           = io_loadRFwaveform(exciteWaveform,'exc',0);     
    refRF           = io_loadRFwaveform(refWaveform,'inv'); % adiabatic RF pulse shaped waveform 
    editRF_A        = io_loadRFwaveform(editWaveform,'inv',0); %ON/OFF
    
    % Construct the editing pulses from the waveforms and defined
    % frequencies
    editRFonA       =rf_freqshift(editRF_A,editTp1,(centreFreq-freq_ppm(2)).*Bfield.*gamma/1e6); %OFF
    editRFonB       =rf_freqshift(editRF_A,editTp1,(centreFreq-freq_ppm(1)).*Bfield.*gamma/1e6); %1.90ppm

    % Load the spin system definitions and pick the spin system of choice
    load spinSystems
    sys=eval(['sys' spinSys]);
    
    % Resample the refocusing RF pulses to 100 pts to reduce computational workload 
%   excRF = rf_resample(excRF,100);
%   refRF = rf_resample(refRF,100);
    
    % Set up gradients
    
    % Gradient pulse shape generated from mpqspy (PPE) and divided by 10 to
    % convert mT/m to G/cm, scnh and mgs

    Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
    
%% Built from sim_MP_Siemens.m and sim_sLASER_shaped.m
   %Initialize structures:

   outA_posxy_rpc = cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(ph1));
   outB_posxy_rpc = cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(ph1));
   outA_posxy_pre = cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
   outB_posxy_pre = cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
   outA_posxy_post= cell(length(x),length(y));
   outB_posxy_post= cell(length(x),length(y));
   outA           = struct([]);
   outB           = struct([]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   
%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
%for X=1:length(x)
parfor X=1:length(x)
    for Y=1:length(y)
        for EP1=1:length(editPhCyc1)
            for EP2=1:length(editPhCyc2)
                for m=1:length(ph1)
                    disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) '; Y-position ' num2str(Y) ' of ' num2str(length(y))...
                        'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
                        'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
                        '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
                    outA_posxy_rpc{X}{Y}{EP1}{EP2}{m}=sim_sLASER_shaped_SH(Npts,sw,Bfield,lw,sys,TE,refRF,refTp,editRFonA,editTp1,x(X),y(Y),Gx(X),Gy(Y),editPhCyc1(EP1),editPhCyc2(EP2),ph1(m),ph2(m),ph3(m),ph4(m),taus);
                    outB_posxy_rpc{X}{Y}{EP1}{EP2}{m}=sim_sLASER_shaped_SH(Npts,sw,Bfield,lw,sys,TE,refRF,refTp,editRFonB,editTp2,x(X),y(Y),Gx(X),Gy(Y),editPhCyc1(EP1),editPhCyc2(EP2),ph1(m),ph2(m),ph3(m),ph4(m),taus);
        
                    if m==1
                        %out_posxy{X}{Y}=out_posxy_rpc{X}{Y}{m};
                        outA_posxy_pre{X}{Y}{EP1}{EP2} = outA_posxy_rpc{X}{Y}{EP1}{EP2}{m};
                        outB_posxy_pre{X}{Y}{EP1}{EP2} = outB_posxy_rpc{X}{Y}{EP1}{EP2}{m};
                      
                    elseif m==2 || m==3
                        %out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m},1);
                        outA_posxy_pre{X}{Y}{EP1}{EP2} = op_addScans(outA_posxy_pre{X}{Y}{EP1}{EP2},outA_posxy_rpc{X}{Y}{EP1}{EP2}{m},1);
                        outB_posxy_pre{X}{Y}{EP1}{EP2} = op_addScans(outB_posxy_pre{X}{Y}{EP1}{EP2},outB_posxy_rpc{X}{Y}{EP1}{EP2}{m},1);
                      
                    else
                        %out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m});
                        outA_posxy_pre{X}{Y}{EP1}{EP2} = op_addScans(outA_posxy_pre{X}{Y}{EP1}{EP2},outA_posxy_rpc{X}{Y}{EP1}{EP2}{m});
                        outB_posxy_pre{X}{Y}{EP1}{EP2} = op_addScans(outB_posxy_pre{X}{Y}{EP1}{EP2},outB_posxy_rpc{X}{Y}{EP1}{EP2}{m});
                      
                    end
                end
                
                if EP1==1 && EP2==1
                    outA_posxy_post{X}{Y} = outA_posxy_pre{X}{Y}{EP1}{EP2};
                    outB_posxy_post{X}{Y} = outB_posxy_pre{X}{Y}{EP1}{EP2};
              
                else
                    outA_posxy_post{X}{Y} = op_addScans(outA_posxy_post{X}{Y},outA_posxy_pre{X}{Y}{EP1}{EP2});
                    outB_posxy_post{X}{Y} = op_addScans(outB_posxy_post{X}{Y},outB_posxy_pre{X}{Y}{EP1}{EP2});
                
                end 
            end %end of 2nd editing phase cycle loop.
        end %end of 1st editing phase cycle loop.
        
        %out=op_addScans(out,out_posxy{X}{Y});
        outA = op_addScans(outA,outA_posxy_post{X}{Y});
        outB = op_addScans(outB,outB_posxy_post{X}{Y});    

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY*length(ph1)*length(editPhCyc1)*length(editPhCyc2));
%out=op_ampScale(out,1/numSims);
outA = op_ampScale(outA,1/numSims);
outB = op_ampScale(outB,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
if fovX>thkX
    voxRatio=(thkX*thkY)/(fovX*fovY);
else
    voxRatio=1;
end
%out=op_ampScale(out,1/voxRatio);
outA = op_ampScale(outA,1/voxRatio);
outB = op_ampScale(outB,1/voxRatio);
    
%Correct residual DC offset
outA = op_dccorr(outA,'p');
outB = op_dccorr(outB,'p');
    
outA.name = metabolite;
outB.name = metabolite;
    
save(out_name,'outA','outB','outA_posxy_post','outB_posxy_post');
    
end
end %END OF MAIN FUNCTION:  Nested functions below.

%%
function out = sim_sLASER_shaped_SH(n,sw,Bfield,linewidth,sys,te,RF,refTp,editPulse,editTp,dx,dy,Gx,Gy,editPh1,editPh2,ph1,ph2,ph3,ph4,taus,centreFreq,flipAngle)

if nargin<23
    centreFreq=3.0; % Change it to GABA 3 ppm -- SH and MGSaleh
    if nargin<22
        flipAngle=180;
    end
end


editFlip{1}=[0 0 180 180 0 0]; %Cell array of edit-on pulse flip angles.

delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-((refTp+editTp)/2);
delays(3)=taus(3)-((editTp+refTp)/2);
delays(4)=taus(4)-((refTp+refTp)/2);
delays(5)=taus(5)-((refTp+refTp)/2);
delays(6)=taus(6)-((refTp+editTp)/2);
delays(7)=taus(7)-(editTp/2);

%  if (te/4)<(tp/1000)
%      error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
%  end

% % initialize evolution times
%  tau1=(te/4-refTp)/2;
%  tau2=te/4-refTp;
%  tp=refTp;


%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

% %BEGIN sLASER PULSE SEQUENCE************  IDEAL editing pulse
% d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
% d=sim_evolve(d,H,delays(1)/1000);                       %Evolve by tau1
% d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph1,dx,Gx);       %1st shaped 180 degree adiabatic refocusing pulse along X gradient
% d=sim_evolve(d,H,delays(2)/1000);                       %Evolve by tau2
% d=sim_evolve(d,H,editTp/2000);                       %Evolve by tau2
% d=sim_rotate(d,H,editFlip,'y');                        %1st editing pulse
% d=sim_evolve(d,H,editTp/2000);                       %Evolve by tau2
% % d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh1);    %editing pulse along y gradient
% d=sim_evolve(d,H,delays(3)/1000);                       %Evolve by tau2
% d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph2,dx,Gx);       %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
% d=sim_evolve(d,H,delays(4)/1000);                       %Evolve by tau2
% d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph3,dy,Gy);       %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
% d=sim_evolve(d,H,delays(5)/1000);                       %Evolve by tau2
% d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph4,dy,Gy);       %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
% d=sim_evolve(d,H,delays(6)/1000);                       %Evolve by tau1
% d=sim_evolve(d,H,editTp/2000);                       %Evolve by tau2
% d=sim_rotate(d,H,editFlip,'y');                        %1st editing pulse
% % d=sim_evolve(d,H,editTp/1000);                       %Evolve by tau2
% d=sim_evolve(d,H,editTp/2000);                       %Evolve by tau2
% % d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh2);    %editing pulse along y gradient
% d=sim_evolve(d,H,delays(7)/1000);                       %Evolve by tau1

% %BEGIN sLASER PULSE SEQUENCE************  Shaped editing pulse 
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d=sim_evolve(d,H,delays(1)/1000);                       %Evolve by tau1
d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph1,dx,Gx);       %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,delays(2)/1000);                       %Evolve by tau2
d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh1);    %editing pulse along y gradient
d=sim_evolve(d,H,delays(3)/1000);                       %Evolve by tau2
d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph2,dx,Gx);       %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,delays(4)/1000);                       %Evolve by tau2
d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph3,dy,Gy);       %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,delays(5)/1000);                       %Evolve by tau2
d=sim_shapedRF(d,H,RF,refTp,flipAngle,ph4,dy,Gy);       %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,delays(6)/1000);                       %Evolve by tau1
d=sim_shapedRF(d,H,editPulse,editTp,180,90+editPh2);    %editing pulse along y gradient
d=sim_evolve(d,H,delays(7)/1000);                       %Evolve by tau1

%Without editing pulse - SH and MGSaleh 08132019
% d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
% d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
% d=sim_shapedRF(d,H,RF,tp,flipAngle,ph1,dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
% d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
% d=sim_shapedRF(d,H,RF,tp,flipAngle,ph2,dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
% d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
% d=sim_shapedRF(d,H,RF,tp,flipAngle,ph3,dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
% d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
% d=sim_shapedRF(d,H,RF,tp,flipAngle,ph4,dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
% d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1

[out,dout]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isISIS=0;

end


