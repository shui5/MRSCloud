function MRS_opt = load_parameters_1_5T(MRS_opt)

spinSys         = MRS_opt.metab;
FieldStr        = MRS_opt.FieldStr;
vendor          = MRS_opt.vendor;
localization    = MRS_opt.localization;
%metabolite      = MRS_opt.metab;
%localization    = MRS_temp.localization;
% Define spatial resolution of simulation grid
fovX            = 4.5;  % size of the full simulation Field of View in the x-direction [cm]
fovY            = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
fovZ            = 4.5;  % size of the full simulation Field of View in the y-direction [cm]
thkX            = 3;    % slice thickness of x refocusing pulse [cm]
thkY            = 3;    % slice thickness of y refocusing pulse [cm]
thkZ            = 3;    % slice thickness of z excitation pulse [cm]
Npts            = 8192; % 2048 / 8192;     % number of spectral points
sw              = 4000;     % spectral width [Hz]
lw              = 1;        % linewidth of the output spectrum [Hz]
gamma           = 42577000; % gyromagnetic ratio (1H = 42.58 MHz/T)

if strcmp(vendor, 'Siemens')||strcmp(vendor, 'Universal_Siemens')
    Bfield = 1.49;    % Siemens magnetic field strength [Tesla]
else
    Bfield = 1.5;     % Philips magnetic field strength [Tesla]
end

MRS_opt.fovX    = fovX;             % size of the full simulation Field of View in the x-direction [cm]
MRS_opt.fovY    = fovY;             % size of the full simulation Field of View in the y-direction [cm]
MRS_opt.fovZ    = fovZ;             % size of the full simulation Field of View in the y-direction [cm]
MRS_opt.thkX    = thkX;             % slice thickness of x refocusing pulse [cm]
MRS_opt.thkY    = thkY;             % slice thickness of y refocusing pulse [cm]
MRS_opt.thkZ    = thkZ;             % slice thickness of z excitation pulse [cm]
MRS_opt.Npts    = Npts;     % number of spectral points
MRS_opt.sw      = sw;     % spectral width [Hz]
MRS_opt.lw      = lw;        % linewidth of the output spectrum [Hz]
MRS_opt.gamma   = gamma; % gyromagnetic ratio
MRS_opt.Bfield  = Bfield;        % Philips magnetic field strength [Tesla]

% Define the pulse waveforms here
switch(vendor{1})
    case 'Philips'
        if strcmp(MRS_opt.localization, 'PRESS')
            excWaveform             = 'Philips_spredrex.pta';
            refocWaveform           = 'gtst1203_sp.pta';     % name of refocusing pulse waveform.
        else % sLASER GOIA pulse
            refocWaveform           = 'GOIA';                % name of refocusing pulse waveform.
        end
        
        if ~strcmp(MRS_opt.seq, 'UnEdited')&&~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
            if strcmp(MRS_opt.seq, 'HERCULES')
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 2nd single editing pulse waveform. [4.18ppm]
                editWaveform3       = 'dl_Philips_4_58_1_90.pta';   % name of 1st dual editing pulse waveform. [4.58ppm 1.90ppm]
                editWaveform4       = 'dl_Philips_4_18_1_90.pta';   % name of 2nd dual editing pulse waveform. [4.18ppm 1.90ppm]
            elseif strcmp(MRS_opt.seq, 'HERMES')
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [4.56ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 2nd single editing pulse waveform. [1.9ppm]
                editWaveform3       = 'dl_Philips_4_56_1_90.pta';   % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
                editWaveform4       = 'sg100_100_0_14ms_88hz.pta';  % name of non-editing pulse waveform. [non-editing]
            elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
                editWaveform1       = 'dl_Philips_4_56_1_90.pta';
                editWaveform2       = 'dl_Philips_3_67_1_90_20ms.pta';
                editWaveform3       = 'dl_Philips_3_67_4_56_20ms.pta';
                editWaveform4       = 'sg100_100_0_14ms_88hz.pta';
            else                                                    % MEGA PRESS
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [1.9ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [7.5ppm]
                if strcmp(MRS_opt.seq, 'Edited_se_MRSI')
                    editWaveform1       = 'am_sg_150_100.pta';  % name of 1st single editing pulse waveform. [1.9ppm]
                    editWaveform2       = 'am_sg_150_100.pta';  % name of 1st single editing pulse waveform. [7.5ppm]
                end
            end
        end
        
    case 'Siemens'
        if strcmp(MRS_opt.localization, 'PRESS')
            refocWaveform           = 'orig_refoc_mao_100_4.pta';      % name of refocusing pulse waveform.
        else % sLASER GOIA pulse
            refocWaveform           = 'GOIA';                          % name of refocusing pulse waveform.
        end
        
        if ~strcmp(MRS_opt.seq, 'UnEdited')
            if strcmp(MRS_opt.seq, 'HERCULES')
                editWaveform1       = 'Siemens_filtered_editing.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
                editWaveform2       = 'Siemens_filtered_editing.pta';  % name of 1st single editing pulse waveform. [4.18ppm]
                editWaveform3       = 'dl_Siemens_4_58_1_90.pta';      % name of 1st single editing pulse waveform. [4.58ppm 1.90ppm]
                editWaveform4       = 'dl_Siemens_4_18_1_90.pta';      % name of 1st single editing pulse waveform. [4.18ppm 1.90ppm]
            elseif strcmp(MRS_opt.seq, 'HERMES')
                editWaveform1       = 'Siemens_filtered_editing.pta';  % name of 1st single editing pulse waveform. [4.56ppm]
                editWaveform2       = 'Siemens_filtered_editing.pta';  % name of 2nd single editing pulse waveform. [1.90ppm]
                editWaveform3       = 'dl_Siemens_4_56_1_90.pta';      % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
                editWaveform4       = 'Siemens_filtered_editing.pta';  % name of non-editing pulse waveform. [non-editing]
            elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
                editWaveform1       = 'dl_Siemens_4_56_1_90.pta';
                editWaveform2       = 'dl_Siemens_3_67_1_90_20ms.pta';
                editWaveform3       = 'dl_Siemens_3_67_4_56_20ms.pta'; 
                editWaveform4       = 'Siemens_filtered_editing.pta';
            else % MEGA PRESS
                editWaveform1       = 'Siemens_filtered_editing.pta';  % name of 1st single editing pulse waveform. [1.9ppm]
                editWaveform2       = 'Siemens_filtered_editing.pta';  % name of 1st single editing pulse waveform. [7.5ppm]
            end
        end
        
    case 'GE'
        if strcmp(MRS_opt.localization, 'PRESS')
            refocWaveform           = 'GE_rfa167.pta';                % name of refocusing pulse waveform.
        else % sLASER GOIA pulse
            refocWaveform           = 'GE_GOIA_WURST';                % name of refocusing pulse waveform.
        end
        
        if ~strcmp(MRS_opt.seq, 'UnEdited')
            if strcmp(MRS_opt.seq, 'HERCULES')
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta'; %'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [4.58ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta'; %'sg100_100_0_14ms_88hz.pta';  % name of 2nd single editing pulse waveform. [4.18ppm]
                editWaveform3       = 'dl_Philips_4_58_1_90.pta';   % name of 1st dual editing pulse waveform. [4.58ppm 1.90ppm]
                editWaveform4       = 'dl_Philips_4_18_1_90.pta';   % name of 2nd dual editing pulse waveform. [4.18ppm 1.90ppm]
            elseif strcmp(MRS_opt.seq, 'HERMES')
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [4.56ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 2nd single editing pulse waveform. [1.9ppm]
                editWaveform3       = 'dl_Philips_4_56_1_90.pta';   % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
                editWaveform4       = 'sg100_100_0_14ms_88hz.pta';  % name of non-editing pulse waveform. [non-editing]
            elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
                editWaveform1       = 'dl_Philips_4_56_1_90.pta';
                editWaveform2       = 'dl_Philips_3_67_1_90_20ms.pta';
                editWaveform3       = 'dl_Philips_3_67_4_56_20ms.pta';
                editWaveform4       = 'sg100_100_0_14ms_88hz.pta';
            else                                                    % MEGA PRESS
                editWaveform1       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [1.9ppm]
                editWaveform2       = 'sg100_100_0_14ms_88hz.pta';  % name of 1st single editing pulse waveform. [7.5ppm]
            end
        end
        
    case 'Universal_Philips'
        if strcmp(MRS_opt.localization, 'PRESS')
            refocWaveform           = 'univ_eddenrefo.pta';        % name of refocusing pulse waveform.
        else % sLASER GOIA pulse
            refocWaveform           = 'GOIA';                      % name of refocusing pulse waveform.
        end
        
        if ~strcmp(MRS_opt.seq, 'UnEdited')
            if strcmp(MRS_opt.seq, 'HERCULES')
                editWaveform1       = 'sl_univ_pulse.pta';         % name of 1st single editing pulse waveform. [4.58ppm]
                editWaveform2       = 'sl_univ_pulse.pta';         % name of 1st single editing pulse waveform. [4.18ppm]
                editWaveform3       = 'dl_Philips_4_58_1_90.pta';  % name of 1st single editing pulse waveform. [4.58ppm 1.90ppm]
                editWaveform4       = 'dl_Philips_4_18_1_90.pta';  % name of 1st single editing pulse waveform. [4.18ppm 1.90ppm]
            elseif strcmp(MRS_opt.seq, 'HERMES')
                editWaveform1       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [4.56ppm]
                editWaveform2       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [1.90ppm]
                editWaveform3       = 'dl_Philips_univ_4_56_1_90.pta'; % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
                editWaveform4       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [non-editing]
             elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
                editWaveform1       = 'dl_Philips_univ_4_56_1_9_20ms.pta';
                editWaveform2       = 'dl_Philips_univ_3_67_1_9_20ms.pta';
                editWaveform3       = 'dl_Philips_univ_3_67_4_56_20ms.pta';
                editWaveform4       = 'sl_univ_pulse.pta';
            else                                                       % MEGA PRESS
                editWaveform1       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [1.9ppm]
                editWaveform2       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [7.5ppm]
            end
        end
        
    case 'Universal_Siemens'
        if strcmp(MRS_opt.localization, 'PRESS')
            refocWaveform           = 'univ_eddenrefo.pta';        % name of refocusing pulse waveform.
        else % sLASER GOIA pulse
            refocWaveform           = 'GOIA';                      % name of refocusing pulse waveform.
        end
        
        if ~strcmp(MRS_opt.seq, 'UnEdited')
            if strcmp(MRS_opt.seq, 'HERCULES')
                editWaveform1       = 'sl_univ_pulse.pta';         % name of 1st single editing pulse waveform. [4.58ppm]
                editWaveform2       = 'sl_univ_pulse.pta';         % name of 1st single editing pulse waveform. [4.18ppm]
                editWaveform3       = 'dl_Siemens_4_58_1_90.pta';  % name of 1st single editing pulse waveform. [4.58ppm 1.90ppm]
                editWaveform4       = 'dl_Siemens_4_18_1_90.pta';  % name of 1st single editing pulse waveform. [4.18ppm 1.90ppm]
            elseif strcmp(MRS_opt.seq, 'HERMES') 
                editWaveform1       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [4.56ppm]
                editWaveform2       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [1.90ppm]
                editWaveform3       = 'dl_Siemens_4_56_1_90.pta'; % name of 1st dual editing pulse waveform. [4.56ppm 1.90ppm]
                editWaveform4       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [non-editing]
             elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
                editWaveform1       = 'dl_Siemens_univ_4_56_1_9_20ms.pta';
                editWaveform2       = 'dl_Siemens_univ_3_67_1_9_20ms.pta';
                editWaveform3       = 'dl_Siemens_univ_3_67_4_56_20ms.pta';
                editWaveform4       = 'sl_univ_pulse.pta';
            else                                                       % MEGA PRESS
                editWaveform1       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [1.9ppm]
                editWaveform2       = 'sl_univ_pulse.pta';             % name of 1st single editing pulse waveform. [7.5ppm]
            end
        end
end

% Define frequency parameters for editing targets

if ~strcmp(MRS_opt.seq, 'UnEdited')&&~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
    editOnFreq1             = MRS_opt.editON{1};    % Center frequency of 1st editing experiment [ppm]
    editOnFreq2             = MRS_opt.editON{2};    % Center frequency of 2nd editing experiment [ppm]
    if ~strcmp(MRS_opt.seq, 'MEGA')&&~strcmp(MRS_opt.seq, 'Edited_se_MRSI')
        editOnFreq3         = MRS_opt.editON{3};    % Center frequency of 3rd HERMES/HERCULES experiment [ppm]
        editOnFreq4         = MRS_opt.editON{4};    % Center frequency of 4th HERMES/HERCULES experiment [ppm]
    end
    MRS_opt.editOnFreq1 = editOnFreq1;
    MRS_opt.editOnFreq2 = editOnFreq2;
    if ~strcmp(MRS_opt.seq, 'MEGA')&&~strcmp(MRS_opt.seq, 'Edited_se_MRSI')
        MRS_opt.editOnFreq3 = editOnFreq3;
        MRS_opt.editOnFreq4 = editOnFreq4;
    end
end

% Define pulse durations and flip angles specific for every vendor
% if strcmp(vendor, 'Siemens')
%     refTp = 7.0;         % Siemens univ_eddenrefo, duration of refocusing pulses [ms], set to zero for hard pulse
% else
%     refTp = 6.8944;      % Philips gtst1203_sp and univ_eddenrefo
% end
%
% MRS_opt.refTp       = refTp;              % duration of refocusing pulses [ms], set to zero for hard pulse

if ~strcmp(MRS_opt.seq, 'UnEdited')&&~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
    switch(vendor{1})
        case 'Philips'
            TE1                 = 7.0*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
        case 'Siemens'
            TE1                 = 7.75*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
        case 'GE'
            TE1                 = 7.35*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
        case {'Universal_Philips','Universal_Siemens'}
            TE1                 = 6.55*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
    end
else
    if strcmp(MRS_opt.seq, 'UnEdited')
        switch(vendor{1})  
        case 'Philips'
            TE1                 = 8.15*2;%7.41*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
            MRS_opt.TE2         = MRS_opt.TEs{1}-TE1; % scnh setup TE2
        case 'Siemens'
            TE1                 = 8.2*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
            MRS_opt.TE2         = MRS_opt.TEs{1}-TE1; % scnh setup TE2
        case 'GE'
            TE1                 = 7.35*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
            MRS_opt.TE2         = MRS_opt.TEs{1}-TE1; % scnh setup TE2
        case {'Universal_Philips','Universal_Siemens'}
            TE1                 = 7.0*2;           % TE1 [ms]
            MRS_opt.TE1         = TE1;              % TE1 [ms]
            MRS_opt.TE2         = MRS_opt.TEs{1}-TE1; % scnh setup TE2
        end
    end
end

if ~strcmp(MRS_opt.seq, 'UnEdited')&&~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
    % TE = str2double(string(MRS_opt.TEs));
    % TE = cell2mat(MRS_opt.TEs); % use this line in old versions of MATLAB
    %    if TE<80
    %     if TE==68
    %         editTp1             = 14;               % If TE<80 for MEGA, editTp is 14 [ms]
    %         editTp2             = 14;               % If TE<80 for MEGA, editTp is 14 [ms]
    %         %    else
    %     elseif TE==80
    %         editTp1             = 20;               % duration of 1st editing pulse [ms]
    %         editTp2             = 20;               % duration of 2nd editing pulse [ms]
    %     else
    if strcmp(MRS_opt.seq, 'MEGA')||strcmp(MRS_opt.seq, 'Edited_se_MRSI')
        editTp1             = MRS_opt.editTp;   % duration of 1st editing pulse [ms]
        editTp2             = MRS_opt.editTp;   % duration of 2nd editing pulse [ms]
        %     end
    else% ~strcmp(MRS_opt.seq, 'MEGA')
        editTp1             = 20;               % duration of 1st editing pulse [ms]
        editTp2             = 20;               % duration of 2nd editing pulse [ms]
        editTp3             = 20;               % duration of 2nd editing pulse [ms]
        editTp4             = 20;               % duration of 2nd editing pulse [ms]
    end
    MRS_opt.editTp      = editTp1;
    MRS_opt.editTp1     = editTp1;               % duration of 1st editing pulse [ms]
    MRS_opt.editTp2     = editTp2;               % duration of 2nd editing pulse [ms]
    if ~strcmp(MRS_opt.seq, 'MEGA')&&~strcmp(MRS_opt.seq, 'Edited_se_MRSI')
        MRS_opt.editTp3     = editTp3;               % duration of 2nd editing pulse [ms]
        MRS_opt.editTp4     = editTp4;               % duration of 2nd editing pulse [ms]
    end
end

if MRS_opt.nX>1
    MRS_opt.x       = linspace(-MRS_opt.fovX/2,MRS_opt.fovX/2,MRS_opt.nX); %X positions to simulate [cm]
else
    MRS_opt.x=0;
end
if MRS_opt.nY>1
    MRS_opt.y       = linspace(-MRS_opt.fovY/2,MRS_opt.fovY/2,MRS_opt.nY); %X positions to simulate [cm]
else
    MRS_opt.y=0;
end
if MRS_opt.nZ>1
    MRS_opt.z       = linspace(-MRS_opt.fovZ/2,MRS_opt.fovZ/2,MRS_opt.nZ); %X positions to simulate [cm]
else
    MRS_opt.z=0;
end

if strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')||strcmp(MRS_opt.seq, 'Edited_se_MRSI')
    switch excWaveform
        case 'Philips_spredrex.pta'
            excRF     = io_loadRFwaveform(excWaveform,'exc',0);
            excTp     = 7.13; % Check duration
            excRF.tbw = 2.277*excTp; %BW99 (kHz) * dur (ms)%1.354
            excRF=rf_resample(excRF,100);
    end
end

switch refocWaveform
    case 'gtst1203_sp.pta'
        refRF     = io_loadRFwaveform(refocWaveform,'ref',0);
        refTp     = 6.89;
        refRF.tbw = 1.354*refTp; %BW50 (kHz) * dur (ms)%1.354 %1.412
        refRF=rf_resample(refRF,100);
    case 'orig_refoc_mao_100_4.pta'
        refRF     = io_loadRFwaveform(refocWaveform,'ref',0);
        refTp     = 4.4;
        %refRF.tbw = 1.31*refTp;  %BW50 (kHz) * dur (ms)
    case 'GE_rfa167.pta'
        refRF     = io_loadRFwaveform(refocWaveform,'ref',0);
        refTp     = 6.5;
        %refRF.tbw = 1.11*refTp;  %BW50 (kHz) * dur (ms)
        refRF=rf_resample(refRF,100);
    case 'univ_eddenrefo.pta'
        refRF     = io_loadRFwaveform(refocWaveform,'ref',0);
        refTp     = 7.0;
        %refRF.tbw = 1.342*refTp; %BW50 (kHz) * dur (ms)
    case 'GOIA'
        load Philips_GOIA_WURST_100pts.mat;
        refRF      = Sweep2;
        refTp      = 4.5; % (ms)
        BW         = 8; % (kHz)
        refRF.tbw  = refTp*BW; % dur (ms) * BW50 (kHz)
        refRF.f0   = 0;
        %refRF.isGM = 1; %is the pulse gradient mdoulated? - 02262020 SH
        refRF.tthk = MRS_opt.thkX*(refTp/1000); %This is the time x sliceThickness product for gradient modulated pulses.  It is in units [cm.s]
    case 'GE_GOIA_WURST'
        load GE_GOIA_WURST_100pts.mat;
        refRF      = Sweep_GE_100;
        refTp      = 4.5; % (ms)
        BW         = 10; % (kHz)
        refRF.tbw  = refTp*BW; % dur (ms) * BW50 (kHz)
        refRF.f0   = 0;
        %refRF.isGM = 1;
        refRF.tthk = MRS_opt.thkX*(refTp/1000); %This is the time x sliceThickness product for gradient modulated pulses.  It is in units [cm.s]
        
end

if strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')||strcmp(MRS_opt.seq, 'Edited_se_MRSI')
    MRS_opt.excRF       = excRF;
    MRS_opt.excTp       = excTp;
end

MRS_opt.refRF       = refRF;
MRS_opt.refTp       = refTp;
% Load RF waveforms for editing pulses, defined up
if strcmp(MRS_opt.seq, 'HERCULES')
    editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
    editRF2             = io_loadRFwaveform(editWaveform2,'inv',0);
    editRF3             = io_loadRFwaveform(editWaveform3,'inv',0);
    editRF4             = io_loadRFwaveform(editWaveform4,'inv',0);
    editRF3.tw1         = editRF2.tw1*2; %1.8107;
    editRF4.tw1         = editRF2.tw1*2; %1.8107;
    MRS_opt.editRF1     = editRF1;
    MRS_opt.editRF2     = editRF2;
    MRS_opt.editRF3     = editRF3;
    MRS_opt.editRF4     = editRF4;
    
elseif strcmp(MRS_opt.seq, 'HERMES')
    editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
    editRF2             = io_loadRFwaveform(editWaveform2,'inv',0);
    editRF3             = io_loadRFwaveform(editWaveform3,'inv',0);
    editRF3.tw1         = editRF2.tw1*2;
    editRF4             = io_loadRFwaveform(editWaveform4,'inv',0);
    MRS_opt.editRF1     = editRF1;
    MRS_opt.editRF2     = editRF2;
    MRS_opt.editRF3     = editRF3;
    MRS_opt.editRF4     = editRF4;
elseif strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
    editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
    editRF2             = io_loadRFwaveform(editWaveform2,'inv',0);
    editRF3             = io_loadRFwaveform(editWaveform3,'inv',0);
    editRF4             = io_loadRFwaveform(editWaveform4,'inv',0);
    editRF1.tw1         = 1.8115; %editRF4.tw1*2;
    editRF2.tw1         = 1.8354;
    editRF3.tw1         = 1.9165;
    editRF4.tw1         = editRF4.tw1;
    MRS_opt.editRF1     = editRF1;
    MRS_opt.editRF2     = editRF2;
    MRS_opt.editRF3     = editRF3;
    MRS_opt.editRF4     = editRF4;
elseif strcmp(MRS_opt.seq, 'MEGA') ||strcmp(MRS_opt.seq, 'Edited_se_MRSI') 
    editRF1             = io_loadRFwaveform(editWaveform1,'inv',0);
    editRF2             = io_loadRFwaveform(editWaveform2,'inv',0);
    MRS_opt.editRF1     = editRF1;
    MRS_opt.editRF2     = editRF2;
end

% Construct the editing pulses from the waveforms and defined
% frequencies
if ~strcmp(MRS_opt.seq, 'UnEdited')&&~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
    MRS_opt.editRFonA   = rf_freqshift(editRF1,editTp1,(MRS_opt.centreFreq-editOnFreq1)*Bfield*gamma/1e6); %1.90 = GABA ON
    MRS_opt.editRFonB   = rf_freqshift(editRF2,editTp2,(MRS_opt.centreFreq-editOnFreq2)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
    if ~strcmp(MRS_opt.seq, 'MEGA')&&~strcmp(MRS_opt.seq, 'Edited_se_MRSI')
        MRS_opt.editRFonC   = rf_freqshift(editRF3,editTp3,(MRS_opt.centreFreq-editOnFreq3)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
        MRS_opt.editRFonD   = rf_freqshift(editRF4,editTp4,(MRS_opt.centreFreq-editOnFreq4)*Bfield*gamma/1e6); %7.5 = GABA OFF or MM supp
    end
    % HERCULES has the same editing pulse duration and timing for all sub-experiments. Valid for Mega-press as well:
        MRS_opt.editTp      = editTp1;
end

% % Load the spin system definitions and pick the spin system of choice
load spinSystems
sys=eval(['sys' spinSys]);

% Set up gradients
if strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
    MRS_opt.Gx          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm] %this is for refoc
    MRS_opt.Gx_exc      = (excRF.tbw/(excTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm] %this is for excitation
elseif strcmp(MRS_opt.seq, 'Edited_se_MRSI')
    MRS_opt.Gx          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm] %this is for refoc
    MRS_opt.Gx_exc      = -(excRF.tbw/(excTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm] %this is for excitation
else
MRS_opt.Gx          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm]
MRS_opt.Gy          = (refRF.tbw/(refTp/1000))/(gamma*MRS_opt.thkY/10000); %[G/cm]
%MRS_opt.Gx_exc      = (excRF.tbw/(excTp/1000))/(gamma*MRS_opt.thkX/10000); %[G/cm] %this is for excitation
end

for k=1:length(sys)
    sys(k).shifts   = sys(k).shifts-MRS_opt.centreFreq;
end
MRS_opt.sys=sys;

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
if ~strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')
taus           = [TE1/2];               %middle 2nd EDITING to the start of readout
%taus                = tausA;
delays         = zeros(size(taus));
%delays(1)      = taus(1)-(refTp/2);
%delays(2)      = 0; %taus(2)-((refTp+editTp)/2);
%delays(3)      = 0; %taus(3)-((editTp+refTp)/2);
%delays(4)      = 0; %taus(4)-((refTp+editTp)/2);
%delays(5)      = 0; %taus(5)-(editTp/2);

MRS_opt.taus   = taus;
MRS_opt.delays = delays;
end
%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian_mgs(sys,Bfield);

MRS_opt.H = H;
MRS_opt.d = d;
%Creating propagators for editing pulse
% if ~strcmp(MRS_opt.seq, 'UnEdited')
%     [MRS_opt.QoutONA]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonA,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
%     [MRS_opt.QoutONB]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonB,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
%     if ~strcmp(MRS_opt.seq, 'MEGA')
%         [MRS_opt.QoutONC]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonC,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
%         [MRS_opt.QoutOND]  = calc_shapedRF_propagator_edit(MRS_opt.H,MRS_opt.editRFonD,MRS_opt.editTp,MRS_opt.edit_flipAngle,0);
%     end
% end

%Creating propagators for Refoc pulse ONLY in the X direction With an assumption x=y, i.e., the voxel is isotropic
%(H,RF,Tp,flipAngle,phase,dfdx,grad)

%% load excitation pulse
%excWaveform = 'univ_spreddenrex.pta';
excWaveform = 'Philips_spredrex.pta';
RF_struct = io_loadRFwaveform (excWaveform,'exc',0);
MRS_opt.excRF = RF_struct;

if strcmp(MRS_opt.seq, 'UnEdited_se_MRSI')||strcmp(MRS_opt.seq, 'Edited_se_MRSI')
    parfor X=1:length(MRS_opt.x)  %Use this if you do have the MATLAB parallel processing toolbox
        [Qexc{X}]  = calc_shapedRF_propagator_exc(MRS_opt.H,MRS_opt.excRF,MRS_opt.excTp,MRS_opt.excflipAngle,0,MRS_opt.y(X),MRS_opt.Gx_exc);
    end
    MRS_opt.Qexc = Qexc;
end

if strcmp(MRS_opt.seq, 'UnEdited')||strcmp(MRS_opt.seq, 'MEGA')||strcmp(MRS_opt.seq, 'HERMES')||strcmp(MRS_opt.seq, 'HERCULES')||strcmp(MRS_opt.seq, 'HERMES_GABA_GSH_EtOH')
    parfor X=1:length(MRS_opt.x)  %Use this if you do have the MATLAB parallel processing toolbox
        [Qrefoc{X}]  = calc_shapedRF_propagator_refoc(MRS_opt.H,MRS_opt.refRF,MRS_opt.refTp,MRS_opt.flipAngle,    0,MRS_opt.y(X),MRS_opt.Gx);
    end
    MRS_opt.Qrefoc = Qrefoc;
end
