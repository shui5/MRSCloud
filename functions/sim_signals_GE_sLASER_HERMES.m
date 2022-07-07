function [MRS_opt, outA, outB, outC, outD] = sim_signals_GE_sLASER_HERMES(MRS_opt)

for ii = 1:length(MRS_opt)

    TE         = MRS_opt(ii).TEs{1}; % echo time [ms]
    refTp      = MRS_opt(ii).refTp;  % duration of the refocusing pulse [ms]
    metabolite = MRS_opt(ii).metab;
    editTp     = MRS_opt(ii).editTp;

    fprintf('\n%s simulation at TE %d ms\n\n', metabolite, TE);
    out_name = [MRS_opt(ii).vendor{1} '_' MRS_opt(ii).seq{1} '_' MRS_opt(ii).localization{1} '_TE' num2str(TE) '_' metabolite '.mat'];
    save_dir = MRS_opt(ii).save_dir;

    editOnFreq1 = MRS_opt(ii).editOnFreq1;
    editOnFreq2 = MRS_opt(ii).editOnFreq2;
    editOnFreq3 = MRS_opt(ii).editOnFreq3;
    editOnFreq4 = MRS_opt(ii).editOnFreq4;

    %taus = [5.0688, abs(5.0688 - 24.3619), (38.3882-24.3619), (43.0007-38.3882), (49.6813-43.0007), (64.3619-49.6813), (80.0-64.3619)];
    taus = [12.8600   14.4940    6.3440    5.5420   14.6320   20.0250    8.1030]; % MM

    % ********************SET UP SIMULATION**********************************

    % Calculate new delays by subtracting the pulse durations from the taus
    % vector;
    delays = zeros(size(taus));

    %delays(1) = taus(1) - (refTp/2);
    %delays(2) = taus(2) - ((refTp+editTp)/2);
    %delays(3) = taus(3) - ((editTp+refTp)/2);
    %delays(4) = taus(4) - ((refTp+refTp)/2);
    %delays(5) = taus(5) - ((refTp+refTp)/2);
    %delays(6) = taus(6) - ((refTp+editTp)/2);
    %delays(7) = taus(7) - (editTp/2);

    % MM
    delays(1) = taus(1) - (editTp/2);
    delays(2) = taus(2) - ((refTp+editTp)/2);
    delays(3) = taus(3) - ((refTp+refTp)/2);
    delays(4) = taus(4) - ((refTp+refTp)/2);
    delays(5) = taus(5) - ((refTp+editTp)/2);
    delays(6) = taus(6) - ((editTp+refTp)/2);
    delays(7) = taus(7) - (refTp/2);

    MRS_opt(ii).delays = delays;
    MRS_opt(ii).taus   = taus;

    MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1, MRS_opt(ii).editTp, (MRS_opt(ii).centreFreq - editOnFreq1) * MRS_opt(ii).Bfield * MRS_opt(ii).gamma/1e6); %editOnFreq1
    MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H, MRS_opt(ii).editRFonA, MRS_opt(ii).editTp, MRS_opt(ii).edit_flipAngle, 0);
    
    MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2, MRS_opt(ii).editTp, (MRS_opt(ii).centreFreq - editOnFreq2) * MRS_opt(ii).Bfield * MRS_opt(ii).gamma/1e6); %editOnFreq2
    MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H, MRS_opt(ii).editRFonB, MRS_opt(ii).editTp, MRS_opt(ii).edit_flipAngle, 0);
    
    MRS_opt(ii).editRFonC  = rf_freqshift(MRS_opt(ii).editRF3, MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq - editOnFreq3) * MRS_opt(ii).Bfield * MRS_opt(ii).gamma/1e6);  %editOnFreq3
    MRS_opt(ii).QoutONC    = calc_shapedRF_propagator_edit(MRS_opt(ii).H, MRS_opt(ii).editRFonC, MRS_opt(ii).editTp, MRS_opt(ii).edit_flipAngle, 0);
    
    MRS_opt(ii).editRFonD  = rf_freqshift(MRS_opt(ii).editRF4, MRS_opt(ii).editTp, (MRS_opt(ii).centreFreq - editOnFreq4) * MRS_opt(ii).Bfield * MRS_opt(ii).gamma/1e6); %editOnFreq4
    MRS_opt(ii).QoutOND    = calc_shapedRF_propagator_edit(MRS_opt(ii).H, MRS_opt(ii).editRFonD, MRS_opt(ii).editTp, MRS_opt(ii).edit_flipAngle, 0);

    %% Start the sequence
    % EXCITE
    d = sim_excite(MRS_opt(ii).d, MRS_opt(ii).H, 'x');
    d = sim_apply_pfilter(d, MRS_opt(ii).H, -1);
    d = sim_evolve(d, MRS_opt(ii).H, delays(1)/1000);

    % EDIT 1
    d1_edit1 = apply_propagator_edit(d, MRS_opt(ii).H, MRS_opt(ii).QoutONA);
    d2_edit1 = apply_propagator_edit(d, MRS_opt(ii).H, MRS_opt(ii).QoutONB);
    d3_edit1 = apply_propagator_edit(d, MRS_opt(ii).H, MRS_opt(ii).QoutONC);
    d4_edit1 = apply_propagator_edit(d, MRS_opt(ii).H, MRS_opt(ii).QoutOND);

    d1_edit1 = sim_apply_pfilter(d1_edit1, MRS_opt(ii).H, -1);
    d2_edit1 = sim_apply_pfilter(d2_edit1, MRS_opt(ii).H, -1);
    d3_edit1 = sim_apply_pfilter(d3_edit1, MRS_opt(ii).H, -1);
    d4_edit1 = sim_apply_pfilter(d4_edit1, MRS_opt(ii).H, -1);

    d1_edit1 = sim_evolve(d1_edit1, MRS_opt(ii).H, delays(2)/1000);
    d2_edit1 = sim_evolve(d2_edit1, MRS_opt(ii).H, delays(2)/1000);
    d3_edit1 = sim_evolve(d3_edit1, MRS_opt(ii).H, delays(2)/1000);
    d4_edit1 = sim_evolve(d4_edit1, MRS_opt(ii).H, delays(2)/1000);

    %% AFP pulse x-direction
    parfor X = 1:length(MRS_opt(ii).x)
    %for X = 1:length(MRS_opt(ii).x) % uncomment this line if parfor is unavailable, scnh

        % RF 1
        d1_x{X} = apply_propagator_refoc(d1_edit1, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X}); %#ok<*PFBNS> 
        d2_x{X} = apply_propagator_refoc(d2_edit1, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});
        d3_x{X} = apply_propagator_refoc(d3_edit1, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});
        d4_x{X} = apply_propagator_refoc(d4_edit1, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});

        d1_x{X} = sim_apply_pfilter(d1_x{X}, MRS_opt(ii).H, +1);
        d2_x{X} = sim_apply_pfilter(d2_x{X}, MRS_opt(ii).H, +1);
        d3_x{X} = sim_apply_pfilter(d3_x{X}, MRS_opt(ii).H, +1);
        d4_x{X} = sim_apply_pfilter(d4_x{X}, MRS_opt(ii).H, +1);

        d1_x{X} = sim_evolve(d1_x{X}, MRS_opt(ii).H, delays(3)/1000);
        d2_x{X} = sim_evolve(d2_x{X}, MRS_opt(ii).H, delays(3)/1000);
        d3_x{X} = sim_evolve(d3_x{X}, MRS_opt(ii).H, delays(3)/1000);
        d4_x{X} = sim_evolve(d4_x{X}, MRS_opt(ii).H, delays(3)/1000);


        % RF 2
        d1_x{X} = apply_propagator_refoc(d1_x{X}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});
        d2_x{X} = apply_propagator_refoc(d2_x{X}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});
        d3_x{X} = apply_propagator_refoc(d3_x{X}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});
        d4_x{X} = apply_propagator_refoc(d4_x{X}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{X});

        d1_x{X} = sim_apply_pfilter(d1_x{X}, MRS_opt(ii).H, -1);
        d2_x{X} = sim_apply_pfilter(d2_x{X}, MRS_opt(ii).H, -1);
        d3_x{X} = sim_apply_pfilter(d3_x{X}, MRS_opt(ii).H, -1);
        d4_x{X} = sim_apply_pfilter(d4_x{X}, MRS_opt(ii).H, -1);

        d1_x{X} = sim_evolve(d1_x{X}, MRS_opt(ii).H, delays(4)/1000);
        d2_x{X} = sim_evolve(d2_x{X}, MRS_opt(ii).H, delays(4)/1000);
        d3_x{X} = sim_evolve(d3_x{X}, MRS_opt(ii).H, delays(4)/1000);
        d4_x{X} = sim_evolve(d4_x{X}, MRS_opt(ii).H, delays(4)/1000);

    end

    d_Ax = struct([]);
    d_Bx = struct([]);
    d_Cx = struct([]);
    d_Dx = struct([]);

    % Calculate the average density matrix (Doing this inside a separate for
    % loop because I couldn't figure out how to do this inside the parfor loop):
    for X = 1:length(MRS_opt(ii).x)
        d_Ax = sim_dAdd(d_Ax, d1_x{X});
        d_Bx = sim_dAdd(d_Bx, d2_x{X});
        d_Cx = sim_dAdd(d_Cx, d3_x{X});
        d_Dx = sim_dAdd(d_Dx, d4_x{X});
    end

    %% AFP pulse y-direction
    parfor Y = 1:length(MRS_opt(ii).y)
    %for Y = 1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox

        % RF 3
        d1_y{Y} = apply_propagator_refoc(d_Ax, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d2_y{Y} = apply_propagator_refoc(d_Bx, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d3_y{Y} = apply_propagator_refoc(d_Cx, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d4_y{Y} = apply_propagator_refoc(d_Dx, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});

        d1_y{Y} = sim_apply_pfilter(d1_y{Y}, MRS_opt(ii).H, +1);
        d2_y{Y} = sim_apply_pfilter(d2_y{Y}, MRS_opt(ii).H, +1);
        d3_y{Y} = sim_apply_pfilter(d3_y{Y}, MRS_opt(ii).H, +1);
        d4_y{Y} = sim_apply_pfilter(d4_y{Y}, MRS_opt(ii).H, +1);

        d1_y{Y} = sim_evolve(d1_y{Y}, MRS_opt(ii).H, delays(5)/1000);
        d2_y{Y} = sim_evolve(d2_y{Y}, MRS_opt(ii).H, delays(5)/1000);
        d3_y{Y} = sim_evolve(d3_y{Y}, MRS_opt(ii).H, delays(5)/1000);
        d4_y{Y} = sim_evolve(d4_y{Y}, MRS_opt(ii).H, delays(5)/1000);


        % EDIT 2
        d1_y{Y} = apply_propagator_edit(d1_y{Y}, MRS_opt(ii).H, MRS_opt(ii).QoutONA);
        d2_y{Y} = apply_propagator_edit(d2_y{Y}, MRS_opt(ii).H, MRS_opt(ii).QoutONB);
        d3_y{Y} = apply_propagator_edit(d3_y{Y}, MRS_opt(ii).H, MRS_opt(ii).QoutONC);
        d4_y{Y} = apply_propagator_edit(d4_y{Y}, MRS_opt(ii).H, MRS_opt(ii).QoutOND);

        d1_y{Y} = sim_apply_pfilter(d1_y{Y}, MRS_opt(ii).H, +1);
        d2_y{Y} = sim_apply_pfilter(d2_y{Y}, MRS_opt(ii).H, +1);
        d3_y{Y} = sim_apply_pfilter(d3_y{Y}, MRS_opt(ii).H, +1);
        d4_y{Y} = sim_apply_pfilter(d4_y{Y}, MRS_opt(ii).H, +1);

        d1_y{Y} = sim_evolve(d1_y{Y}, MRS_opt(ii).H, delays(6)/1000);
        d2_y{Y} = sim_evolve(d2_y{Y}, MRS_opt(ii).H, delays(6)/1000);
        d3_y{Y} = sim_evolve(d3_y{Y}, MRS_opt(ii).H, delays(6)/1000);
        d4_y{Y} = sim_evolve(d4_y{Y}, MRS_opt(ii).H, delays(6)/1000);


        % RF 4
        d1_y{Y} = apply_propagator_refoc(d1_y{Y}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d2_y{Y} = apply_propagator_refoc(d2_y{Y}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d3_y{Y} = apply_propagator_refoc(d3_y{Y}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});
        d4_y{Y} = apply_propagator_refoc(d4_y{Y}, MRS_opt(ii).H, MRS_opt(ii).Qrefoc{Y});

        d1_y{Y} = sim_apply_pfilter(d1_y{Y}, MRS_opt(ii).H, -1);
        d2_y{Y} = sim_apply_pfilter(d2_y{Y}, MRS_opt(ii).H, -1);
        d3_y{Y} = sim_apply_pfilter(d3_y{Y}, MRS_opt(ii).H, -1);
        d4_y{Y} = sim_apply_pfilter(d4_y{Y}, MRS_opt(ii).H, -1);

        d1_y{Y} = sim_evolve(d1_y{Y}, MRS_opt(ii).H, delays(7)/1000);
        d2_y{Y} = sim_evolve(d2_y{Y}, MRS_opt(ii).H, delays(7)/1000);
        d3_y{Y} = sim_evolve(d3_y{Y}, MRS_opt(ii).H, delays(7)/1000);
        d4_y{Y} = sim_evolve(d4_y{Y}, MRS_opt(ii).H, delays(7)/1000);

    end % end of spatial loop (parfor) in y-direction

    d_Ay = struct([]);
    d_By = struct([]);
    d_Cy = struct([]);
    d_Dy = struct([]);

    for Y = 1:length(MRS_opt(ii).y)
        d_Ay = sim_dAdd(d_Ay, d1_y{Y});
        d_By = sim_dAdd(d_By, d2_y{Y});
        d_Cy = sim_dAdd(d_Cy, d3_y{Y});
        d_Dy = sim_dAdd(d_Dy, d4_y{Y});
    end

    %% Now running the tail end of the seqeunce
    [outA, outB, outC, outD] = sim_mega_slaser_shaped_ultrafast_readout(MRS_opt(ii), d_Ay, d_By, d_Cy, d_Dy);

    numSims = MRS_opt(ii).nX * MRS_opt(ii).nY;

    outA = op_ampScale(outA, 1/numSims);
    outB = op_ampScale(outB, 1/numSims);
    outC = op_ampScale(outC, 1/numSims);
    outD = op_ampScale(outD, 1/numSims);

    % Scale by the total size of the simulated region, relative to the size
    % of the voxel.
    if MRS_opt(ii).fovX > MRS_opt(ii).thkX
        voxRatio = (MRS_opt(ii).thkX * MRS_opt(ii).thkY) / (MRS_opt(ii).fovX * MRS_opt(ii).fovY);
    else
        voxRatio = 1;
    end

    outA = op_ampScale(outA, 1/voxRatio);
    outB = op_ampScale(outB, 1/voxRatio);
    outC = op_ampScale(outC, 1/voxRatio);
    outD = op_ampScale(outD, 1/voxRatio);

    % Correct residual DC offset
    outA = op_dccorr_FID_A(outA,'p');
    outB = op_dccorr_FID_A(outB,'p');
    outC = op_dccorr_FID_A(outC,'p');
    outD = op_dccorr_FID_A(outD,'p');

    outA.name = metabolite;
    outB.name = metabolite;
    outC.name = metabolite;
    outD.name = metabolite;
    save(fullfile(save_dir, out_name),'outA','outB','outC','outD');

end

end


function [out1, out2, out3, out4] = sim_mega_slaser_shaped_ultrafast_readout(MRS_opt, d1, d2, d3, d4)

% Readout along y (90 degree phase)
out1 = sim_readout(d1, MRS_opt.H, MRS_opt.Npts, MRS_opt.sw, MRS_opt.lw, 90);
out2 = sim_readout(d2, MRS_opt.H, MRS_opt.Npts, MRS_opt.sw, MRS_opt.lw, 90);
out3 = sim_readout(d3, MRS_opt.H, MRS_opt.Npts, MRS_opt.sw, MRS_opt.lw, 90);
out4 = sim_readout(d4, MRS_opt.H, MRS_opt.Npts, MRS_opt.sw, MRS_opt.lw, 90);

% Correct the ppm scale:
out1.ppm = out1.ppm - (4.68 - MRS_opt.centreFreq);
out2.ppm = out2.ppm - (4.68 - MRS_opt.centreFreq);
out3.ppm = out3.ppm - (4.68 - MRS_opt.centreFreq);
out4.ppm = out4.ppm - (4.68 - MRS_opt.centreFreq);

% Fill in structure header fields:
out1.seq = MRS_opt.seq;
out1.te  = sum(MRS_opt.taus);
out1.sim = 'shaped';

out2.seq = MRS_opt.seq;
out2.te  = sum(MRS_opt.taus);
out2.sim = 'shaped';

out3.seq = MRS_opt.seq;
out3.te  = sum(MRS_opt.taus);
out3.sim = 'shaped';

out4.seq = MRS_opt.seq;
out4.te  = sum(MRS_opt.taus);
out4.sim = 'shaped';

% Additional fields for compatibility with FID-A processing tools.
out1.sz=size(out1.specs);
out1.date=date;
out1.dims.t=1;
out1.dims.coils=0;
out1.dims.averages=0;
out1.dims.subSpecs=0;
out1.dims.extras=0;
out1.averages=1;
out1.rawAverages=1;
out1.subspecs=1;
out1.rawSubspecs=1;
out1.flags.writtentostruct=1;
out1.flags.gotparams=1;
out1.flags.leftshifted=0;
out1.flags.filtered=0;
out1.flags.zeropadded=0;
out1.flags.freqcorrected=0;
out1.flags.phasecorrected=0;
out1.flags.averaged=1;
out1.flags.addedrcvrs=1;
out1.flags.subtracted=1;
out1.flags.writtentotext=0;
out1.flags.downsampled=0;
out1.flags.isISIS=0;
out1.seq = char(MRS_opt.seq); % scnh

out2.sz=size(out2.specs);
out2.date=date;
out2.dims.t=1;
out2.dims.coils=0;
out2.dims.averages=0;
out2.dims.subSpecs=0;
out2.dims.extras=0;
out2.averages=1;
out2.rawAverages=1;
out2.subspecs=1;
out2.rawSubspecs=1;
out2.flags.writtentostruct=1;
out2.flags.gotparams=1;
out2.flags.leftshifted=0;
out2.flags.filtered=0;
out2.flags.zeropadded=0;
out2.flags.freqcorrected=0;
out2.flags.phasecorrected=0;
out2.flags.averaged=1;
out2.flags.addedrcvrs=1;
out2.flags.subtracted=1;
out2.flags.writtentotext=0;
out2.flags.downsampled=0;
out2.flags.isISIS=0;
out2.seq = char(MRS_opt.seq); % scnh

out3.sz=size(out3.specs);
out3.date=date;
out3.dims.t=1;
out3.dims.coils=0;
out3.dims.averages=0;
out3.dims.subSpecs=0;
out3.dims.extras=0;
out3.averages=1;
out3.rawAverages=1;
out3.subspecs=1;
out3.rawSubspecs=1;
out3.flags.writtentostruct=1;
out3.flags.gotparams=1;
out3.flags.leftshifted=0;
out3.flags.filtered=0;
out3.flags.zeropadded=0;
out3.flags.freqcorrected=0;
out3.flags.phasecorrected=0;
out3.flags.averaged=1;
out3.flags.addedrcvrs=1;
out3.flags.subtracted=1;
out3.flags.writtentotext=0;
out3.flags.downsampled=0;
out3.flags.isISIS=0;
out3.seq = char(MRS_opt.seq); % scnh

out4.sz=size(out4.specs);
out4.date=date;
out4.dims.t=1;
out4.dims.coils=0;
out4.dims.averages=0;
out4.dims.subSpecs=0;
out4.dims.extras=0;
out4.averages=1;
out4.rawAverages=1;
out4.subspecs=1;
out4.rawSubspecs=1;
out4.flags.writtentostruct=1;
out4.flags.gotparams=1;
out4.flags.leftshifted=0;
out4.flags.filtered=0;
out4.flags.zeropadded=0;
out4.flags.freqcorrected=0;
out4.flags.phasecorrected=0;
out4.flags.averaged=1;
out4.flags.addedrcvrs=1;
out4.flags.subtracted=1;
out4.flags.writtentotext=0;
out4.flags.downsampled=0;
out4.flags.isISIS=0;
out4.seq = char(MRS_opt.seq); % scnh

end



