function [MRS_opt,outA,outB] = sim_signals_MRSI(MRS_opt)
%run the loop for every metabolite
for ii = 1:length(MRS_opt)
    TE             = MRS_opt(ii).TEs{1}; % Echo time [ms]
    metabolite     = MRS_opt(ii).metab;
    
    sprintf('%s simulation at TE %d ms',metabolite, TE)
    out_name       = [MRS_opt(ii).vendor{1} '_' MRS_opt(ii).seq{1} '_' MRS_opt(ii).localization{1} '_TE' num2str(TE) '_' metabolite '.mat'];
    save_dir       = MRS_opt(ii).save_dir;
    mega_or_hadam  = MRS_opt(ii).seq{1};
    
    switch (mega_or_hadam)
        case 'Edited_se_MRSI'
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
    end
    
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            taus      = TE;
            factor    = (1-170/200); % for the fraction of the excitation pulse happens after the zero-phase point
            delays(1) = (taus - MRS_opt(ii).refTp - 2*MRS_opt(ii).excTp*factor)/2;
            delays(2) = (taus - MRS_opt(ii).refTp)/2;
            
        case 'Edited_se_MRSI'
            taus      = [TE/4 TE/4 TE/4 TE/4]; % assuming they are equally splitted
            %fully slice-selective sequence
            factor    = (1-170/200); % for the fraction of the excitation pulse happens after the zero-phase point
            delays(1) = taus(1) - MRS_opt(ii).editTp/2 - MRS_opt(ii).excTp*factor;
            delays(2) = taus(2) - MRS_opt(ii).editTp/2 - MRS_opt(ii).refTp/2;
            delays(3) = taus(3) - MRS_opt(ii).editTp/2 - MRS_opt(ii).refTp/2;
            delays(4) = taus(4) - MRS_opt(ii).editTp/2;
    end
    MRS_opt(ii).delays     = delays;
    MRS_opt(ii).taus       = taus;
    MRS_opt(ii).d_eqm      = MRS_opt(ii).d;
    
    %% Start the sequence
    
    % Initialize structures:
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            outA=struct([]);
            outA_temp=cell(length(MRS_opt(ii).y));
        case 'Edited_se_MRSI'
            outA=struct([]);
            outB=struct([]);
            outA_temp=cell(length(MRS_opt(ii).y));
            outB_temp=cell(length(MRS_opt(ii).y));
    end
    
    % fully slice-selective sequence
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
             parfor Y=1:length(MRS_opt(ii).y) %Use this if you do have the MATLAB parallel processing toolbox
                [outA_temp{Y}]     = sim_press_shaped_ultrafast_exc(MRS_opt(ii),MRS_opt(ii).Qexc{Y},MRS_opt(ii).d_eqm);
 %               d_A = MRS_opt(ii).d;
                outA_temp{Y}=sim_apply_pfilter(outA_temp{Y},MRS_opt(ii).H,+1);     %Apply p_filter
                outA_temp{Y}=sim_evolve(outA_temp{Y},MRS_opt(ii).H,(delays(1))/1000);
                [outA_temp{Y}]     = sim_press_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},outA_temp{Y});
            end  %end of spatial loop (parfor) in y direction.
            
        case 'Edited_se_MRSI'
            MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6);
            MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6);
            MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            
            parfor Y=1:length(MRS_opt(ii).y)
                %for Y=1:length(MRS_opt(ii).y)
                [d{Y}]  = sim_press_shaped_ultrafast_exc(MRS_opt(ii),MRS_opt(ii).Qexc{Y},MRS_opt(ii).d_eqm);
                d{Y}    = sim_apply_pfilter(d{Y},MRS_opt(ii).H,+1);
                d{Y}    = sim_evolve(d{Y},MRS_opt(ii).H,(delays(1))/1000);
                [d_out1{Y},d_out2{Y}]       = sim_MRSI_ultrafast_edit1(d{Y},MRS_opt(ii).H,MRS_opt(ii).QoutONA,MRS_opt(ii).QoutONB,+1,delays(2));
                [outA_temp{Y},outB_temp{Y}]   = sim_edited_MRSI_ultrafast_Ref2(MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y},d_out1{Y},d_out2{Y});
            end
    end
    
    %Now combine the outputs;  Again, doing this inside a separate for loop
    %becuase I can't figure out how to do this inside the parfor loop:
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            d_Ay=struct([]);
            
        case {'MEGA','Edited_se_MRSI'}
            d_Ay=struct([]);
            d_By=struct([]);
    end
    
    for Y=1:length(MRS_opt(ii).y)
        switch (mega_or_hadam)
            case 'UnEdited_se_MRSI'
                d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
            case 'Edited_se_MRSI'
                d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
                d_By         =sim_dAdd(d_By,outB_temp{Y});
        end
    end
    
    % Now running the tail end of the seqeunce
        switch (mega_or_hadam)
         case 'UnEdited_se_MRSI'
            d_Ay=sim_apply_pfilter(d_Ay,MRS_opt(ii).H,-1);
            d_Ay=sim_evolve(d_Ay,MRS_opt(ii).H,(delays(2))/1000);                     %Evolve by (delays(1)+delays(2))/2
            [outA]                  = sim_readout(d_Ay,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,270);            %Readout along y (90 degree phase);
            outA.ppm                = outA.ppm-(4.68-MRS_opt.centreFreq); % sim_readout used 4.65, override ppm scale from 4.65 to 3.0, % scnh
            case 'Edited_se_MRSI'
                [outA, outB] = sim_edited_MRSI_ultrafast_readout(MRS_opt(ii),d_Ay,d_By);
        end
    
    % op_ampScale
    numSims=(MRS_opt(ii).nX*MRS_opt(ii).nY);
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            outA=op_ampScale(outA,1/numSims);
        case 'Edited_se_MRSI'
            outA=op_ampScale(outA,1/numSims);
            outB=op_ampScale(outB,1/numSims);
    end
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    if MRS_opt(ii).fovX>MRS_opt(ii).thkX
        voxRatio=(MRS_opt(ii).thkX*MRS_opt(ii).thkY)/(MRS_opt(ii).fovX*MRS_opt(ii).fovY);
    else
        voxRatio=1;
    end
    
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            outA=op_ampScale(outA,1/voxRatio);
        case 'Edited_se_MRSI'
            outA=op_ampScale(outA,1/voxRatio);
            outB=op_ampScale(outB,1/voxRatio);
    end
    
    % Correct residual DC offset
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            outA = op_dccorr_FID_A(outA,'p');
        case 'Edited_se_MRSI'
            outA = op_dccorr_FID_A(outA,'p');
            outB = op_dccorr_FID_A(outB,'p');
    end
    
    switch (mega_or_hadam)
        case 'UnEdited_se_MRSI'
            outA.name=metabolite;
            save(fullfile(save_dir,out_name),'outA');
        case 'Edited_se_MRSI'
            outA.name=metabolite;
            outB.name=metabolite;
            save(fullfile(save_dir,out_name),'outA','outB');
    end
end
end

%Nested Function #1
function [d_out1] = sim_press_shaped_ultrafast_exc(MRS_opt,Qexc,d1)
d_out1   = apply_propagator_exc(d1,MRS_opt.H,Qexc);
end

%Nested Function #2
function [d1,d2] = sim_MRSI_ultrafast_edit1(density_matrix,hamiltonian,Q_A,Q_B,p_select,delay)
% d_out1  = apply_propagator_edit(density_matrix,MRS_opt.H,Q_A);                          %Apply pulse propagator
% d_out2  = apply_propagator_edit(density_matrix,MRS_opt.H,Q_B);                          %Apply pulse propagator
% d_out1      = sim_apply_pfilter(d_out1,MRS_opt.H,p_select);
% d_out2      = sim_apply_pfilter(d_out2,MRS_opt.H,p_select);
% d1      = sim_evolve(d_out1,MRS_opt.H,delay/1000);
% d2      = sim_evolve(d_out2,MRS_opt.H,delay/1000);

d_out1  = apply_propagator_edit(density_matrix,hamiltonian,Q_A);                          %Apply pulse propagator
d_out2  = apply_propagator_edit(density_matrix,hamiltonian,Q_B);                          %Apply pulse propagator
d_out1  = sim_apply_pfilter(d_out1,hamiltonian,p_select);
d_out2  = sim_apply_pfilter(d_out2,hamiltonian,p_select);
d1      = sim_evolve(d_out1,hamiltonian,delay/1000);
d2      = sim_evolve(d_out2,hamiltonian,delay/1000);
end

%Nested Function #3
function [d_out1] = sim_press_shaped_ultrafast_Ref2(MRS_opt,Qrefoc,d1)
%CONTINUE PULSE SEQUENCE************
d_out1   = apply_propagator_refoc(d1,MRS_opt.H,Qrefoc);
end

%Nested Function #4
function [d_out1,d_out2] = sim_edited_MRSI_ultrafast_Ref2(hamiltonian,Qrefoc,d1,d2)
d_out1   = apply_propagator_refoc(d1,hamiltonian,Qrefoc);
d_out2   = apply_propagator_refoc(d2,hamiltonian,Qrefoc);
end

%Nested Function #5
function [out1,out2] = sim_edited_MRSI_ultrafast_readout(MRS_opt,d1,d2)
%CONTINUE PULSE SEQUENCE************

d1       = sim_apply_pfilter(d1,MRS_opt.H,-1);                                         %Apply pfilter
d2       = sim_apply_pfilter(d2,MRS_opt.H,-1);                                         %Apply pfilter

d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(3)/1000);                            %Evolve by delays(4)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(3)/1000);                            %Evolve by delays(4)

d_out1   = apply_propagator_edit(d1,MRS_opt.H,MRS_opt.QoutONA);                        %Apply pulse propagator
d_out2   = apply_propagator_edit(d2,MRS_opt.H,MRS_opt.QoutONB);                        %Apply pulse propagator

d1       = sim_apply_pfilter(d_out1,MRS_opt.H,-1);                                     %Apply pfilter
d2       = sim_apply_pfilter(d_out2,MRS_opt.H,-1);                                     %Apply pfilter

d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(5)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(5)

[out1,~] = sim_readout(d1,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,270);            %Readout along y (90 degree phase);
[out2,~] = sim_readout(d2,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,270);            %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out1.ppm =out1.ppm-(4.68-MRS_opt.centreFreq);
out2.ppm =out2.ppm-(4.68-MRS_opt.centreFreq);

%Fill in structure header fields:
out1.seq=MRS_opt.seq;
out1.te=sum(MRS_opt.taus);
out1.sim='shaped';

out2.seq=MRS_opt.seq;
out2.te=sum(MRS_opt.taus);
out2.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
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

end
