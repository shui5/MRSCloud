function [MRS_opt,outA,outB,outC,outD] = sim_signals_sLASER(MRS_opt)
%run the loop for every metabolite
for ii = 1:length(MRS_opt)
    TE             = MRS_opt(ii).TEs{1}; % Echo time [ms]
    refTp          = MRS_opt(ii).refTp; % duration of the refocusing pulse [ms]
    metabolite     = MRS_opt(ii).metab;
    if ~strcmp(MRS_opt(ii).seq, 'UnEdited')
        editTp         = MRS_opt(ii).editTp;
    end
    
    sprintf('%s simulation at TE %d ms',metabolite, TE)
    out_name       = [MRS_opt(ii).vendor{1} '_' MRS_opt(ii).seq{1} '_' MRS_opt(ii).localization{1} '_TE' num2str(TE) '_' metabolite '.mat'];
    save_dir       = MRS_opt(ii).save_dir;
    
    mega_or_hadam = MRS_opt(ii).seq{1};
    switch (mega_or_hadam)
        case 'MEGA'
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
        case 'HERMES'
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
            editOnFreq3 = MRS_opt(ii).editOnFreq3;
            editOnFreq4 = MRS_opt(ii).editOnFreq4;
        case 'HERCULES'
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
            editOnFreq3 = MRS_opt(ii).editOnFreq3;
            editOnFreq4 = MRS_opt(ii).editOnFreq4;
    end
    
    switch (mega_or_hadam)
        case 'UnEdited'
            %initialize evolution times
            %            tau1=(TE/4-refTp)/2;
            %            tau2=TE/4-refTp;
            tau1=(TE/4)/2;
            tau2=TE/4;
            tau3=TE/4;
            tau4=TE/4;
            tau5=(TE/4)/2;
            taus = [tau1,tau2,tau3,tau4,tau5];
        case 'MEGA'
            % for TE 80 ms
            taus = [4.7523, (23.9861-4.7523), (38.1735-23.9861), (42.8285-38.1735), (49.4073-42.8285), (63.9861-49.4073), (80.0-63.9861)];
        case {'HERMES', 'HERCULES'}
            taus = [5.0688, abs(5.0688 - 24.3619), (38.3882-24.3619), (43.0007-38.3882), (49.6813-43.0007), (64.3619-49.6813), (80.0-64.3619)];
    end
    
    % ********************SET UP SIMULATION**********************************
    
    %Calculate new delays by subtracting the pulse durations from the taus
    %vector;
    mega_or_hadam = MRS_opt(ii).seq{1}; % scnh
    switch (mega_or_hadam)
        case 'UnEdited'
            delays=zeros(size(taus));
            delays(1)=taus(1)-(refTp/2);
            delays(2)=taus(2)-((refTp+refTp)/2);
            delays(3)=taus(3)-((refTp+refTp)/2);
            delays(4)=taus(4)-((refTp+refTp)/2);
            delays(5)=taus(5)-(refTp/2);
            
        case {'MEGA','HERMES','HERCULES'}
            delays=zeros(size(taus));
            delays(1)=taus(1)-(refTp/2);
            delays(2)=taus(2)-((refTp+editTp)/2);
            delays(3)=taus(3)-((editTp+refTp)/2);
            delays(4)=taus(4)-((refTp+refTp)/2);
            delays(5)=taus(5)-((refTp+refTp)/2);
            delays(6)=taus(6)-((refTp+editTp)/2);
            delays(7)=taus(7)-(editTp/2);
    end
    
    MRS_opt(ii).delays     = delays;
    MRS_opt(ii).taus       = taus;
    
    %% Start the sequence, sLASER x-direction
    switch (mega_or_hadam)
        case 'UnEdited'
            d=sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');
            d=sim_apply_pfilter(d,MRS_opt(ii).H,-1);
            d=sim_evolve(d,MRS_opt(ii).H,MRS_opt(ii).delays(1)/1000);
            parfor X=1:length(MRS_opt(ii).x)
                % for X=1:length(MRS_opt(ii).x) % uncomment this line if parfor is unavailable, scnh
                d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});
                d_temp     = sim_apply_pfilter(d_temp,MRS_opt(ii).H,+1);
                d_temp     = sim_evolve(d_temp,MRS_opt(ii).H,(delays(2)/1000));
                d_temp     = apply_propagator_refoc(d_temp,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});
                d_temp     = sim_apply_pfilter(d_temp,MRS_opt(ii).H,-1);
                d_temp2{X} = sim_evolve(d_temp,MRS_opt(ii).H,(delays(3)/1000));
            end
            
            d_A=struct([]);
            for X=1:length(MRS_opt(ii).x)
                d_A         =sim_dAdd(d_A,d_temp2{X});
            end
            
            MRS_opt(ii).d = d_A;
            
        case {'MEGA','HERMES','HERCULES'}
            %Create the propagators for editON
            if strcmp(MRS_opt(ii).seq, 'HERCULES') || strcmp(MRS_opt(ii).seq, 'HERMES')
                MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq2
                MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq3
                MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                MRS_opt(ii).editRFonC  = rf_freqshift(MRS_opt(ii).editRF3,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq3)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %mean([editOnFreq2,editOnFreq4])
                MRS_opt(ii).QoutONC    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonC,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                MRS_opt(ii).editRFonD  = rf_freqshift(MRS_opt(ii).editRF4,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq4)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %mean([editOnFreq3,editOnFreq4])
                MRS_opt(ii).QoutOND    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonD,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            else
                MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq1
                MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq2
                MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            end
            
            d=sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');                                          %EXCITE
            d=sim_apply_pfilter(d,MRS_opt(ii).H,-1);
            d=sim_evolve(d,MRS_opt(ii).H,MRS_opt(ii).delays(1)/1000);
            %parfor X=1:length(MRS_opt(ii).x)
            for X=1:length(MRS_opt(ii).x) % uncomment this line if parfor is unavailable, scnh
                d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});  %Refocusing in the X-direction
                d_temp     = sim_apply_pfilter(d_temp,MRS_opt(ii).H,+1);                   %Apply p_filter % added on Dec 9
                
                %                 MRS_opt(ii).d = d_temp; %From here, transfer to the fminsearch loop
                if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                    [outA_temp{X},outB_temp{X},outC_temp{X},outD_temp{X}] = sim_mega_slaser_shaped_ultrafast_edit1(d_temp,MRS_opt(ii));
                    [outA_temp{X},outB_temp{X},outC_temp{X},outD_temp{X}] = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{X},outA_temp{X},outB_temp{X},outC_temp{X},outD_temp{X});
                else
                    [outA_temp{X},outB_temp{X},~,~]       = sim_mega_slaser_shaped_ultrafast_edit1(d_temp,MRS_opt(ii));
                    [outA_temp{X},outB_temp{X},~,~]       = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{X},outA_temp{X},outB_temp{X});
                end
                                
                switch (mega_or_hadam)
                    case 'MEGA'
                        outA_temp{X}       = sim_apply_pfilter(outA_temp{X},MRS_opt(ii).H,-1);
                        outA_temp{X}       =sim_evolve(outA_temp{X},MRS_opt(ii).H,delays(4)/1000);
                        outB_temp{X}       = sim_apply_pfilter(outB_temp{X},MRS_opt(ii).H,-1);
                        outB_temp{X}        =sim_evolve(outB_temp{X},MRS_opt(ii).H,delays(4)/1000);
                    case {'HERMES', 'HERCULES'}
                        outA_temp{X}      = sim_apply_pfilter(outA_temp{X},MRS_opt(ii).H,-1);
                        outA_temp{X}       =sim_evolve(outA_temp{X},MRS_opt(ii).H,delays(4)/1000);
                        outB_temp{X}       = sim_apply_pfilter(outB_temp{X},MRS_opt(ii).H,-1);
                        outB_temp{X}        =sim_evolve(outB_temp{X},MRS_opt(ii).H,delays(4)/1000);
                        outC_temp{X}       = sim_apply_pfilter(outC_temp{X},MRS_opt(ii).H,-1);
                        outC_temp{X}       =sim_evolve(outC_temp{X},MRS_opt(ii).H,delays(4)/1000);
                        outD_temp{X}       = sim_apply_pfilter(outD_temp{X},MRS_opt(ii).H,-1);
                        outD_temp{X}        =sim_evolve(outD_temp{X},MRS_opt(ii).H,delays(4)/1000);
                end
            end
            %calculate the average density matrix (Doing this inside a separate for
            %loop because I couldn't figure out how to do this inside the parfor loop):
            
            d_Ax=struct([]);
            d_Bx=struct([]);
            d_Cx=struct([]);
            d_Dx=struct([]);
            for X=1:length(MRS_opt(ii).x)
                switch (mega_or_hadam)
                    case 'MEGA'
                        d_Ax         =sim_dAdd(d_Ax,outA_temp{X});
                        d_Bx         =sim_dAdd(d_Bx,outB_temp{X});
                    case {'HERMES', 'HERCULES'}
                        d_Ax         =sim_dAdd(d_Ax,outA_temp{X});
                        d_Bx         =sim_dAdd(d_Bx,outB_temp{X});
                        d_Cx         =sim_dAdd(d_Cx,outC_temp{X});
                        d_Dx         =sim_dAdd(d_Dx,outD_temp{X});
                end
            end
    end
    
    %% AFP pulse y-direction
    %Now loop through y direction (second refoc pulse only);
    switch (mega_or_hadam)
        case 'UnEdited'
            parfor Y=1:length(MRS_opt(ii).y)
                % for X=1:length(MRS_opt(ii).x) % scnh
                d_temp     = apply_propagator_refoc(d_A,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});  %Refocusing in the X-direction
                d_temp     = sim_apply_pfilter(d_temp,MRS_opt(ii).H,+1);                   %Apply p_filter % added on Dec 9
                d_temp     = sim_evolve(d_temp,MRS_opt(ii).H,(delays(4)/1000));                     %Evolve by (delays(1)+delays(2))/2
                d_temp     = apply_propagator_refoc(d_temp,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});  %Refocusing in the X-direction
                d_temp     = sim_apply_pfilter(d_temp,MRS_opt(ii).H,-1);                   %Apply p_filter % added on Dec 9
                d_temp2{Y} = sim_evolve(d_temp,MRS_opt(ii).H,(delays(5)/1000));                     %Evolve by (delays(1)+delays(2))/2
            end
            
            d_Ay=struct([]);
            for Y=1:length(MRS_opt(ii).y)
                d_Ay         =sim_dAdd(d_Ay,d_temp2{Y});
            end
            
            MRS_opt(ii).d = d_Ay; %From here, transfer to the fminsearch loop
            
        case {'MEGA','HERMES','HERCULES'}
            %parfor Y=1:length(MRS_opt(ii).y)
            for Y=1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox
                if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                    [outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y}] = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_Ax,d_Bx,d_Cx,d_Dx);
                else
                    [outA_temp{Y},outB_temp{Y},~,~]                       = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_Ax,d_Bx);
                end
                switch (mega_or_hadam)
                    case 'MEGA'
                        outA_temp{Y}       = sim_apply_pfilter(outA_temp{Y},MRS_opt(ii).H,+1);
                        outA_temp{Y}       =sim_evolve(outA_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                        outB_temp{Y}       = sim_apply_pfilter(outB_temp{Y},MRS_opt(ii).H,+1);
                        outB_temp{Y}        =sim_evolve(outB_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                    case {'HERMES', 'HERCULES'}
                        outA_temp{Y}       = sim_apply_pfilter(outA_temp{Y},MRS_opt(ii).H,+1);
                        outA_temp{Y}       =sim_evolve(outA_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                        outB_temp{Y}       = sim_apply_pfilter(outB_temp{Y},MRS_opt(ii).H,+1);
                        outB_temp{Y}        =sim_evolve(outB_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                        outC_temp{Y}       = sim_apply_pfilter(outC_temp{Y},MRS_opt(ii).H,+1);
                        outC_temp{Y}       =sim_evolve(outC_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                        outD_temp{Y}       = sim_apply_pfilter(outD_temp{Y},MRS_opt(ii).H,+1);
                        outD_temp{Y}        =sim_evolve(outD_temp{Y},MRS_opt(ii).H,delays(5)/1000);
                end
                % 2nd y-direction AFP
                if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                    [outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y}] = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y});
                else
                    [outA_temp{Y},outB_temp{Y},~,~]                       = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},outA_temp{Y},outB_temp{Y});
                end %end of spatial loop (parfor) in y direction.
                % end of 2nd y-direction AFP
                if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                    [outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y}] = sim_mega_slaser_shaped_ultrafast_edit2(MRS_opt(ii),outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y});
                else
                    [outA_temp{Y},outB_temp{Y},~,~]                       = sim_mega_slaser_shaped_ultrafast_edit2(MRS_opt(ii),outA_temp{Y},outB_temp{Y});
                end %end of spatial loop (parfor) in y direction.
            end
    end
    
    switch (mega_or_hadam)
        case 'MEGA'
            d_Ay=struct([]);
            d_By=struct([]);
        case {'HERMES', 'HERCULES'}
            d_Ay=struct([]);
            d_By=struct([]);
            d_Cy=struct([]);
            d_Dy=struct([]);
    end
    
    for Y=1:length(MRS_opt(ii).y)
        switch (mega_or_hadam)
            case 'MEGA'
                d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
                d_By         =sim_dAdd(d_By,outB_temp{Y});
            case {'HERMES', 'HERCULES'}
                d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
                d_By         =sim_dAdd(d_By,outB_temp{Y});
                d_Cy         =sim_dAdd(d_Cy,outC_temp{Y});
                d_Dy         =sim_dAdd(d_Dy,outD_temp{Y});
        end
    end
    
    %% Now running the tail end of the seqeunce: %pfilter-delay4-editing-delay5-RO
    switch (mega_or_hadam)
        case 'UnEdited'
            %d_Ay=sim_apply_pfilter(d_Ay,MRS_opt(ii).H,-1);
            %d_Ay=sim_evolve(d_Ay,MRS_opt(ii).H,(delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
            
            [outA]                  = sim_readout(d_Ay,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
            %plot(outA.ppm,real(outA.specs),'r') % scnh
            outA.ppm                = outA.ppm-(4.68-MRS_opt.centreFreq); % sim_readout used 4.65, override ppm scale from 4.65 to 3.0, % scnh
        case 'MEGA'
            [outA,outB,~,~]         = sim_mega_slaser_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By);%,d_Cy,d_Dy);
        case {'HERMES', 'HERCULES'}
            [outA,outB,outC,outD]   = sim_mega_slaser_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By,d_Cy,d_Dy);
    end
    
    numSims=(MRS_opt(ii).nX*MRS_opt(ii).nY);
    switch (mega_or_hadam)
        case 'UnEdited'
            outA=op_ampScale(outA,1/numSims);
        case 'MEGA'
            outA=op_ampScale(outA,1/numSims);
            outB=op_ampScale(outB,1/numSims);
        case {'HERMES', 'HERCULES'}
            outA=op_ampScale(outA,1/numSims);
            outB=op_ampScale(outB,1/numSims);
            outC=op_ampScale(outC,1/numSims);
            outD=op_ampScale(outD,1/numSims);
    end
    
    %2.  Scale by the total size of the simulated region, relative to the size
    %    of the voxel.
    if MRS_opt(ii).fovX>MRS_opt(ii).thkX
        voxRatio=(MRS_opt(ii).thkX*MRS_opt(ii).thkY)/(MRS_opt(ii).fovX*MRS_opt(ii).fovY);
    else
        voxRatio=1;
    end
    
    switch (mega_or_hadam)
        case 'UnEdited'
            outA=op_ampScale(outA,1/voxRatio);
        case 'MEGA'
            outA=op_ampScale(outA,1/voxRatio);
            outB=op_ampScale(outB,1/voxRatio);
        case {'HERMES', 'HERCULES'}
            outA=op_ampScale(outA,1/voxRatio);
            outB=op_ampScale(outB,1/voxRatio);
            outC=op_ampScale(outC,1/voxRatio);
            outD=op_ampScale(outD,1/voxRatio);
    end
    
    % Correct residual DC offset
    switch (mega_or_hadam)
        case 'UnEdited'
            outA = op_dccorr_FID_A(outA,'p');
        case 'MEGA'
            outA = op_dccorr_FID_A(outA,'p');
            outB = op_dccorr_FID_A(outB,'p');
        case {'HERMES', 'HERCULES'}
            outA = op_dccorr_FID_A(outA,'p');
            outB = op_dccorr_FID_A(outB,'p');
            outC = op_dccorr_FID_A(outC,'p');
            outD = op_dccorr_FID_A(outD,'p');
    end
    
    switch (mega_or_hadam)
        case 'UnEdited'
            outA.name=metabolite;
            %save(out_name,'outA');
            save(fullfile(save_dir,out_name),'outA');
        case 'MEGA'
            outA.name=metabolite;
            outB.name=metabolite;
            %save(out_name,'outA','outB');
            save(fullfile(save_dir,out_name),'outA','outB');
        case {'HERMES', 'HERCULES'}
            outA.name=metabolite;
            outB.name=metabolite;
            outC.name=metabolite;
            outD.name=metabolite;
            %save(out_name,'outA','outB','outC','outD');
            save(fullfile(save_dir,out_name),'outA','outB','outC','outD');
    end
end
end


%Nested Function #1
function [d1,d2,d3,d4] = sim_mega_slaser_shaped_ultrafast_edit1(d_temp,MRS_opt)
%CONTINUE PULSE SEQUENCE************
% if not PRESDS
%d       = sim_evolve(MRS_opt.d,MRS_opt.H,MRS_opt.delays(2)/1000);                      %Evolve by delays(2)
d       = sim_evolve(d_temp,MRS_opt.H,MRS_opt.delays(2)/1000);                      %Evolve by delays(2)
d_out1  = apply_propagator_edit(d,MRS_opt.H,MRS_opt.QoutONA);                          %Apply pulse propagator
d_out2  = apply_propagator_edit(d,MRS_opt.H,MRS_opt.QoutONB);                          %Apply pulse propagator
if ~strcmp(MRS_opt.seq, 'MEGA')
    d_out3  = apply_propagator_edit(d,MRS_opt.H,MRS_opt.QoutONC);                          %Apply pulse propagator
    d_out4  = apply_propagator_edit(d,MRS_opt.H,MRS_opt.QoutOND);
end
%Apply pulse propagator
d1      = sim_apply_pfilter(d_out1,MRS_opt.H,+1);                                      %Apply pfilter
d2      = sim_apply_pfilter(d_out2,MRS_opt.H,+1);                                      %Apply pfilter
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3      = sim_apply_pfilter(d_out3,MRS_opt.H,+1);                                      %Apply pfilter
    d4      = sim_apply_pfilter(d_out4,MRS_opt.H,+1);                                      %Apply pfilter
end
d1      = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(3)/1000);                             %Evolve by delays(3)
d2      = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(3)/1000);                             %Evolve by delays(3)
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3      = sim_evolve(d3,MRS_opt.H,MRS_opt.delays(3)/1000);                             %Evolve by delays(3)
    d4      = sim_evolve(d4,MRS_opt.H,MRS_opt.delays(3)/1000);                             %Evolve by delays(3)
else
    d3      = 0;
    d4      = 0;
end
%SEND IT TO THE NEXT DIRECTION**************
end

%Nested Function #2
function [d_out1,d_out2,d_out3,d_out4] = sim_mega_slaser_shaped_ultrafast_Ref2(MRS_opt,Qrefoc,d1,d2,d3,d4)
%CONTINUE PULSE SEQUENCE************
d_out1   = apply_propagator_refoc(d1,MRS_opt.H,Qrefoc);
d_out2   = apply_propagator_refoc(d2,MRS_opt.H,Qrefoc);
if ~strcmp(MRS_opt.seq, 'MEGA')
    d_out3   = apply_propagator_refoc(d3,MRS_opt.H,Qrefoc);
    d_out4   = apply_propagator_refoc(d4,MRS_opt.H,Qrefoc);
else
    d_out3   = 0;
    d_out4   = 0;
end
%SEND IT TO THE TAIL OF THE SEQUENCE**************
end

%Nested Function #3
function [d1,d2,d3,d4] = sim_mega_slaser_shaped_ultrafast_edit2(MRS_opt,d1,d2,d3,d4)
%CONTINUE PULSE SEQUENCE************
d1       = sim_apply_pfilter(d1,MRS_opt.H,-1);                                         %Apply pfilter
d2       = sim_apply_pfilter(d2,MRS_opt.H,-1);                                         %Apply pfilter
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_apply_pfilter(d3,MRS_opt.H,-1);                                         %Apply pfilter
    d4       = sim_apply_pfilter(d4,MRS_opt.H,-1);                                         %Apply pfilter
end
d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(6)/1000);                            %Evolve by delays(4)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(6)/1000);                            %Evolve by delays(4)
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_evolve(d3,MRS_opt.H,MRS_opt.delays(6)/1000);                            %Evolve by delays(4)
    d4       = sim_evolve(d4,MRS_opt.H,MRS_opt.delays(6)/1000);                            %Evolve by delays(4)
end
d_out1   = apply_propagator_edit(d1,MRS_opt.H,MRS_opt.QoutONA);                        %Apply pulse propagator
d_out2   = apply_propagator_edit(d2,MRS_opt.H,MRS_opt.QoutONB);                        %Apply pulse propagator
if ~strcmp(MRS_opt.seq, 'MEGA')
    d_out3   = apply_propagator_edit(d3,MRS_opt.H,MRS_opt.QoutONC);                        %Apply pulse propagator
    d_out4   = apply_propagator_edit(d4,MRS_opt.H,MRS_opt.QoutOND);                        %Apply pulse propagator
end
d1       = sim_apply_pfilter(d_out1,MRS_opt.H,-1);                                     %Apply pfilter
d2       = sim_apply_pfilter(d_out2,MRS_opt.H,-1);                                     %Apply pfilter
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_apply_pfilter(d_out3,MRS_opt.H,-1);                                     %Apply pfilter
    d4       = sim_apply_pfilter(d_out4,MRS_opt.H,-1);                                     %Apply pfilter
end
d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(7)/1000);                            %Evolve by delays(5)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(7)/1000);                            %Evolve by delays(5)
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_evolve(d3,MRS_opt.H,MRS_opt.delays(7)/1000);                            %Evolve by delays(5)
    d4       = sim_evolve(d4,MRS_opt.H,MRS_opt.delays(7)/1000);                            %Evolve by delays(5)
else
    d3      = 0;
    d4      = 0;
end
end

%Nested Function #4
function [out1,out2,out3,out4] = sim_mega_slaser_shaped_ultrafast_readout(MRS_opt,d1,d2,d3,d4)
[out1,~] = sim_readout(d1,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
[out2,~] = sim_readout(d2,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
if ~strcmp(MRS_opt.seq, 'MEGA')
    [out3,~] = sim_readout(d3,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
    [out4,~] = sim_readout(d4,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
else
    out3     = 0;
    out4     = 0;
end
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out1.ppm =out1.ppm-(4.68-MRS_opt.centreFreq);
out2.ppm =out2.ppm-(4.68-MRS_opt.centreFreq);
if ~strcmp(MRS_opt.seq, 'MEGA')
    out3.ppm =out3.ppm-(4.68-MRS_opt.centreFreq);
    out4.ppm =out4.ppm-(4.68-MRS_opt.centreFreq);
end
%Fill in structure header fields:
out1.seq=MRS_opt.seq;
out1.te=sum(MRS_opt.taus);
out1.sim='shaped';

out2.seq=MRS_opt.seq;
out2.te=sum(MRS_opt.taus);
out2.sim='shaped';

if ~strcmp(MRS_opt.seq, 'MEGA')
    out3.seq=MRS_opt.seq;
    out3.te=sum(MRS_opt.taus);
    out3.sim='shaped';
    
    out4.seq=MRS_opt.seq;
    out4.te=sum(MRS_opt.taus);
    out4.sim='shaped';
end

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

if ~strcmp(MRS_opt.seq, 'MEGA') && ~strcmp(MRS_opt.seq, 'UnEdited')
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
end