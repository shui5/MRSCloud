function [MRS_opt,outA,outB,outC,outD] = sim_signals(MRS_opt)
%run the loop for every metabolite
for ii = 1:length(MRS_opt)
    TE             = MRS_opt(ii).TEs{1}; % Echo time [ms]
    metabolite     = MRS_opt(ii).metab;

    sprintf('%s simulation at TE %d ms',metabolite, TE)
    out_name       = [MRS_opt(ii).vendor{1} '_' MRS_opt(ii).seq{1} '_' MRS_opt(ii).localization{1} '_TE' num2str(TE) '_' metabolite '.mat'];   
    %out_name       = [MRS_opt(ii).seq{1} '_' metabolite '_' num2str(TE) '.mat'];   
    save_dir       = MRS_opt(ii).save_dir; % scnh
    
    mega_or_hadam = MRS_opt(ii).seq{1}; % scnh
    switch (mega_or_hadam)
        case 'MEGA'
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            editOnFreq1 = MRS_opt(ii).editOnFreq1;
            editOnFreq2 = MRS_opt(ii).editOnFreq2;
            editOnFreq3 = MRS_opt(ii).editOnFreq3;
            editOnFreq4 = MRS_opt(ii).editOnFreq4;
%         case 'HERMES'
%             editOnFreq1 = MRS_opt(ii).editOnFreq1;
%             editOnFreq2 = MRS_opt(ii).editOnFreq2;
%             editOnFreq3 = MRS_opt(ii).editOnFreq3;
%             editOnFreq4 = MRS_opt(ii).editOnFreq4;
%         case 'HERMES_GABA_GSH_EtOH'
%             editOnFreq1 = MRS_opt(ii).editOnFreq1;
%             editOnFreq2 = MRS_opt(ii).editOnFreq2;
%             editOnFreq3 = MRS_opt(ii).editOnFreq3;
%             editOnFreq4 = MRS_opt(ii).editOnFreq4;
%         case 'HERCULES'
%             editOnFreq1 = MRS_opt(ii).editOnFreq1;
%             editOnFreq2 = MRS_opt(ii).editOnFreq2;
%             editOnFreq3 = MRS_opt(ii).editOnFreq3;
%             editOnFreq4 = MRS_opt(ii).editOnFreq4;
    end
    
    %     if ~strcmp(mega_or_hadam, 'un_edited')
    %         if strcmp(MRS_opt(1).seq, 'HERC')
    %             editOnFreq1 = MRS_opt(ii).editOnFreq1;
    %             dl_split1   = abs(diff([MRS_opt(ii).editOnFreq2,1.9]))*MRS_opt(ii).Bfield*MRS_opt(1).gamma/1e6;
    %             editOnFreq2 = MRS_opt(ii).editOnFreq2;
    %             editOnFreq3 = MRS_opt(ii).editOnFreq3;
    %             dl_split2   = abs(diff([MRS_opt(ii).editOnFreq4,1.9]))*MRS_opt(ii).Bfield*MRS_opt(1).gamma/1e6;
    %             editOnFreq4 = MRS_opt(ii).editOnFreq4;
    %
    %         elseif strcmp(MRS_opt(1).seq, 'HERM')
    %             editOnFreq1 = MRS_opt(ii).editOnFreq1;
    %             dl_split1   = abs(diff([MRS_opt(ii).editOnFreq2,MRS_opt(ii).editOnFreq3]))*MRS_opt(ii).Bfield*MRS_opt(1).gamma/1e6;
    %             editOnFreq2 = MRS_opt(ii).editOnFreq2;
    %             editOnFreq3 = MRS_opt(ii).editOnFreq3;
    %             editOnFreq4 = MRS_opt(ii).editOnFreq4;
    %
    %         else strcmp(MRS_opt(1).seq, 'MEGA');
    %             editOnFreq1 = MRS_opt(ii).editOnFreq1;
    %             editOnFreq2 = MRS_opt(ii).editOnFreq2;
    %         end
    %     end
    
    % scnh
    switch (mega_or_hadam)
        case 'UnEdited'
            %if strcmp(mega_or_hadam, 'UnEdited')
            
            taus    = [MRS_opt.TE1 MRS_opt.TE2];
            %Calculate new delays by subtracting the pulse duration from tau1 and tau2;
            delays  = MRS_opt(ii).delays; %zeros(size(taus));
            %delays(1)=taus(1)-tp;
            delays(1)=taus(1)-(MRS_opt(ii).refTp);
            %delays(2)=taus(2)-tp;
            delays(2)=taus(2)-(MRS_opt(ii).refTp);
            
            if sum(delays<0)
                error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
            end
            if taus(1)<MRS_opt(ii).refTp/1000
                error('ERROR:  Echo-time 1 cannot be less than duration of refocusing pulse! ABORTING!!');
            end
            if taus(2)<MRS_opt(ii).refTp/1000
                error('ERROR:  Echo-time 2 cannot be less than duration of refocusing pulse! ABORTING!!');
            end
            
            MRS_opt(ii).delays     = delays;
            MRS_opt(ii).taus       = taus;
        otherwise
            % ********************PARAMETERS**********************************
            % set up exact timing - based on 'normal' pulse set on Philips 3T - SH 07252019
            taus = [MRS_opt(ii).TE1/2,...                                  %middleEXC pulse to middleREFOC1
                ((TE/2 + MRS_opt(ii).TE1/2)/2)-MRS_opt(ii).TE1/2,...        %middleREFOC1 to middle 1st EDITING
                (TE/2+MRS_opt(ii).TE1/2)-((TE/2 + MRS_opt(ii).TE1/2)/2),... %middle 1st EDITING to middle REFOC2
                (3*TE/4 + MRS_opt(ii).TE1/4)-((TE/2+MRS_opt(ii).TE1/2)),... %middle REFOC2 to middle 2nd EDITING
                TE-(3*TE/4 + MRS_opt(ii).TE1/4)];                       %middle 2nd EDITING to the start of readout
            % ********************SET UP SIMULATION**********************************
            
            %Calculate new delays by subtracting the pulse durations from the taus
            %vector;
            delays    = MRS_opt(ii).delays; %zeros(size(taus));
            delays(1) = taus(1)-(MRS_opt(ii).refTp/2);
            delays(2) = taus(2)-((MRS_opt(ii).refTp+MRS_opt(ii).editTp)/2);
            delays(3) = taus(3)-((MRS_opt(ii).editTp+MRS_opt(ii).refTp)/2);
            delays(4) = taus(4)-((MRS_opt(ii).refTp+MRS_opt(ii).editTp)/2);
            delays(5) = taus(5)-(MRS_opt(ii).editTp/2);
            if sum(delays<0)
                error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
            end
            MRS_opt(ii).delays     = delays;
            MRS_opt(ii).taus       = taus;
    end % scnh
    
    %% Start the sequence
    switch (mega_or_hadam)
        case 'UnEdited'
            d=sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');                                          %EXCITE
            d=sim_apply_pfilter(d,MRS_opt(ii).H,-1);
            d=sim_evolve(d,MRS_opt(ii).H,MRS_opt(ii).delays(1)/2000); % scnh /1000 or /2000 ???
            parfor X=1:length(MRS_opt(ii).x)
            % for X=1:length(MRS_opt(ii).x) % scnh
                d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});  %Refocusing in the X-direction
                d_temp2{X} = d_temp;              
            end
            %calculate the average density matrix (Doing this inside a separate for
            %loop because I couldn't figure out how to do this inside the parfor loop):
            d_A=struct([]);
            for X=1:length(MRS_opt(ii).x)
                d_A         =sim_dAdd(d_A,d_temp2{X});
            end
            d_A       = sim_apply_pfilter(d_A,MRS_opt(ii).H,+1);                   %Apply p_filter
            d_A=sim_evolve(d_A,MRS_opt(ii).H,(delays(1)+delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
            
            MRS_opt(ii).d = d_A; %From here, transfer to the fminsearch loop
            % No editON for PRESS
        otherwise
            
            %if ~strcmp(mega_or_hadam, 'UnEdited') % scnh
            
            d=sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');                                          %EXCITE
            d=sim_apply_pfilter(d,MRS_opt(ii).H,-1);
            d=sim_evolve(d,MRS_opt(ii).H,MRS_opt(ii).delays(1)/1000);
            parfor X=1:length(MRS_opt(ii).x)
                d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});  %Refocusing in the X-direction
                d_temp2{X} = d_temp;
            end
            %calculate the average density matrix (Doing this inside a separate for
            %loop because I couldn't figure out how to do this inside the parfor loop):
            d_A=struct([]);
            for X=1:length(MRS_opt(ii).x)
                d_A         =sim_dAdd(d_A,d_temp2{X});
            end
            d_A       = sim_apply_pfilter(d_A,MRS_opt(ii).H,+1);                   %Apply p_filter
            MRS_opt(ii).d = d_A; %From here, transfer to the fminsearch loop
            
            %Create the propagators for editON
            if strcmp(MRS_opt(ii).seq, 'HERCULES') || strcmp(MRS_opt(ii).seq, 'HERMES') || strcmp(MRS_opt(ii).seq, 'HERMES_GABA_GSH_EtOH')

                    MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq2
                    MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);

                    MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq3
                    MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);

                    MRS_opt(ii).editRFonC  = rf_freqshift(MRS_opt(ii).editRF3,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq3)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %mean([editOnFreq2,editOnFreq4])
                    MRS_opt(ii).QoutONC    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonC,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);

                    MRS_opt(ii).editRFonD  = rf_freqshift(MRS_opt(ii).editRF4,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq4)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %mean([editOnFreq3,editOnFreq4])
                    MRS_opt(ii).QoutOND    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonD,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                
%             elseif strcmp(MRS_opt(ii).seq, 'HERMES')
%                 % MRS_opt(ii).editRF11   = dual_lobe_split_fn(MRS_opt(ii).editRF1,polyval(MRS_opt(ii).split_fit,dl_split1),polyval(MRS_opt(ii).B1power_fit,dl_split1));
%                 % MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF11,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6);
%                 % MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
%                 MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq1
%                 MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
%                 
%                 MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq2
%                 MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
%                 
%                 % MRS_opt(ii).editRFonC  = rf_freqshift(MRS_opt(ii).editRF3,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq3)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq3
%                 % MRS_opt(ii).QoutONC    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonC,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
%                 %MRS_opt(ii).editRF33   = dual_lobe_split_fn(MRS_opt(ii).editRF3,polyval(MRS_opt(ii).split_fit,dl_split1),polyval(MRS_opt(ii).B1power_fit,dl_split1));
%                 MRS_opt(ii).editRFonC  = rf_freqshift(MRS_opt(ii).editRF3,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq3)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq3
%                 MRS_opt(ii).QoutONC    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonC,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
%                 
%                 MRS_opt(ii).editRFonD  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq4)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq4
%                 MRS_opt(ii).QoutOND    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonD,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            else
                MRS_opt(ii).editRFonA  = rf_freqshift(MRS_opt(ii).editRF1,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq1)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq1
                MRS_opt(ii).QoutONA    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonA,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
                
                MRS_opt(ii).editRFonB  = rf_freqshift(MRS_opt(ii).editRF2,MRS_opt(ii).editTp,(MRS_opt(ii).centreFreq-editOnFreq2)*MRS_opt(ii).Bfield*MRS_opt(ii).gamma/1e6); %editOnFreq2
                MRS_opt(ii).QoutONB    = calc_shapedRF_propagator_edit(MRS_opt(ii).H,MRS_opt(ii).editRFonB,MRS_opt(ii).editTp,MRS_opt(ii).edit_flipAngle,0);
            end
    end
    %     else % scnh for PRESS
    %         d=sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');                                          %EXCITE
    %         d=sim_apply_pfilter(d,MRS_opt(ii).H,-1);
    %         d=sim_evolve(d,MRS_opt(ii).H,MRS_opt(ii).delays(1)/2000); % scnh /1000 or /2000 ???
    %         parfor X=1:length(MRS_opt(ii).x)
    %             d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});  %Refocusing in the X-direction
    %             d_temp2{X} = d_temp;
    %         end
    %         %calculate the average density matrix (Doing this inside a separate for
    %         %loop because I couldn't figure out how to do this inside the parfor loop):
    %         d_A=struct([]);
    %         parfor X=1:length(MRS_opt(ii).x)
    %             d_A         =sim_dAdd(d_A,d_temp2{X});
    %         end
    %         d_A       = sim_apply_pfilter(d_A,MRS_opt(ii).H,+1);                   %Apply p_filter
    %         d_A=sim_evolve(d_A,MRS_opt(ii).H,(delays(1)+delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
    %
    %         MRS_opt(ii).d = d_A; %From here, transfer to the fminsearch loop
    %
    %     end % scnh
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start the sequence from delay 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Initialize structures:
    d_A=struct([]);
    
    switch (mega_or_hadam)
        case 'UnEdited'
%            d_A = MRS_opt(ii).d;
%            outA_temp=cell(length(MRS_opt(ii).y));
%            outA=struct([]);
            
%             for Y=1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox
%                 %inA_temp{Y} = d_A;
%                 %outA_temp{Y}   = apply_propagator_refoc(inA_temp{Y},MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});
%                 [outA_temp{Y}]                   = sim_press_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_A);
%                 
%             end %end of spatial loop (parfor) in y direction.
            
        otherwise
            d_A=struct([]);
            d_B=struct([]);
            d_C=struct([]);
            d_D=struct([]);
            % if ~strcmp(MRS_opt(ii).seq, 'UnEdited') % scnh
            %if not PRESS
            if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                [d_A,d_B,d_C,d_D]   = sim_megapress_shaped_ultrafast_edit1(MRS_opt(ii));
            else
                [d_A,d_B,~,~]       = sim_megapress_shaped_ultrafast_edit1(MRS_opt(ii));
            end
            %     else
            % write eveolve functions
    end
    
    % %Initialize structures:
    switch (mega_or_hadam)
        case 'UnEdited'
            outA=struct([]);
            outA_temp=cell(length(MRS_opt(ii).y));
        case 'MEGA'
            outA=struct([]);
            outB=struct([]);
            outA_temp=cell(length(MRS_opt(ii).y));
            outB_temp=cell(length(MRS_opt(ii).y));
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            outA=struct([]);
            outB=struct([]);
            outC=struct([]);
            outD=struct([]);
            outA_temp=cell(length(MRS_opt(ii).y));
            outB_temp=cell(length(MRS_opt(ii).y));
            outC_temp=cell(length(MRS_opt(ii).y));
            outD_temp=cell(length(MRS_opt(ii).y));
    end
    
    %Now loop through y direction (second refoc pulse only);
    switch (mega_or_hadam)
        case 'UnEdited'
            d_A = MRS_opt(ii).d;
            %outA_temp=cell(length(MRS_opt(ii).y));
            %outA=struct([]);
            
            parfor Y=1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox
                %inA_temp{Y} = d_A;
                %outA_temp{Y}   = apply_propagator_refoc(inA_temp{Y},MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});
                [outA_temp{Y}]                   = sim_press_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_A);
            end %end of spatial loop (parfor) in y direction.
            
        otherwise
            parfor Y=1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox
                %                 disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(MRS_opt(ii).y))])
                if ~strcmp(MRS_opt(ii).seq, 'MEGA')
                    [outA_temp{Y},outB_temp{Y},outC_temp{Y},outD_temp{Y}] = sim_megapress_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_A,d_B,d_C,d_D);
                else
                    [outA_temp{Y},outB_temp{Y},~,~]                       = sim_megapress_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_A,d_B);
                    
                end %end of spatial loop (parfor) in y direction.
            end
    end
    
    %     else    % scnh for PRESS
    %                 d_A = MRS_opt(ii).d;
    %                 outA_temp=cell(length(MRS_opt(ii).y));
    %                 outA=struct([]);
    %
    %         for Y=1:length(MRS_opt(ii).y)  %Use this if you do have the MATLAB parallel processing toolbox
    %             %inA_temp{Y} = d_A;
    %             %outA_temp{Y}   = apply_propagator_refoc(inA_temp{Y},MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});
    %                [outA_temp{Y}]                   = sim_press_shaped_ultrafast_Ref2(MRS_opt(ii),MRS_opt(ii).Qrefoc{Y},d_A);
    %
    %         end %end of spatial loop (parfor) in y direction.
    %     end % scnh
    
    %Now combine the outputs;  Again, doing this inside a separate for loop
    %becuase I can't figure out how to do this inside the parfor loop:
    
    %Now combine the outputs;  Again, doing this inside a separate for loop
    %becuase I can't figure out how to do this inside the parfor loop:
    switch (mega_or_hadam)
        case 'UnEdited'
            d_Ay=struct([]);
        case 'MEGA'
            d_Ay=struct([]);
            d_By=struct([]);
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            d_Ay=struct([]);
            d_By=struct([]);
            d_Cy=struct([]);
            d_Dy=struct([]);
    end
    %              d_Ay=struct([]);
    %             d_By=struct([]);
    %         if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         d_Cy=struct([]);
    %         d_Dy=struct([]);
    %     end
 for Y=1:length(MRS_opt(ii).y)   
    switch (mega_or_hadam)
        case 'UnEdited'
             %d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
              d_Ay         =sim_dAdd(d_Ay,outA_temp{Y}); % scnh
              %outA_temp{Y} =sim_apply_pfilter(outA_temp{Y},MRS_opt(ii).H,-1); % scnh
              %outA_temp{Y}=sim_evolve(outA_temp{Y},MRS_opt(ii).H,(delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2
              %[outY]                  = sim_readout(outA_temp{Y},MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);
              %figure (7)
              %hold on
              %plot(outY.ppm,real(outY.specs),'k')
              %Y
        case 'MEGA'
            d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
            d_By         =sim_dAdd(d_By,outB_temp{Y});
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
            d_By         =sim_dAdd(d_By,outB_temp{Y});
            d_Cy         =sim_dAdd(d_Cy,outC_temp{Y});
            d_Dy         =sim_dAdd(d_Dy,outD_temp{Y});
    end
 end   
    %     for Y=1:length(MRS_opt(ii).y)
    %         d_Ay         =sim_dAdd(d_Ay,outA_temp{Y});
    %         d_By         =sim_dAdd(d_By,outB_temp{Y});
    %         if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %             d_Cy         =sim_dAdd(d_Cy,outC_temp{Y});
    %             d_Dy         =sim_dAdd(d_Dy,outD_temp{Y});
    %         end
    %     end
    
    %Now running the tail end of the seqeunce: %pfilter-delay4-editing-delay5-RO
    switch (mega_or_hadam)
        case 'UnEdited'
            d_Ay=sim_apply_pfilter(d_Ay,MRS_opt(ii).H,-1);
            d_Ay=sim_evolve(d_Ay,MRS_opt(ii).H,(delays(2))/2000);                     %Evolve by (delays(1)+delays(2))/2

            [outA]                  = sim_readout(d_Ay,MRS_opt.H,MRS_opt.Npts,MRS_opt.sw,MRS_opt.lw,90);            %Readout along y (90 degree phase);
                          %plot(outA.ppm,real(outA.specs),'r') % scnh
            outA.ppm                = outA.ppm-(4.68-MRS_opt.centreFreq); % sim_readout used 4.65, override ppm scale from 4.65 to 3.0, % scnh
        case 'MEGA'
            [outA,outB,~,~]         = sim_megapress_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By);%,d_Cy,d_Dy);
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            [outA,outB,outC,outD]   = sim_megapress_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By,d_Cy,d_Dy);
    end
    
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         [outA,outB,outC,outD] = sim_megapress_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By,d_Cy,d_Dy);
    %     else
    %         [outA,outB,~,~]       = sim_megapress_shaped_ultrafast_readout(MRS_opt(ii),d_Ay,d_By);
    %     end
    
    %For consistent scaling across different shaped simulations, we need to :
    %1.  Scale down by the total number of simulations run (since these were
    %    all added together.
    numSims=(MRS_opt(ii).nX*MRS_opt(ii).nY);
    switch (mega_or_hadam)
        case 'UnEdited'
            outA=op_ampScale(outA,1/numSims);
        case 'MEGA'
            outA=op_ampScale(outA,1/numSims);
            outB=op_ampScale(outB,1/numSims);
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            outA=op_ampScale(outA,1/numSims);
            outB=op_ampScale(outB,1/numSims);
            outC=op_ampScale(outC,1/numSims);
            outD=op_ampScale(outD,1/numSims);
    end
    
    %     numSims=(MRS_opt(ii).nX*MRS_opt(ii).nY);
    %     outA=op_ampScale(outA,1/numSims);
    %     outB=op_ampScale(outB,1/numSims);
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         outC=op_ampScale(outC,1/numSims);
    %         outD=op_ampScale(outD,1/numSims);
    %     end
    
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
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            outA=op_ampScale(outA,1/voxRatio);
            outB=op_ampScale(outB,1/voxRatio);
            outC=op_ampScale(outC,1/voxRatio);
            outD=op_ampScale(outD,1/voxRatio);
    end
    
    %     outA=op_ampScale(outA,1/voxRatio);
    %     outB=op_ampScale(outB,1/voxRatio);
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         outC=op_ampScale(outC,1/voxRatio);
    %         outD=op_ampScale(outD,1/voxRatio);
    %     end
    
    % Correct residual DC offset
    switch (mega_or_hadam)
        case 'UnEdited'
            outA = op_dccorr_FID_A(outA,'p');
        case 'MEGA'
            outA = op_dccorr_FID_A(outA,'p');
            outB = op_dccorr_FID_A(outB,'p');
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            outA = op_dccorr_FID_A(outA,'p');
            outB = op_dccorr_FID_A(outB,'p');
            outC = op_dccorr_FID_A(outC,'p');
            outD = op_dccorr_FID_A(outD,'p');
    end
    
    %     outA = op_dccorr(outA,'p');
    %     outB = op_dccorr(outB,'p');
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA') && ~strcmp(MRS_opt(ii).seq, 'UnEdited') % scnh
    %         outC = op_dccorr(outC,'p');
    %         outD = op_dccorr(outD,'p');
    %     end
    
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         if strcmp(MRS_opt(ii).metab,'GABA')
    %             h = figure(1);
    %             subplot(4,1,1),plot(outA.ppm,outA.specs),xlim([MRS_opt(ii).ppm_min,MRS_opt(ii).ppm_max]),ylim([-0.3 0.3]);ylabel('subexp A'),title(sprintf('TE of %0.5g ms & editON @ %0.5g, %0.5g, %0.5g, and %0.5g',TE,MRS_opt(ii).editOnFreq1,MRS_opt(ii).editOnFreq2,MRS_opt(ii).editOnFreq3,MRS_opt(ii).editOnFreq4));
    %             subplot(4,1,2),plot(outB.ppm,outB.specs),xlim([MRS_opt(ii).ppm_min,MRS_opt(ii).ppm_max]),ylim([-0.3 0.3]);ylabel('subexp B')
    %             subplot(4,1,3),plot(outC.ppm,outC.specs),xlim([MRS_opt(ii).ppm_min,MRS_opt(ii).ppm_max]),ylim([-0.3 0.3]);ylabel('subexp C')
    %             subplot(4,1,4),plot(outD.ppm,outD.specs),xlim([MRS_opt(ii).ppm_min,MRS_opt(ii).ppm_max]),ylim([-0.3 0.3]);ylabel('subexp D')
    %             pause(0.1)
    %         end
    %     else
    %         edited_signal = op_addScans(outA,outB,1);
    %         b     = -1 * op_integrate(edited_signal,MRS_opt.ppm_min,MRS_opt.ppm_max,'re')/2;
    %         %plots
    %         h = figure(1);
    %         subplot(3,1,1),plot(outA.ppm,outA.specs),xlim([MRS_opt.ppm_min,MRS_opt.ppm_max]),ylim([-0.1 0.15]);ylabel('edit ON'), title(sprintf('TE %0.5g ms & edit ON/OFF @ %0.5g/%0.5g',TE,MRS_opt(ii).editOnFreq1,MRS_opt(ii).editOnFreq2));
    %         subplot(3,1,2),plot(outB.ppm,outB.specs),xlim([MRS_opt.ppm_min,MRS_opt.ppm_max]),ylim([-0.1 0.1]);ylabel('edit OFF')
    %         subplot(3,1,3),plot(outB.ppm,outA.specs-outB.specs),xlim([MRS_opt.ppm_min,MRS_opt.ppm_max]),ylim([-0.1 0.1]);ylabel('DIFF')
    %         pause(0.1)
    %     end
    
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
        case {'HERMES', 'HERCULES', 'HERMES_GABA_GSH_EtOH'}
            outA.name=metabolite;
            outB.name=metabolite;
            outC.name=metabolite;
            outD.name=metabolite;
            %save(out_name,'outA','outB','outC','outD');
            save(fullfile(save_dir,out_name),'outA','outB','outC','outD');
    end
    
    %     if ~strcmp(MRS_opt(ii).seq, 'MEGA')
    %         outA.name=metabolite;
    %         outB.name=metabolite;
    %         outC.name=metabolite;
    %         outD.name=metabolite;
    %         save(out_name,'outA','outB','outC','outD');
    %     else
    %         outA.name=metabolite;
    %         outB.name=metabolite;
    %         save(out_name,'outA','outB');
    %     end
    
end
end

%Nested Function #1, create a dual-lobe pulse with a fixed 20 ms pulse
function [dl_struct] = dual_lobe_split_fn(sl_struct, dl_split, tw1)
%Load sl editing RF waveform
SG_sl               = sl_struct;
SG_sl_dur           = 20; %ms - duration of the editing pulse
time2               = 1:200; %200 points
time2               = time2/200*SG_sl_dur/1000;
%CALCULATE THE SPLITTING
split               = dl_split;
COSCOS              = 2*cos((time2-0.01005)*pi*split);
%Add the COS modulation
SG_dl               = SG_sl;
SG_dl.waveform(:,2) = SG_sl.waveform(:,2).*COSCOS.';
SG_dl.tw1           = tw1;
dl_struct           = SG_dl;
end

%Nested Function #2
function [d1,d2,d3,d4] = sim_megapress_shaped_ultrafast_edit1(MRS_opt)
%CONTINUE PULSE SEQUENCE************
% if not PRESDS
d       = sim_evolve(MRS_opt.d,MRS_opt.H,MRS_opt.delays(2)/1000);                      %Evolve by delays(2)
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
% % else
% d       = sim_evolve(MRS_opt.d,MRS_opt.H,MRS_opt.delays(2)/1000);                      %Evolve by delays(2)
% d1      = sim_apply_pfilter(d_out1,MRS_opt.H,+1);                                      %Apply pfilter
% end

%SEND IT TO THE NEXT DIRECTION**************
end

%Nested Function #3
function [d_out1,d_out2,d_out3,d_out4] = sim_megapress_shaped_ultrafast_Ref2(MRS_opt,Qrefoc,d1,d2,d3,d4)
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

%Nested Function #4
function [d_out1] = sim_press_shaped_ultrafast_Ref2(MRS_opt,Qrefoc,d1)
%CONTINUE PULSE SEQUENCE************
d_out1   = apply_propagator_refoc(d1,MRS_opt.H,Qrefoc);

%SEND IT TO THE TAIL OF THE SEQUENCE**************
end

function [out1,out2,out3,out4] = sim_megapress_shaped_ultrafast_readout(MRS_opt,d1,d2,d3,d4)
%CONTINUE PULSE SEQUENCE************
d1       = sim_apply_pfilter(d1,MRS_opt.H,-1);                                         %Apply pfilter
d2       = sim_apply_pfilter(d2,MRS_opt.H,-1);                                         %Apply pfilter
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_apply_pfilter(d3,MRS_opt.H,-1);                                         %Apply pfilter
    d4       = sim_apply_pfilter(d4,MRS_opt.H,-1);                                         %Apply pfilter
end
d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(4)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(4)
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_evolve(d3,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(4)
    d4       = sim_evolve(d4,MRS_opt.H,MRS_opt.delays(4)/1000);                            %Evolve by delays(4)
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
d1       = sim_evolve(d1,MRS_opt.H,MRS_opt.delays(5)/1000);                            %Evolve by delays(5)
d2       = sim_evolve(d2,MRS_opt.H,MRS_opt.delays(5)/1000);                            %Evolve by delays(5)
if ~strcmp(MRS_opt.seq, 'MEGA')
    d3       = sim_evolve(d3,MRS_opt.H,MRS_opt.delays(5)/1000);                            %Evolve by delays(5)
    d4       = sim_evolve(d4,MRS_opt.H,MRS_opt.delays(5)/1000);                            %Evolve by delays(5)
end
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