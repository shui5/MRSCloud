function [MRS_opt,out] = sim_signals_STEAM(MRS_opt)

for ii = 1:length(MRS_opt)
    
    TE             = MRS_opt.TEs{1}; % Echo time [ms]
    metabolite     = MRS_opt.metab;
    sprintf('%s simulation at TE %d ms',metabolite, TE)
    out_name       = [MRS_opt.vendor{1} '_' MRS_opt.seq{1} '_' MRS_opt.localization{1} '_TE' num2str(TE) '_' metabolite '.mat'];
    save_dir       = MRS_opt.save_dir;
    mega_or_hadam  = MRS_opt.seq;
    Bfield         = MRS_opt.Bfield;
    n              = MRS_opt.Npts;       % number of spectral points
    sw             = MRS_opt.sw;         % spectral width [Hz]
    linewidth      = MRS_opt.lw;         % linewidth of the output spectrum [Hz]
    tm             = MRS_opt.tm;         % mixing time
    centreFreq     = 3.0;
    out            = struct();
    
    %% Start the sequence
    
    d = sim_excite(MRS_opt(ii).d,MRS_opt(ii).H,'x');                                        %EXCITE
    d = sim_apply_pfilter(d,MRS_opt(ii).H,1);
    d = sim_evolve(d,MRS_opt(ii).H,TE/2000);
    d_1 = sim_rotate(d,MRS_opt(ii).H,-90,'x'); % 90 degree pulse about x' axis.
    %             parfor X = 1:length(MRS_opt(ii).x)
    %                 d_temp     = apply_propagator_refoc(d,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{X});  %Refocusing in the X-direction
    %                 d_temp2{X} = d_temp;
    %             end
    %             d_A = struct([]);
    %             for X = 1:length(MRS_opt(ii).x)
    %                 d_A = sim_dAdd(d_A,d_temp2{X});
    %             end
    d_1 = sim_apply_pfilter(d_1,MRS_opt(ii).H,0);
    d_1 = sim_evolve(d_1,MRS_opt(ii).H,tm/1000);
    d_2 = sim_rotate(d_1,MRS_opt(ii).H,tm,'x');                          %Second 90 degree pulse about x' axis.
    %             parfor Y=1:length(MRS_opt(ii).y)
    %                 d_temp     = apply_propagator_refoc(d_y,MRS_opt(ii).H,MRS_opt(ii).Qrefoc{Y});  %Refocusing in the Y-direction
    %                 d_temp2{Y} = d_temp;
    %             end
    %calculate the average density matrix (Doing this inside a separate for
    %loop because I couldn't figure out how to do this inside the parfor loop):
    %             d_A=struct([]);
    %             for X=1:length(MRS_opt(ii).x)
    %                 d_A         =sim_dAdd(d_A,d_temp2{X});
    %             end
    d_2 = sim_apply_pfilter(d_2,MRS_opt(ii).H,-1);                   %Apply p_filter
    d_2 = sim_evolve(d_2,MRS_opt(ii).H,TE/2000);
    [out,dout] = sim_readout(d_2,MRS_opt(ii).H,n,sw,linewidth,90);     %Readout along +y' (90 degree phase);
    
%     numSims=(MRS_opt(ii).nX*MRS_opt(ii).nY);
%     out=op_ampScale(out,1/numSims);
%     
%     if MRS_opt(ii).fovX>MRS_opt(ii).thkX
%         voxRatio=(MRS_opt(ii).thkX*MRS_opt(ii).thkY)/(MRS_opt(ii).fovX*MRS_opt(ii).fovY);
%     else
%         voxRatio=1;
%     end
%     
%     out         = op_ampScale(out,1/voxRatio);
%     out         = op_dccorr_FID_A(out,'p');
    
    out.name    = metabolite;
    out.ppm     = out.ppm-(4.68-centreFreq);
    out.tm      = tm;
    out.seq     = mega_or_hadam;
    out.Bfield  = Bfield;
    
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
    
    save(fullfile(save_dir,out_name),'out');
end
end
