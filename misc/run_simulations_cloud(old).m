function basis = run_simulations_cloud(json_input)

addpath(genpath('/Users/steve/Desktop/sim_ultrafast'));
json_input      = '/Users/steve/Desktop/sim_ultrafast/simMRS.json';
sim_paras_json  = loadjson(json_input);

% write parameters to .json file at data path
metablist       = sim_paras_json.userInput.metablist;
vendor          = sim_paras_json.userInput.vendor;
mega_or_hadam   = sim_paras_json.userInput.mega_or_hadam;
localization    = sim_paras_json.userInput.localization;
TE              = sim_paras_json.userInput.TE;   % TE for UnEdited and MEGA only, HERMES and HERCULES are internally fixed at TE=80 ms

% initialize metab parameters

% metablist = {'Ala','Asc','Asp','Cit','Cr','EA','EtOH','GABA','GPC',...
%             'GSH','Gln','Glu','Gly','H2O','Ins','Lac','NAA','NAAG','PCh',...
%             'PCr','PE','Phenyl','Ref0ppm','Scyllo','Ser','Tau','Tyros','bHB','bHG'};

% metablist       = {'GABA'};
% mega_or_hadam   = {'MEGA'}; %Options: UnEdited, MEGA, HERMES or HERCULES
% localization    = {'PRESS'};    %Options: PRESS or sLASER
% vendor          = {'Philips'};  %Options: Philips, Philips_universal or Siemens

flipAngle       = 180;  % Flip angle
centreFreq      = 3.0;  % Center frequency of MR spectrum [ppm]
edit_flipAngle  = 180;  % Edited flip angle
spatial_points  = 101;  % number of spatial points to simulate

for iii = 1:length(metablist)
    %if strcmp(mega_or_hadam, 'un_edited') %PRESS
    MRS_temp.flipAngle      = flipAngle;        % Flip angle degrees
    MRS_temp.centreFreq     = centreFreq;       % Center frequency of MR spectrum [ppm]
    MRS_temp.edit_flipAngle = edit_flipAngle;   % Edited flip angle
    MRS_temp.nX             = spatial_points;               % number of spatial points to simulate in x direction
    MRS_temp.nY             = spatial_points;               % number of spatial points to simulate in y direction
    MRS_temp.nZ             = spatial_points;               % number of spatial points to simulate in z direction
    
    switch(mega_or_hadam{1})
        case 'UnEdited'
            metab                  = metablist(iii);
            Nmetab                 = length(metab);   %Only metabolite to be run
            ppm_min                = [1.1]; %check
            ppm_max                = [1.5]; %check
            %     create a struct variable
            MRS_temp.seq           = mega_or_hadam; %MRS_temp.HERC_or_HERM    = HERC_or_HERM;
            MRS_temp.metab         = metab{Nmetab}; %{ii};
            MRS_temp.Nmetab        = Nmetab; %ii;
            MRS_temp.localization  = localization;
            MRS_temp.vendor        = vendor;
            MRS_temp.TEs           = num2cell(TE);
            MRS_temp.ppm_min       = ppm_min(Nmetab); %(ii);
            MRS_temp.ppm_max       = ppm_max(Nmetab); %(ii);
            MRS_temp               = load_parameters(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
            MRS_opt                = MRS_temp; % Creating a struct variable with dimens >=1;
        case 'MEGA'
            metab                  = metablist(iii);
            Nmetab                 = length(metab); %1;         %Only metabolite to be run
            ppm_min                = [1.1];
            ppm_max                = [1.5];
            A                      = 1.9;       %single-lobe pulse
            B                      = 7.50;      %single-lobe pulse
            %     create a struct variable
            MRS_temp.seq           = mega_or_hadam; %MRS_temp.HERC_or_HERM    = HERC_or_HERM;
            MRS_temp.metab         = metab{Nmetab};
            MRS_temp.Nmetab        = Nmetab;
            MRS_temp.localization  = localization;
            MRS_temp.vendor        = vendor;
            MRS_temp.TEs           = num2cell(TE);
            MRS_temp.ppm_min       = ppm_min(Nmetab);
            MRS_temp.ppm_max       = ppm_max(Nmetab);
            MRS_temp.editON        = num2cell([A B]);
            MRS_temp               = load_parameters(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
            MRS_opt                = MRS_temp; % Creating a struct variable with dimens >=1;
        case 'HERMES'
            metab                  = metablist(iii); %{'GABA'};  %{'GABA','GSH','Lac'};
            TE                     = 80;
            ppm_min                = [2.8 2.75];%
            ppm_max                = [3.2 3.15];%
            A                      = 4.56;          %single-lobe pulse
            B                      = 1.90;          %single-lobe pulse
            C                      = (4.56+1.9)/2;  %dual-lobe pulse
            D                      = 7.50;          %single-lobe pulse
            Nmetab                 = length(metab); %Only metabolite to be run
        case 'HERCULES'
            metab                  = metablist(iii); %{'GABA'};  %{'GABA','GSH','Lac'};
            TE                     = 80;
            ppm_min                = [2.8 2.75];
            ppm_max                = [3.2 3.15];
            A                      = 4.58;          %single-lobe pulse
            B                      = 4.18;          %single-lobe pulse
            C                      = (4.58+1.9)/2;  %dual-lobe pulse
            D                      = (4.18+1.9)/2;  %dual-lobe pulse
            Nmetab                 = length(metab); %Only metabolite to be run
    end
    
    if ~strcmp(mega_or_hadam, 'MEGA') && ~strcmp(mega_or_hadam, 'UnEdited')
        %       for ii = 1:Nmetab
        %           clear MRS_temp
        %             MRS_temp.flipAngle       = flipAngle;        % Flip angle degrees
        %             MRS_temp.centreFreq      = centreFreq;       % Center frequency of MR spectrum [ppm]
        %             MRS_temp.edit_flipAngle  = edit_flipAngle;   % Edited flip angle
        %             MRS_temp.nX              = spatial_points;               % number of spatial points to simulate in x direction
        %             MRS_temp.nY              = spatial_points;               % number of spatial points to simulate in y direction
        %             MRS_temp.nZ              = spatial_points;               % number of spatial points to simulate in z direction
        MRS_temp.seq             = mega_or_hadam; %MRS_temp.HERC_or_HERM    = HERC_or_HERM;
        %MRS_temp.metab           = metab{ii};
        %MRS_temp.Nmetab          = ii;
        MRS_temp.metab         = metab{Nmetab}; %{ii}; %scnh
        MRS_temp.Nmetab        = Nmetab; %ii; % scnh
        MRS_temp.localization    = localization;
        MRS_temp.vendor          = vendor;
        MRS_temp.TEs             = num2cell(TE);
        MRS_temp.editON          = num2cell([A B C D]);
        MRS_temp.ppm_min         = ppm_min;%(ii); %(ii);
        MRS_temp.ppm_max         = ppm_max;%(ii); %(ii);
        MRS_temp                 = load_parameters(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
        %            MRS_opt(ii)              = MRS_temp; % Creating a struct variable with dimens >=1;
        MRS_opt              = MRS_temp; % Creating a struct variable with dimens >=1;
        
        %        end
    end
    
    %Simulate
    if strcmp(mega_or_hadam, 'HERMES') || strcmp(mega_or_hadam, 'HERCULES')
        [MRS_opt,outA, outB, outC, outD]  = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
    elseif strcmp(mega_or_hadam, 'MEGA')
        [MRS_opt,outA, outB]              = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
    else
        [MRS_opt,outA]                    = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
    end
end % for metablist
%% plot figure, remove in cloud version
% switch (mega_or_hadam{1})
%     case 'MEGA'
%         figure(1),
%         hold on
%         plot(outA.ppm, outA.specs), set(gca,'XDir','reverse'), xlim ([1 5]), title([MRS_opt.seq{1} ' ' MRS_temp.metab ' TE' num2str(TE)]),legend
%         hold off
%         x_lim = [1 5];
%     case 'HERMES' %'HERMES') || strcmp(mega_or_hadam, 'HERCULES')
        % y_lim = [-0.2 0.4]
%         x_lim = [1 5];
%         figure(1),subplot(4,1,1),plot(outA.ppm,outA.specs,'r','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('HERMES A: 4.56ppm')
%         figure(1),subplot(4,1,2),plot(outB.ppm,outB.specs,'b','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('B: 1.9ppm')
%         figure(1),subplot(4,1,3),plot(outC.ppm,outC.specs,'g','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('C: 4.56ppm 1.90ppm')
%         figure(1),subplot(4,1,4),plot(outD.ppm,outD.specs,'k','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('D: non-editing')
%         figure(2), plot(outA.ppm,outA.specs+outB.specs+outC.specs+outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('Sum')
%         figure(3), plot(outA.ppm,outA.specs-outB.specs+outC.specs-outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('GSH spec')
%         figure(4), plot(outA.ppm,-outA.specs+outB.specs+outC.specs-outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('GABAGlx spec')
%     case 'HERCULES'
%         x_lim = [1 5];
%         figure(1),subplot(4,1,1),plot(outA.ppm,outA.specs,'r','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('HERCULES A: 4.58ppm')
%         figure(1),subplot(4,1,2),plot(outB.ppm,outB.specs,'b','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('B: 4.18')
%         figure(1),subplot(4,1,3),plot(outC.ppm,outC.specs,'g','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('C: 4.58 ppm 1.90ppm')
%         figure(1),subplot(4,1,4),plot(outD.ppm,outD.specs,'k','linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('D: 4.18 1.9ppm')
%         figure(2), plot(outA.ppm,outA.specs+outB.specs+outC.specs+outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('Sum')
%         figure(3), plot(outA.ppm,outA.specs-outB.specs+outC.specs-outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('GSH spec')
%         figure(4), plot(outA.ppm,-outA.specs-outB.specs+outC.specs+outD.specs,'linewidth',2),set(gca,'xdir','reverse'),xlim(x_lim), xlabel('ppm'),title('GABAGlx spec')
%         
% end

%% create basis set
folder          = sim_paras_json.working_dir;
folder          = ('/Users/steve/Desktop/sim_ultrafast/Philips_MEGA_PRESS_TE68');
addMMFlag       = 1;
vendor          = vendor{1,1};
sequence        = mega_or_hadam{1,1};
localization    = localization{1,1};
editTarget      = sim_paras_json.userInput.editTarget{1:1};

% create a basis set in .mat for Osprey
[basis]         = fit_makeBasis(folder, addMMFlag, sequence, editTarget, TE, localization, vendor);

% create a basis set in .BASIS format for LCModel
switch sequence
    case 'UnEdited'
        subspec = 1;
        subspecName = {''};
    case 'MEGA'
        subspec = [1 2 3 4];
        subspecName = {'off','on','diff1','sum'};
    case 'HERMES'
        subspec = [1 2 3 4 5 6 7];
        subspecName = {'a','b','c','d','diff1','diff2','sum'};
    case 'HERCULES'
        subspec = [1 2 3 4 5 6 7];
        subspecName = {'a','b','c','d','diff1','diff2','sum'};
end

for ss = 1:length(subspec)
    outfile   = [folder '/LCModel_' vendor '_' sequence '_' localization '_' editTarget '' num2str(TE) '_' subspecName{ss} '' '.BASIS'];
    RF        = io_writelcmBASIS(basis,outfile,vendor,sequence,metablist,subspec(ss));
end

% generate plot of metabolite signal from basis set
out = fit_plotBasis(basis, 1, 1);
saveas(out,fullfile(folder,'basis_set'),'bmp');

% end of basis set creation

end % end of the function
      