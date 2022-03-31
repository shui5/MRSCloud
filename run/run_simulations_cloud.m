function basis = run_simulations_cloud(json_input)
% addpath(genpath('/Users/steve/Desktop/sim_ultrafast'));
% "metablist": ["Cystat","EA","HCr","Lys","PE","Thr","Val"],
% "metablist": ["Ala","Asc","Asp","Cit","Cr","EtOH","GABA","GPC","GSH","Gln","Glu","Gly","Ins","Lac","NAA","NAAG","PCh","PCr","Phenyl","Scyllo","Ser","Tau","Tyros","bHB","bHG"],
% "metablist": ["Ace","AcAc","Asc","Asp","GABA","GPC","GSH","Gln","Glu","mI","Lac","NAAG","PCh","PCr","PE","sI","Ser","Tau","bHB"],
% "metablist": ["GABA","GPC","GSH","Gln","Glu","mI","Lac","NAAG","PCh","PCr","sI"],
tic
json_input      = '/Users/steve/Documents/My_Studies/MRSCloud/simMRS.json';
sim_paras_json  = loadjson(json_input);

metab_default   = sim_paras_json.private.metab_default;
metab_default_2 = sim_paras_json.private.metab_default_2;
metablist       = horzcat(metab_default,sim_paras_json.userInput.metablist);
vendor          = sim_paras_json.userInput.vendor;          % Options: GE/Philips/Siemens/Universal_Philips/Universal_Siemens
mega_or_hadam   = sim_paras_json.userInput.mega_or_hadam;   % Options: UnEdited/MEGA/HERMES/HERCULES
localization    = sim_paras_json.userInput.localization;    % Options: PRESS/sLASER
editTarget      = sim_paras_json.userInput.editTarget;      % Options: GABA/GSH/Lac/PE
TE              = sim_paras_json.userInput.TE;              % TE for UnEdited and MEGA only, HERMES and HERCULES are internally fixed at TE=80 ms
editOn          = sim_paras_json.userInput.editOn;          % For MEGA only, HERMES and HERCULES are internally fixed
editOff         = sim_paras_json.userInput.editOff;         % For MEGA only, HERMES and HERCULES are internally fixed
spatial_points  = sim_paras_json.userInput.spatial_points;  % Number of spatial points to simulate

if strcmp(mega_or_hadam, 'HERMES') || strcmp(mega_or_hadam, 'HERCULES')
    metablist       = horzcat(metab_default,metab_default_2,sim_paras_json.userInput.metablist);
end

flipAngle       = sim_paras_json.private.flipAngle;
centreFreq      = sim_paras_json.private.centreFreq;
edit_flipAngle  = sim_paras_json.private.edit_flipAngle;

work_dir        = sim_paras_json.DIRECTORY.work_dir;
save_dir        = sim_paras_json.DIRECTORY.save_dir;
outputFile      = sim_paras_json.DIRECTORY.outputFile;
delete([save_dir,'/*']);

% initialize metab parameters
% save_dir          = ('/Users/steve/Desktop/sim_ultrafast/save_dir');
% metablist = {'Ala','Asc','Asp','Cit','Cr','EA','EtOH','GABA','GPC',...
%             'GSH','Gln','Glu','Gly','H2O','Ins','Lac','NAA','NAAG','PCh',...
%             'PCr','PE','Phenyl','Ref0ppm','Scyllo','Ser','Tau','Tyros','bHB','bHG'};

% metablist       = {'GABA'};
% vendor          = {'Philips'};  %Options: Philips, Philips_universal, Siemens or GE
% mega_or_hadam   = {'HERCULES'}; %Options: UnEdited, MEGA, HERMES or HERCULES
% localization    = {'PRESS'};    %Options: PRESS or sLASER

% % for UnEdited and MEGA
% TE = 80; % GABA: 68/80, Lac: 140
% % for MEGA only
% editOn = 1.9; % GABA: 1.9, Lac: 4.1;
% editOff = 7.5;
% spatial_points  = 101;  % number of spatial points to simulate

%%%%%%%%%%%%%%%%%%%%
% % run this part for MATLAB R2013a or before
% c = parcluster('local'); % build the 'local' cluster object
% nw = c.NumWorkers;       % get the number of workers
% matlabpool (nw);         % assign the maximum number of available cores
%%%%%%%%%%%%%%%%%%%%

for iii = 1:length(metablist)
    switch(mega_or_hadam{1})
        case 'UnEdited'
            ppm_min                = [1.1];
            ppm_max                = [1.5];
        case 'MEGA'
%            ppm_min                = [1.1];
%            ppm_max                = [1.5];
            A                      = editOn;       %single-lobe pulse
            B                      = editOff;      %single-lobe pulse
            MRS_temp.editON        = num2cell([A B]);

        case 'HERMES'
            TE                     = 80;
%            ppm_min                = [2.8 2.75];
%            ppm_max                = [3.2 3.15];
            A                      = 4.56;          %single-lobe pulse
            B                      = 1.90;          %single-lobe pulse
            C                      = (4.56+1.9)/2;  %dual-lobe pulse
            D                      = 7.50;          %single-lobe pulse
            MRS_temp.editON        = num2cell([A B C D]);

        case 'HERCULES'
            TE                     = 80;
%            ppm_min                = [2.8 2.75];
%            ppm_max                = [3.2 3.15];
            A                      = 4.58;          %single-lobe pulse
            B                      = 4.18;          %single-lobe pulse
            C                      = (4.58+1.9)/2;  %dual-lobe pulse
            D                      = (4.18+1.9)/2;  %dual-lobe pulse
            MRS_temp.editON        = num2cell([A B C D]);
    end

    MRS_temp.flipAngle             = flipAngle;        % Flip angle degrees
    MRS_temp.centreFreq            = centreFreq;       % Center frequency of MR spectrum [ppm]
    MRS_temp.edit_flipAngle        = edit_flipAngle;   % Edited flip angle
    MRS_temp.nX                    = spatial_points;   % number of spatial points to simulate in x direction
    MRS_temp.nY                    = spatial_points;   % number of spatial points to simulate in y direction
    MRS_temp.nZ                    = spatial_points;   % number of spatial points to simulate in z direction
    MRS_temp.localization          = localization;
    MRS_temp.vendor                = vendor;
    MRS_temp.seq                   = mega_or_hadam;
    metab                          = metablist(iii);
    Nmetab                         = length(metab);
    MRS_temp.metab                 = metab{Nmetab}; %{ii};
    MRS_temp.Nmetab                = Nmetab; %ii;
    MRS_temp.TEs                   = num2cell(TE);
%    MRS_temp.ppm_min               = ppm_min;%(ii);
%    MRS_temp.ppm_max               = ppm_max;%(ii);
    MRS_temp.save_dir              = save_dir; % scnh
    MRS_temp                       = load_parameters(MRS_temp); % This is the function you need to edit to change the simulation parameters (not specified above this line)
    MRS_opt                        = MRS_temp; % Creating a struct variable with dimens >=1;

    %Simulate
    switch(localization{1})
        case 'PRESS'
            if strcmp(mega_or_hadam, 'HERMES') || strcmp(mega_or_hadam, 'HERCULES')
                [MRS_opt,outA, outB, outC, outD]  = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
            elseif strcmp(mega_or_hadam, 'MEGA')
                [MRS_opt,outA, outB]              = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
            else
                [MRS_opt,outA]                    = sim_signals(MRS_opt); %This function saves mat files for each sub-spectra simuation
            end
        case 'sLASER'
            if strcmp(mega_or_hadam, 'HERMES') || strcmp(mega_or_hadam, 'HERCULES')
                [MRS_opt,outA, outB, outC, outD]  = sim_signals_sLASER(MRS_opt); %This function saves mat files for each sub-spectra simuation
            elseif strcmp(mega_or_hadam, 'MEGA')
                [MRS_opt,outA, outB]              = sim_signals_sLASER(MRS_opt); %This function saves mat files for each sub-spectra simuation
            else
                [MRS_opt,outA]                    = sim_signals_sLASER(MRS_opt); %This function saves mat files for each sub-spectra simuation
            end
    end
end % for metablist

% matlabpool close % only for MATLAB R2013a or before
toc
%% create basis set

% create a basis set in .BASIS format for LCModel
addMMFlag       = 0;
vendor          = vendor{1,1};
sequence        = mega_or_hadam{1,1};
localization    = localization{1,1};
editTarget      = sim_paras_json.userInput.editTarget{1:1};

[basis]         = fit_makeBasis(save_dir, addMMFlag, sequence, editTarget, TE, localization, vendor);
basis.vendor    = sim_paras_json.userInput.vendor;
basis.seq       = sim_paras_json.userInput.mega_or_hadam;
basis.localization = sim_paras_json.userInput.localization;

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
    outfile   = [save_dir 'LCModel_' vendor '_' sequence '_' localization '_' editTarget '' num2str(TE) '_' subspecName{ss} '' '.BASIS'];
    %basis.fids = conj(BASIS_MRSCloud2.fids(:,:)); % scnh, Jan 11, 2022
    RF        = io_writelcmBASIS(basis,outfile,vendor,sequence,metablist,subspec(ss));

    % generate plot of metabolite signal from basis set
    out = fit_plotBasis(basis, ss, 1);
    saveas(out,fullfile(save_dir,['basis-set' '_' subspecName{ss}]),'mfig');
    saveas(out,fullfile(save_dir,['basis-set' '_' subspecName{ss}]),'pdf');
    close;
end

% create a basis set in .mat for Osprey
addMMFlag       = 1;
delete([save_dir,'/BASIS_*']); % Remove the previous BASIS with MMFlag off
[basis]         = fit_makeBasis(save_dir, addMMFlag, sequence, editTarget, TE, localization, vendor);

% zip outputfile
zipname = outputFile;
zip(zipname,save_dir);
end % end of the function
