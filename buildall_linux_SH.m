% compile the fmri_preprocess 
% 
% how to compile/build executable of matlab code that calls spm toolbox? 
% 
% e.g. 
% fmri_preprocess_4.m calls both spm.m and spm_jobman, you have to build it by yourself including the spm toolbox. 
% by default, 
% mcc -m fmri_preprocess_4.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\spm12' -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a fisher_r2z.m
% will generate a huge .exe , but when runnin it, there will be errors like unknown variable or function cfg_mlbatch_appcfg_master and failed to get revision .etc, 
% it's because you didn't run the spm_make_standalone.m code, you have to run spm_make_standalone.m, it makes some changes to the spm12 folder, e.g. adds the cfg_mlbatch_appcfg_master.m file, copies Contents.m to Contents.txt, ( in spm.m it explains the reason: % in deployed mode, M-files are encrypted % (even if first two lines of Contents.m "should" be preserved)  vfile = fullfile(spm('Dir'),'Contents.txt');) ,.etc. 
% so first run spm_make_standalone.m and then run your own mcc building command, actually the compiled spm.exe executable is not used, we only need the small changes to the spm folder that spm_make_standalone makes. 
% you have spm toolbox downloaded and 'installed' in your matlab folder, but you have to build it by running spm_make_standalone.m , it will build a standalone program of spm.exe that can be called from commandline. 
% 
% in summary: 1. run spm_make_standalone.m spm_make_standalone ('E:\spm\bin')
% 2. mcc -m fmri_preprocess_4.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\spm12' -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a fisher_r2z.m
% mcc -m fmri_preprocess_4.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a fisher_r2z.m
% it seems that by adding the option -C, the included toolbox will not be saved in the .exe, but in a cft archieve file,
% it's good that every time it's executed, it won't unzip the huge archive of toolbox into a temp cache folder, like in windows C:\Users\Anew\AppData\Local\Temp\Anew\mcrCache8.5 , for spm12, it means unzipping 4000+ files into temp folder before running the real matlab code. starting time is much shorter now. 
% with a .ctf file, 
% without -C option, a huge exe will be generated and unzipped to temp folder 4000+ files each time it's called. 
% but only for the first time it runs, the cft release contents of archive into a folder the same place as the executable. 
% mcc -C -m fmri_preprocess_4_for_deploy.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\spm12' -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a fisher_r2z.m 
% 3. mcc -m time_courses_outrej_4.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\spm12' -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a 'E:\Program Files\MATLAB\R2015a\toolbox\art' -a fisher_r2z.m
% 
% mcc -m time_courses_outrej_4.m -a 'E:\Program Files\MATLAB\R2015a\toolbox\NIfTI_20140122' -a 'E:\Program Files\MATLAB\R2015a\toolbox\art' -a fisher_r2z.m

% first check the machine it's running on. 
% this only works on windows
% machineid = get(com.sun.security.auth.module.NTSystem,'DomainSID');

% machinename = getmachinename();

% switch machinename
%     case {'localhost'} % local machine
%         jsonlab_folder   = '/Users/steve/Documents/My_Studies/MRSCloud/functions/jsonlab-2.0/jsonlab-2.0';
%         waveform_folder  = '/Users/steve/Documents/My_Studies/MRSCloud/pulses';
%         functions_folder = '/Users/steve/Documents/My_Studies/MRSCloud/functions';
%         addpath(genpath('/Users/steve/Documents/My_Studies/MRSCloud'));
        jsonlab_folder   = '/Users/steve/Documents/My_Studies/MRSCloud/functions/jsonlab-2.0/jsonlab-2.0';
        waveform_folder_1  = '/Users/steve/Documents/My_Studies/MRSCloud/pulses_universal';
        waveform_folder_2  = '/Users/steve/Documents/My_Studies/MRSCloud/pulses_vendors';
        functions_folder = '/Users/steve/Documents/My_Studies/MRSCloud/functions/';
        addpath(genpath('/Users/steve/Documents/My_Studies/MRSCloud'));
        
        % artfolder = '/root/code/matlab/toolbox/art';
%     otherwise
%         error('unknown machinename, can''t build');
% end

%spm_make_standalone(spm_output_folder);
main_mfile = 'run_simulations_cloud.m';
% mcc ( '-C', ...
%     '-m', main_mfile, ...
%     '-a', jsonlab_folder); 

mcc ('-C', ...
     '-m', main_mfile, ...
     '-a', jsonlab_folder, '-a', waveform_folder_1, '-a', waveform_folder_2, '-a', functions_folder);

% mcc -m io_loadRFwaveform.m -a RF_GOIA_20200506_100pts.mat -a am_sg_150_100.pta -a am_sg_150_100_200pts.pta -a am_sg_200_100.pta...
%      -a am_sg_300_100.pta -a am_sg_400_150.pta -a am_sg_475_200.pta...
%      -a dl_Philips_4_04_1_90.pta -a dl_Philips_4_18_1_90.pta -a dl_Philips_4_56_1_90.pta -a dl_Philips_4_58_1_90.pta...
%      -a dl_Siemens_4_56_1_90.pta -a dl_Siemens_4_58_1_90.pta -a dl_Siemens_GABA_LAC.pta...
%      -a gtst1203_sp.pta -a orig_refoc_mao_400_4.pta -a sg100.pta -a sg100_100_0_14ms_88hz.pta -a sg100_200pts.pta...
%      -a sl_univ_pulse.pta -a univ_eddenrefo -a univ_spreddenrex.pta);

% D = rdir(spm12folder);
% addfiles = [];
% for i = 1:length(D)
%     k = findstr(D(i).name, '.m');
%     if ~isempty(k)
%         addfiles = [addfiles ,' ', D(i).name ];
%     end
% end
% addfiles
% eval( ['mcc -m fmri_preprocess_4.m -a ', addfiles ] );
