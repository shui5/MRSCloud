function success = ASL_pipeline_jobman(path_code, cfgfn)
% ASL processing pipeline job manager
% yli20161118

% test
% path_code  = 'C:\Users\yli199\Documents\MATLAB\ASLcloud\asl_mricloud_yli';
% cfgfn      = 'E:\asl_mricloud_test\test_multidelay\asl_paras_pcasl_multidelay.json';

asl_paras_json  = loadjson(cfgfn);
path_data       = asl_paras_json.DIRECTORY.path_data;
path_output     = asl_paras_json.DIRECTORY.path_output;

% creat a .log file to record output info/error msg
datetimestr     = sprintf('%4d%02d%02d_%02d%02d%02d',int16(clock));
diaryFile       = [path_data,filesep,'ASLMRICloud_' datetimestr '.log'];
if ~isempty(dir([path_data,filesep,'*.log']))
    delete([path_data,filesep,'*.log']);
end
diary(diaryFile); % start the diary

success = false;
    
try % if error happens, provide a .log for user to download
    %% single-subject processing
    if ~asl_paras_json.BATCHPROC.flag
        disp('ASLMRICloud: Single dataset...');
        asl_paras_json_subj = rmfield(asl_paras_json,'BATCHPROC');
        
        % write parameters to .json file at data path
        name_tag = asl_paras_json_subj.SINGLPROC.name_asl;
        if iscell(name_tag) % in case cell array
            name_tag = [char(name_tag{1}) '_asl_all'];
        end
        asl_paras_subj  = ['asl_paras_' name_tag '.json'];
        cfgfn_subj      = [path_data,filesep,asl_paras_subj];
        savejson('',asl_paras_json_subj,cfgfn_subj);
        
        tic
        ASL_pipeline(path_code, cfgfn_subj);
        toc
        
        % zip the folder for download
        delete([path_output,filesep,'*.zip']);
        zipname = [path_output,filesep,'ASLMRICloud_Output_',name_tag,'.zip'];
        zip(zipname,'ASLMRICloud_Output_*',path_output);
        disp(['ASLMRICloud: (' name_tag ') zipped for download...']);
        
        % write .zip name into text file
        fileID = fopen([path_output,filesep,'result.txt'],'w');
        fprintf(fileID,'%s',['ASLMRICloud_Output_',name_tag,'.zip']);
        fclose(fileID);
        
        success = true;
    end
    
    
    %% batch processing
    if asl_paras_json.BATCHPROC.flag
        disp('ASLMRICloud: Batch mode...');
        
        path_zip        = asl_paras_json.BATCHPROC.path_zip;
        name_zip_asl    = asl_paras_json.BATCHPROC.name_zip_asl;
        flag_zip_m0     = asl_paras_json.BATCHPROC.flag_zip_m0;
        name_zip_m0     = asl_paras_json.BATCHPROC.name_zip_m0;
        flag_zip_t1     = asl_paras_json.BATCHPROC.flag_zip_t1;
        name_zip_t1     = asl_paras_json.BATCHPROC.name_zip_t1;
        
        % unzip the data packages to zip path
        if exist([path_zip filesep 'alldata_asl'],'dir')
            rmdir([path_zip filesep 'alldata_asl'],'s');
        end
        unzip([path_zip filesep name_zip_asl],[path_zip filesep 'alldata_asl']);
        P_asl = spm_select('FPList',[path_zip filesep 'alldata_asl'],['.*\.img$']);
        N_asl = size(P_asl,1);
        [flist_asl, ~] = ASL_batchMatchFiles(P_asl, P_asl);
        
        if flag_zip_m0
            if exist([path_zip filesep 'alldata_m0'],'dir')
                rmdir([path_zip filesep 'alldata_m0'],'s');
            end
            unzip([path_zip filesep name_zip_m0],[path_zip filesep 'alldata_m0']);
            P_m0 = spm_select('FPList',[path_zip filesep 'alldata_m0'],['.*\.img$']);
            N_m0 = size(P_m0,1);
            if N_m0 ~= N_asl
                % ErrMsg
                error('Error: The number of M0 datasets does not match that of ASL datasets.');
            end
            [flist_asl, flist_m0] = ASL_batchMatchFiles(P_asl, P_m0);
        end
        
        if flag_zip_t1
            if exist([path_zip filesep 'alldata_t1'],'dir')
                rmdir([path_zip filesep 'alldata_t1'],'s');
            end
            unzip([path_zip filesep name_zip_t1],[path_zip filesep 'alldata_t1']);
            P_t1 = spm_select('FPList',[path_zip filesep 'alldata_t1'],['.*\.zip$']);
            N_t1 = size(P_t1,1);
            if N_t1 ~= N_asl
                % ErrMsg
                error('Error: The number of T1-multiatlas datasets does not match that of ASL datasets.');
            end
            [flist_asl, flist_t1] = ASL_batchMatchFiles(P_asl, P_t1);
        end
        
        % process the data one by one
        fdrlist = cell(N_asl,1);
        for ii = 1:N_asl
            
            asl_paras_json_subj = rmfield(asl_paras_json,'BATCHPROC');
            
            asl_paras_json_subj.SINGLPROC.name_asl    = [flist_asl{ii}];
            copyfile([path_zip filesep 'alldata_asl' filesep flist_asl{ii} '.*'],path_data);
            
            asl_paras_json_subj.SINGLPROC.flag_m0     = flag_zip_m0;
            asl_paras_json_subj.SINGLPROC.name_m0     = '';
            if flag_zip_m0
                asl_paras_json_subj.SINGLPROC.name_m0 = [flist_m0{ii}];
                copyfile([path_zip filesep 'alldata_m0' filesep flist_m0{ii} '.*'],path_data);
            end
            asl_paras_json_subj.SINGLPROC.flag_t1    = flag_zip_t1;
            asl_paras_json_subj.SINGLPROC.name_t1    = '';
            if flag_zip_t1
                asl_paras_json_subj.SINGLPROC.name_t1 = [flist_t1{ii}];
                copyfile([path_zip filesep 'alldata_t1' filesep flist_t1{ii} '.zip'],path_data);
            end
            
            % write parameters to .json file at data path
            asl_paras_subj  = ['asl_paras_' flist_asl{ii} '.json'];
            cfgfn_subj      = [path_data,filesep,asl_paras_subj];
            savejson('',asl_paras_json_subj,cfgfn_subj);
            
            tic
            ASL_pipeline(path_code, cfgfn_subj);
            toc
            
            % add the output folder of the subject to the zip list
            fdrlist{ii} = ['ASLMRICloud_Output_',flist_asl{ii}];
        end
        
        % zip the folder in the list for download
        delete([path_output,filesep,'*.zip']);
        zipname = [path_output,filesep,'ASLMRICloud_Output_batch.zip'];
        zip(zipname,'ASLMRICloud_Output_*',path_output);
        disp(['ASLMRICloud: (batch) zipped for download...']);
        
        % write .zip name into text file
        fileID = fopen([path_output,filesep,'result.txt'],'w');
        fprintf(fileID,'%s','ASLMRICloud_Output_batch.zip');
        fclose(fileID);
        
        success = true;
    end
    
    
catch ME % if error happens, provide a .log for user to download
    msgErr = getReport(ME);
    fileID = fopen(diaryFile,'a');
    fprintf(fileID,'\n');
    fprintf(fileID,'ASLMRICloud: Error message...\n');
    fprintf(fileID,'%s\n',msgErr);
    fclose(fileID);
    if ~isempty(dir([path_output,filesep,'*.log'])) % delete .log in path_output if any before copy
        delete([path_output,filesep,'*.log']);
    end
    copyfile(diaryFile, path_output);
    
    % zip .log for download
    zipname = [path_output,filesep,'ASLMRICloud_Output_errmsg.zip'];
    zip(zipname,[path_output,filesep,'ASLMRICloud_',datetimestr,'.log']);
    
    % write .zip name into text file
    fileID = fopen([path_output,filesep,'result.txt'],'w');
    fprintf(fileID,'%s','ASLMRICloud_Output_errmsg.zip');
    fclose(fileID);
    
    % rethrow(ME);
end

diary off; % end the diary


