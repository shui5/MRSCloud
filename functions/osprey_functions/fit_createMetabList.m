% fit_createMetabList.m
% Georg Oeltzschner, Johns Hopkins University 2019.
%
% USAGE:
% metabList = fit_createMetabList;
% 
% DESCRIPTION:
% Creates a list of metabolite basis functions that are to be included in 
% a fit.
% 
% OUTPUTS:
% metabList = structure including flags (1 = included, 0 = excluded) for
%             each metabolite included in the FID-A spin system definition,
%             plus various MM basis functions that may have been included
%             with fit_makeBasis.
%
% INPUTS:
% NONE

function metabList = fit_createMetabList(includeMetabs)

% Define the set of available metabolites
% all_mets = {'AcAc','Ace','Ala','Asc','Asp','bHB','bHG','Cit','Cr','CrCH2','EA','EtOH','GABA','GPC','GSH','Glc','Gln' ...
%     ,'Glu','Gly','HCar','H2O','Ins','Lac','NAA','NAAG','PCh','PCr','PE','Phenyl' ...
%     ,'Scyllo','Ser','Tau','Tyros','MM09','MM12','MM14','MM17','MM20' ...
%     ,'Lip09','Lip13','Lip20','MM37','MM38','MM40','MM42','MMexp','mI','sI'};

% Changed Ins and Scyllo to mI and sI and added special metabolites
all_mets = {'AcAc','Ace','AcO','Ala','Asc','Asp','bHB','bHG','Cit','Cr','CrCH2','Cystat','EA','EtOH','GABA','GPC','GSH','Glc','Gln' ...
    ,'Glu','Gly','HCar','H2O','Lac','Lys','NAA','NAAG','PCh','PCr','PE','Phenyl' ...
    ,'Ser','Tau','Thr','Tyros','Val','MM09','MM12','MM14','MM17','MM20' ...
    ,'Lip09','Lip13','Lip20','MM37','MM38','MM40','MM42','MMexp','mI','sI'};

for rr = 1:length(all_mets)
    metabList.(all_mets{rr}) = 0;
end

% Select metabolites to include in basis set depending on user input
% If 'default' or 'full' are input, fill appropriately...
% Changed Ins and Scyllo to mI and sI
if length(includeMetabs) == 1
    if strcmpi(includeMetabs{1}, 'default')
        % Define the default set
%         defaultMets = {'Asc','Asp','Cr','CrCH2' ...
%                      ,'GABA','GPC','GSH','Gln','Glu' ...
%                      ,'Lac','mI','NAA','NAAG','PCh','PCr','PE' ...
%                      ,'sI','Tau','MM09' ...
%                      ,'MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20'};
         defaultMets = {'Asc','Asp','Cr','CrCH2' ...
                     ,'GABA','GPC','GSH','Gln','Glu' ...
                     ,'Lac','NAA','NAAG','PCh','PCr','PE' ...
                     ,'Tau','MM09' ...
                     ,'MM12','MM14','MM17','MM20','Lip09','Lip13','Lip20','mI','sI'};                
                 
        for ll = 1:length(defaultMets)
            metabList.(defaultMets{ll}) = 1;
        end
    elseif strcmpi(includeMetabs{1}, 'full')
        % Define the full set
        for ll = 1:length(all_mets)
            metabList.(all_mets{ll}) = 1;
        end
    end
else
    % ... otherwise, if a list of metabolite names is provided, use it.
    for ll = 1:length(includeMetabs)
        metabList.(includeMetabs{ll}) = 1;
    end
end

end