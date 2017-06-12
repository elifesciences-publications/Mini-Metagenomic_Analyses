function display_kegg_genes_in_pathway( keggData, genomeName, modulePath, verbose)
% Takes in an array of structures each looks like this:
% 
%         Pathway_module: [1x1 struct]
%     Structural_complex: [1x1 struct]
%         Functional_set: [1x1 struct]
%       Signature_module: [1x1 struct]
% 
% These should be the output of extract_KO_module_profile.m
% This functions takes in one moduleID and returns all the genes found in
% that module (Gaxxxx_xxxx). The function also prints to terminal the total number of K
% terms in the moduleID, the genes that exist in the genome, their Ga code
% and the functions from annotation. Also doing this for ALL organisms
%
% Input: 
% keggData is an array of structures
% moduleID is a Mxxxxx number or a string at the D level
% genomeName is a list of genomes in the same order as keggData
%
% Output:
% genearray is a list of Gaxxxxx_xxx genes for all genome bins. Dimension is
% length of KtermLabel X number of genome bins
%
% 2016.12.09 Brian Yu Created

clc;
assert(length(keggData) == length(genomeName));
Karray = cell(length(keggData),1);

switch length(modulePath)
    
    case 1
        
        % For every genome bin
        for g = 1:length(keggData)
            % Go to the appropriate module
            tempData = keggData(g).(modulePath{1});
            fprintf('\n%s\n',genomeName{g});
            Blevel = fieldnames(tempData);
            for b = 1:length(Blevel)
                Clevel = fieldnames(tempData.(Blevel{b}));
                for c = 1:length(Clevel)
                    Dlevel = fieldnames(tempData.(Blevel{b}).(Clevel{c}));
                    for d = 1:length(Dlevel)
                        % For each module
                        if verbose
                            fprintf('%s\t%s\n',Dlevel{d},tempData.(Blevel{b}).(Clevel{c}).(Dlevel{d}).description);
                        end
                        % Collect all K terms
                        Elevel = fieldnames(tempData.(Blevel{b}).(Clevel{c}).(Dlevel{d}));
                        if strcmpi(Elevel{1},'description')
                            Elevel = Elevel(2:end);
                        end
                        for e = 1:length(Elevel)
                            % For each K term. If there is a gene
                            if iscell(tempData.(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename)
                                if verbose
                                    % Print out the K term and description
                                    fprintf('\t%s: %s\n',Elevel{e},...
                                        tempData.(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).description);
                                else
                                    % Print out the K term only
                                    fprintf('%s\n',Elevel{e});
                                end
                                % Add the genelist to output array
                                Karray{g} = [Karray{g}; tempData.(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename];
                            end
                        end
                    end
                end
            end
        end
        
    case 2
        
        % For every genome bin
        for g = 1:length(keggData)
            % Go to the appropriate module
            tempData = keggData(g).(modulePath{1}).(modulePath{2});
            fprintf('\n%s\n',genomeName{g});
            Clevel = fieldnames(tempData);
            for c = 1:length(Clevel)
                Dlevel = fieldnames(tempData.(Clevel{c}));
                for d = 1:length(Dlevel)
                    % For each module
                    if verbose
                        fprintf('%s\t%s\n',Dlevel{d},tempData.(Clevel{c}).(Dlevel{d}).description);
                    end
                    % Collect all K terms
                    Elevel = fieldnames(tempData.(Clevel{c}).(Dlevel{d}));
                    if strcmpi(Elevel{1},'description')
                        Elevel = Elevel(2:end);
                    end
                    for e = 1:length(Elevel)
                        % For each K term. If there is a gene
                        if iscell(tempData.(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename)
                            if verbose
                                % Print out the K term and description
                                fprintf('\t%s: %s\n',Elevel{e},...
                                    tempData.(Clevel{c}).(Dlevel{d}).(Elevel{e}).description);
                            else
                                % Print out the K term only
                                fprintf('%s\n',Elevel{e});
                            end
                            % Add the genelist to output array
                            Karray{g} = [Karray{g}; tempData.(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename];
                        end
                    end
                end
            end
        end
        
    case 3
        
        % For every genome bin
        for g = 1:length(keggData)
            % Go to the appropriate module
            tempData = keggData(g).(modulePath{1}).(modulePath{2}).(modulePath{3});
            Dlevel = fieldnames(tempData);
            fprintf('\n%s\n',genomeName{g});
            for d = 1:length(Dlevel)
                % For each module
                if verbose
                    fprintf('%s\t%s\n',Dlevel{d},tempData.(Dlevel{d}).description);
                end
                % Collect all K terms
                Elevel = fieldnames(tempData.(Dlevel{d}));
                if strcmpi(Elevel{1},'description')
                    Elevel = Elevel(2:end);
                end
                for e = 1:length(Elevel)
                    % For each K term. If there is a gene
                    if iscell(tempData.(Dlevel{d}).(Elevel{e}).genename)
                        if verbose
                            % Print out the K term and description
                            fprintf('\t%s: %s\n',Elevel{e},...
                                tempData.(Dlevel{d}).(Elevel{e}).description);
                        else
                            % Print out the K term only
                            fprintf('%s\n',Elevel{e});
                        end
                        % Add the genelist to output array
                        Karray{g} = [Karray{g}; tempData.(Dlevel{d}).(Elevel{e}).genename];
                    end
                end
            end
        end
        
    otherwise
        error('modulePath variable length is wrong.')
end


end

