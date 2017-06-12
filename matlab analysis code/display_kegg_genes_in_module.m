function [genearray, KtermLabel] = display_kegg_genes_in_module( keggData, genomeName, moduleID)
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
keggModules = keggData(1);

% Find modulePath to moduleID (must be on level d)
Alevel = fieldnames(keggModules);
for a = 1:length(Alevel)
    if exist('modulePath','var')
        break;
    else
        Blevel = fieldnames(keggModules.(Alevel{a}));
        for b = 1:length(Blevel)
            if exist('modulePath','var')
                break;
            else
                Clevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}));
                for c = 1:length(Clevel)
                    if exist('modulePath','var')
                        break;
                    else
                        Dlevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}));
                        for d = 1:length(Dlevel)
                            if strcmpi(moduleID(1),'M') && length(moduleID) == 6 && strcmpi(Dlevel{d}, moduleID)
                                modulePath{1}=Alevel{a}; modulePath{2}=Blevel{b}; modulePath{3}=Clevel{c}; modulePath{4}=Dlevel{d};
                                break;
                            elseif  length(moduleID) > 6 && sum(strfind(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).description, moduleID))
                                modulePath{1}=Alevel{a}; modulePath{2}=Blevel{b}; modulePath{3}=Clevel{c}; modulePath{4}=Dlevel{d};
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
end
if ~exist('modulePath','var')
    error(sprintf('ModuleID for %s not found.',moduleID))
else
    % Collect all the K terms
    Elevel = fieldnames(keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}).(modulePath{4})); % This is all the possible K terms
    if strcmpi(Elevel{1},'description')
        Elevel = Elevel(2:end);
    end
    % Tabulate all K term names
    KtermLabel = cell(length(Elevel),1);
    for e = 1:length(Elevel)
        KtermLabel{e} = [Elevel{e} ': ' ...
            keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}).(modulePath{4}).(Elevel{e}).description];
    end
end
% keyboard;

% Print out total number of K terms in the module
fprintf('\nThe number of K terms in module %s: %s\tis\t%d\n',modulePath{4},...
    keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}).(modulePath{4}).description,...
    length(KtermLabel));
for e = 1:length(Elevel)
    fprintf('%s: %s\n',Elevel{e},...
        keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}).(modulePath{4}).(Elevel{e}).description);
end

% For every genome bin
genearray = cell(length(KtermLabel),length(keggData));
for g = 1:length(keggData)
    % Go to the appropriate module
    tempData = keggData(g).(modulePath{1}).(modulePath{2}).(modulePath{3}).(modulePath{4});
    fprintf('\n%s\n',genomeName{g});
    % For each K term
    for e = 1:length(Elevel)
        % If a gene exist,
        if iscell(tempData.(Elevel{e}).genename)
            % Print out the K term description
            fprintf('%s: %s',Elevel{e},...
                tempData.(Elevel{e}).description);
            % Print out all the genenames associated with that K term
            for i = 1:length(tempData.(Elevel{e}).genename)
                fprintf('\t%s',tempData.(Elevel{e}).genename{i});
            end
            % Add the genelist to output array
            genearray{e,g} = tempData.(Elevel{e}).genename;
            fprintf('\n');
        end
    end
end

% Make a figure showing the number of genes for each K term
geneCnt = zeros(size(genearray));
for g = 1:length(keggData)
    for e = 1:length(Elevel)
        geneCnt(e,g) = length(genearray{e,g});
    end
end
Xlabel = cell(size(genomeName));
for i = 1:length(Xlabel)
    tmp = textscan(genomeName{i},'%s','delimiter',' '); tmp = tmp{1};
    Xlabel{i} = tmp{1};
end
fid = figure(10); 
% pcolor behaves like mesh so you need to add the extra row and col
pcolor([geneCnt zeros(size(geneCnt,1),1); zeros(1,size(geneCnt,2)+1)]); 
set(gca,'fontsize',11,'xtick',1.5:(length(Xlabel)+0.5),'xticklabel',Xlabel,...
    'ytick',1.5:(length(KtermLabel)+0.5),'yticklabel',KtermLabel);
axis([1 length(Xlabel)+1 1 length(KtermLabel)+1]); axis tight; axis ij;
colorbar;

genearray = geneCnt;

end

