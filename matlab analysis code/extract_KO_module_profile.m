function keggData = ...
    extract_KO_module_profile( contigs, ko_list, keggModules, showprogress )
% Given a list of contigs (already extracted using
% extract_structarray_entries) extract a list of genes on those contigs,
% find the corresponding KO ids and then find the abundance of KO
% categories.
%
% Input:
% contigs is a structure
% ko_list is a cell array of 2 columns with gene id and cog ids
% keggModules is structure containing all the module and pathway
% informations.
%
% Output:
% keggData is a structure containing gene names Gaxxx_xxx where KO exist or
% nan when it does not exist
%
% It's possible that genes don't have corresponding KO ids
% genes could also have multiple ko ids.
%
% 2016.02.04 Brian Yu
% 2016.04.27 Brian Yu the genelist returned is now a redundant list
%                     where if a geneid maps to multiple COG ids, it will
%                     return multiple lines
% 2016.04.28 Brian Yu Adapted to KO terms
% 2016.12.07 Brian Yu Adapted to find KO pathway modules. 
%                     each Kxxxxx term can match to multiple Mxxxxx terms.

if ~exist('showprogress','var')
    showprogress = 0;
end

genelist = {};
for i = 1:length(contigs)
    genelist = [genelist; fieldnames(contigs(i).genes)];
end
numgenes = length(genelist);

% filling in KO ids (ie. Kxxxx) ko_id contains all the gene KO numbers
geneid = ko_list(:,1);
kolist = cell(size(ko_list,1), 1);
for i = 1:size(ko_list,1)
    kolist{i} = ko_list{i,2}(4:end); % removes the KO: start
end
ko_id = cell(numgenes,1);
genenames_with_redundancy = cell(numgenes,1);
for i = 1:length(genelist)
    if showprogress && mod(i,1000)==1
        fprintf('.');
    end
    tmp = kolist(ismember(geneid,genelist{i}));
    if isempty(tmp)
        ko_id{i} = 'Unknown_Gene';
        genenames_with_redundancy{i} = genelist{i};
    elseif length(tmp) == 1
        ko_id{i} = tmp{1};
        genenames_with_redundancy{i} = genelist{i};
    else
        ko_id{i} = tmp{1};
        genenames_with_redundancy{i} = genelist{i};
        fprintf('\nGene %s hits multiple Kxxxx tags, taking only the first one.\n',ko_id{i});
    end
end
if  showprogress
    fprintf('\n');
end

keggData = keggModules;

% For each Kxxxxx, see if a gene exist. If the gene exists, add all gene
% names to the structure replacing gene description. Otherwise add nan
Alevel = fieldnames(keggModules);
for a = 1:length(Alevel)
    Blevel = fieldnames(keggModules.(Alevel{a}));
    for b = 1:length(Blevel)
        Clevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}));
        for c = 1:length(Clevel)
            Dlevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}));
            for d = 1:length(Dlevel)
                Elevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}));
                for e = 1:length(Elevel)
                    if ~strcmpi(Elevel{e},'description')
                        ind = ismember(ko_id,Elevel{e});
                        if sum(ind) == 0
                            keggData.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename ...
                                = nan;
                        else
                            keggData.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename ...
                                = genenames_with_redundancy(ind);
                        end
                    end
                end
            end
        end
    end
end

end

