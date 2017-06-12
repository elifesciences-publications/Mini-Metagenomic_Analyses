function [ genenames_with_redundancy, ko_id, kegg_cat, function_cat ] = extract_KO_gene_profile( contigs, ko_list, kegg2cat, showprogress )
% Given a list of contigs (already extracted using
% extract_structarray_entries) extract a list of genes on those contigs,
% find the corresponding KO ids and then find the abundance of KO
% categories.
%
% Input:
% contigs is a structure
% ko_list is a cell array of 2 columns with gene id and cog ids
% kegg2cat is 3 columns with cog id, category, and cat names
%
% Output:
% genelist and cogid are single column lists
% kegg_cat is a profile of ko categories. KO cat has 4 columns,
% designation, name, counts, and percentage
%
% It's possible that genes don't have corresponding KO ids
% genes could also have multiple ko ids.
%
% 2016.02.04 Brian Yu
% 2016.04.27 Brian Yu the genelist returned is now a redundant list
%                     where if a geneid maps to multiple COG ids, it will
%                     return multiple lines
% 2016.04.28 Brian Yu Adapted to KO terms

if ~exist('showprogress','var')
    showprogress = 0;
end

genelist = {};
for i = 1:length(contigs)
    genelist = [genelist; fieldnames(contigs(i).genes)];
end
numgenes = length(genelist);

% filling in KO ids (ie. Kxxxx)
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

% Filling in kegg categories
% assuming all KO ids have exactly one corresponding kegg category
% THIS ASSUMPTION IS WRONG
kegg_cat = cell(size(ko_id));
function_cat = cell(size(ko_id));
for i = 1:length(ko_id)
    if showprogress && mod(i,1000)==1
        fprintf('.');
    end
    if strcmpi(ko_id{i},'Unknown_Gene')
        kegg_cat{i} = 'Unknown_Gene';
        function_cat{i} = 'Unknown_Gene';
    else
        for j = 1:size(kegg2cat,1)
            if ~isempty(strfind(kegg2cat{j,4},ko_id{i}))
                kegg_cat{i} = kegg2cat{j,1};
                function_cat{i} = kegg2cat{j,2};
                break;
            end
        end
        % if no function found, call unknown function
        if isempty(kegg_cat{i})
            kegg_cat{i} = 'Unknown_Function';
            function_cat{i} = 'Unknown_Function';
        end
    end
end
if  showprogress
    fprintf('\n');
end

% converting kegg cat to profile
% h = tabulate(kegg_cat);
% cogcat_dict = [unique(kegg2cat(:,2)) unique(kegg2cat(:,3))];
% cogcat_dict = [cogcat_dict; {'Unknown' 'Unknown'}];
% kegg_cat = [cogcat_dict cell(size(cogcat_dict,1),2)];
% for i = 1:size(kegg_cat,1)
%     kegg_cat{i,3} = 0;
%     kegg_cat{i,4} = 0;
% end
% for i = 1:size(h,1)
%     indr = find(ismember(kegg_cat(:,1),h{i,1}));
%     kegg_cat{indr,3} = h{i,2};
%     kegg_cat{indr,4} = h{i,3};
% end

end

