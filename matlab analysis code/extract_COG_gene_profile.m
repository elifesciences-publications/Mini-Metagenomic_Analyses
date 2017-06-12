function [ genenames_with_redundancy, cog_id, cog_cat ] = extract_COG_gene_profile( contigs, cog_list, cog2cat )
% Given a list of contigs (already extracted using
% extract_structarray_entries) extract a list of genes on those contigs,
% find the corresponding COG ids and then find the abundance of COG
% categories.
%
% Input:
% contigs is a structure
% cog_list is a cell array of 2 columns with gene id and cog ids
% cog2cat is 3 columns with cog id, category, and cat names
%
% Output:
% genelist and cogid are single column lists
% cogcat is a profile of cog categories. Cog cat has 4 columns,
% designation, name, counts, and percentage
%
% It's possible that genes don't have corresponding COG ids
% genes could also have multiple cog ids.
%
% 2016.02.04 Brian Yu
% 2016.04.27 Brian Yu the genelist returned is now a redundant list
%                     where if a geneid maps to multiple COG ids, it will
%                     return multiple lines

genelist = {};
for i = 1:length(contigs)
    genelist = [genelist; fieldnames(contigs(i).genes)];
end

% filling in cog ids
geneid = cog_list(:,1);
coglist = cog_list(:,2);
cog_id = {};
genenames_with_redundancy = {};
for i = 1:length(genelist)
    tmp = coglist(ismember(geneid,genelist{i}));
    if isempty(tmp)
        cog_id = [cog_id; 'Unknown'];
        genenames_with_redundancy = [genenames_with_redundancy; genelist{i}];
    else
        for j = 1:length(tmp)
            cog_id = [cog_id; tmp{j}];
            genenames_with_redundancy = [genenames_with_redundancy; genelist{i}];
        end
    end
end

% Filling in cog categories
% assuming all cog ids have exactly one corresponding cog category
cog_cat = cell(size(cog_id));
for i = 1:length(cog_id)
    if strcmpi(cog_id{i},'Unknown')
        cog_cat{i} = 'Unknown';
    else
        cog_cat{i} = cog2cat{ismember(cog2cat(:,1),cog_id{i}),2};
    end
end

% converting cog cat to profile
h = tabulate(cog_cat);
cogcat_dict = [unique(cog2cat(:,2)) unique(cog2cat(:,3))];
cogcat_dict = [cogcat_dict; {'Unknown' 'Unknown'}];
cog_cat = [cogcat_dict cell(size(cogcat_dict,1),2)];
for i = 1:size(cog_cat,1)
    cog_cat{i,3} = 0;
    cog_cat{i,4} = 0;
end
for i = 1:size(h,1)
    indr = find(ismember(cog_cat(:,1),h{i,1}));
    cog_cat{indr,3} = h{i,2};
    cog_cat{indr,4} = h{i,3};
end

end

