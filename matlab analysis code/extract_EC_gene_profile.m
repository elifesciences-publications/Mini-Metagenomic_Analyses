function [ genenames_with_redundancy, ec_id ] = extract_EC_gene_profile( contigs, ec_list )
% Given a list of contigs (already extracted using
% extract_structarray_entries) extract a list of genes on those contigs,
% find the corresponding EC ids and then find the abundance of COG
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
% 2016.11.22 Brian Yu updated for EC

genelist = {};
for i = 1:length(contigs)
    genelist = [genelist; fieldnames(contigs(i).genes)];
end

% filling in cog ids
geneid = ec_list(:,1);
eclist = ec_list(:,2);
ec_id = {};
genenames_with_redundancy = {};
for i = 1:length(genelist)
    tmp = eclist(ismember(geneid,genelist{i}));
    if isempty(tmp)
        ec_id = [ec_id; 'Unknown'];
        genenames_with_redundancy = [genenames_with_redundancy; genelist{i}];
    else
        for j = 1:length(tmp)
            ec_id = [ec_id; tmp{j}];
            genenames_with_redundancy = [genenames_with_redundancy; genelist{i}];
        end
    end
end

end

