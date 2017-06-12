function [purity, completeness] = cluster_purity_completeness( img_contigs, ...
    distmat, new_header, linkage_cutoff_range, levels, method, metric, stats_required )
% Sweeps an array of linkage cutoff values and levels to generate 2 matrices
% purity and completeness. 
% This function can only be called on contigs
% Does not include unclassified or unassigned contigs
%
% inputs:
% img_contigs: array of contigs
% distmat; distance matrix
% new_header: header of original contigs corresponding to distmat
% linkage_cutoff_range: array of cutoff ranges ie. 0.2:0.5:0.45
% levels: taxonomic levels to represent (nx1 cell)
% method: ie. average or complete
% metric: ie. euclidean or spearman
% stats_required: string used for computing diversity ie. 'mean' 'median'
%
% output:
% purity: number items represented in each cluster
% completeness: number of clusters each item is represented in.
%
% 2015.11.13 Brian Yu 

if ~exist('include_none','var')
    include_none = false;
end

fprintf('Metric is %s and method is %s\n',metric,method);
numcontigs = length(img_contigs);
purity = zeros(length(linkage_cutoff_range),length(levels));
completeness = zeros(length(linkage_cutoff_range),length(levels));

for linkage_ind = 1:length(linkage_cutoff_range)
    
    % Clustering
    linkage_cutoff = linkage_cutoff_range(linkage_ind);
    contig_clusters = clusterdata(distmat,'criterion','distance',...
        'cutoff',linkage_cutoff,'distance',metric,...
        'linkage',method);
    
    % Find singleton and non-singleton groups
    [c,~,ic] = unique(contig_clusters);
    c = hist(ic,length(c));
    singleton_group = find(c == 1);
    non_singleton_group = find(c > 1);
    assert(length(c) == (length(singleton_group)+length(non_singleton_group)));
    
    % order contigs
    contig_pairwise_dist = pdist(distmat,metric);
    contig_linkage = linkage(contig_pairwise_dist,method);
    sample_order = optimalleaforder(contig_linkage,contig_pairwise_dist);
    
    fprintf('For linkage_cutoff = %0.2f The number of contig clusters is %d\n',linkage_cutoff,max(contig_clusters))
    fprintf('The number of singleton clusters is %d\n',length(singleton_group))
    fprintf('The number of non singleton clusters is %d\n',length(non_singleton_group))
    
    % adding in imgName, contigName, grouping information etc.
    for i = 1:numcontigs
        % find index in new_header that contains, also add distmat entry
        s = strfind(new_header, img_contigs(i).imgName);
        s = find(~cellfun(@isempty, s));
        if isempty(s)
            fprintf('IMG Contig with imgName %s is not found in new_header names.\n',img_contigs(i).imgName);
        else
            img_contigs(i).coverage_binomial_dist = distmat(:, s);
        end
        % add group name and graph index and sequence
        % the position of nth contig in new_header in clustergram is the
        % position at which sample_order == n
        img_contigs(i).groupNumber = contig_clusters(s);
        if length(find(sample_order == s)) == 1
            img_contigs(i).graphIndex = find(sample_order == s);
        else
            fprintf('The variable sample_order is either missing an contig or has multiple of the same contig.\n');
        end
    end
    
    %% finding cluster species matrix for each level and linkage value
    tic;
    for j = 1:length(levels)

        item_name = extract_field_profile(img_contigs,levels{j},0);
        item_name = item_name(:,2);
        tmpmat = false(length(item_name),max(contig_clusters));
        
        % tabulating purity and completeness figures for one linkage cutoff
        % and level
        for m = 1:max(contig_clusters)
            p = extract_field_profile(extract_structarray_entries(img_contigs,'groupNumber',m),levels{j},0);
            for n = 1:size(p,1)
                tmpmat(strcmp(item_name,p{n,2}), m) = true;
            end
        end
        fprintf('Calculation for level %s completed. \n',levels{j});

        switch stats_required
            case {'mean','average'}
                purity(linkage_ind,j) = mean(sum(tmpmat,1));
                completeness(linkage_ind,j) = mean(sum(tmpmat,2));
            case 'median'
                purity(linkage_ind,j) = median(sum(tmpmat,1));
                completeness(linkage_ind,j) = median(sum(tmpmat,2));
            otherwise
                error('Unknown stats_required argument %s',stats_required);
        end
        
    end
    toc;
end

end
