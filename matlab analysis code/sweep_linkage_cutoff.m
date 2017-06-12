function diversity_mat = sweep_linkage_cutoff( img_contigs, distmat, new_header, linkage_cutoff_range, ...
    levels, method, metric, algorithm, identity, include_none )
% Sweeps an array of linkage cutoff values and generates a matrix for each
% linkage_cutoff. The matrix includes diversity scores for at each
% taxonomic level according to algorithm.
%
% inputs:
% distmat: distance matrix
% new_header: name of the contigs in order of distmat
% contig_props: must contain img_contigs. (CHANGED TO img_contigs)
% linkage_cutoff_range: array of cutoff ranges ie. 0.2:0.5:0.45
% levels: taxonomic levels to represent (nx1 cell)
% method: ie. average or complete
% metric: ie. euclidean or spearman
% algorithm: string used for computing diversity ie. 'shannon'
% identity: contig or genes (things to tabulate)
% include_none: true or false
%
% output:
% diversity_mat: actually a cell array of different diversity matrices
% representing each linkage cutoff value
%
% 2015.11.11 Brian Yu
% 2015.11.13 contig_props is no longer a file name. Instead it's an array
% of structures containing the contigs. and now it's called img_contigs

if ~exist('include_none','var')
    include_none = false;
end

% import files
% load(datafile);
% This file need to include the follow variables
% contig_coverage,
% contig_length,
% coverage_thresh,
% distmat,
% k4count,
% k5count,
% l_th,
% new_header,
% new_seq,
% numseq

% load(contig_props);
% % This files contains
% % img_contigs

fprintf('Metric is %s and method is %s\n',metric,method);
numcontigs = length(img_contigs);
diversity_mat = cell(size(linkage_cutoff_range));

for linkage_ind = 1:length(linkage_cutoff_range)
    
    %% Block 4
    %  trying different linkage_cutoff thresholds
    
    % figure out linkage map and where should draw cutoff
    contig_pairwise_dist = pdist(distmat,metric);
    contig_linkage = linkage(contig_pairwise_dist,method);
    
    % This is the threshold that can be changed.
    linkage_cutoff = linkage_cutoff_range(linkage_ind);
    
    %  figure out optimal sample order, and cluster organization
    %  contig_clusters have the same order as distmat and new_header
    contig_clusters = clusterdata(distmat,'criterion','distance',...
        'cutoff',linkage_cutoff,'distance',metric,...
        'linkage',method);
    
    % Find singleton and non-singleton groups
    [c,~,ic] = unique(contig_clusters);
    c = hist(ic,length(c));
    singleton_group = find(c == 1);
    non_singleton_group = find(c > 1);
    assert(length(c) == (length(singleton_group)+length(non_singleton_group)));
    
    % sample_order re_orders new_header to correspond to the clustergram
    sample_order = optimalleaforder(contig_linkage,contig_pairwise_dist);
    
    fprintf('For linkage_cutoff = %0.2f The number of contig clusters is %d\n',linkage_cutoff,max(contig_clusters))
    fprintf('The number of singleton clusters is %d\n',length(singleton_group))
    fprintf('The number of non singleton clusters is %d\n',length(non_singleton_group))
    
    %% Block 6
    %  2015.11.10
    %  clustered group numbers
    %  This is only working on contigs >= 10kbp
    
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
    
    
    %% Block 7
    %  2015.11.11 Plotting abundance index of each group
    
    cluster_diversity_index = ones(max(contig_clusters),length(levels));
    
    for j = 1:length(levels)
        for i = 1:max(contig_clusters)
            switch identity
                % gene diversity by identity
                case {'gene','Gene','genes','Genes'}
                    cluster_diversity_index(i,j) = ...
                        contig_cluster_gene_diversity(img_contigs,i,levels{j},...
                        algorithm,include_none); % 1 means include 'none'
                % contig diversity by annotation including lengths
                case {'contig','contigs','Contig','Contigs'}
                    cluster_diversity_index(i,j) = ...
                        contig_cluster_contig_diversity(img_contigs,i,levels{j},...
                        algorithm,include_none);
                otherwise
                    fprintf('Unknown tabulation identity %s.\n',identity);
            end
        end
        fprintf('Calculation for level %s completed. \n',levels{j});
    end
    
    diversity_mat{linkage_ind} = cluster_diversity_index;
    
end

end
