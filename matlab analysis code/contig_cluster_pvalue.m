function pvalues = contig_cluster_pvalue( linkage_cutoff, distmat, method, metric, type )
% based on a linkage and a method, clusters the distmat and then returns
% either the worst or median pvalue in each cluster
%
% linkage_cutoff: double between 0:1
% distmat; distance matrix
% method and metric are for clusterdata function
% type is the pvalue type returned
%
% output pvalues is an array of doubles
%
% 2015.11.13 Brian Yu

contig_clusters = clusterdata(distmat,'criterion','distance',...
    'cutoff',linkage_cutoff,'distance',metric,...
    'linkage',method);
pvalues = ones(max(contig_clusters),1);
for j = 1:max(contig_clusters)
    tmp_pvalues = distmat(contig_clusters == j,contig_clusters == j);
    switch type
        case {'worst','max','highest'}
            pvalues(j) = max(tmp_pvalues(:)); % for worst p_value
        case {'mean','average'}
            pvalues(j) = mean(tmp_pvalues(:));
        case {'median'}
            pvalues(j) = median(tmp_pvalues(:));
        otherwise
            error('Unkown type %s',type);
    end
end


end

