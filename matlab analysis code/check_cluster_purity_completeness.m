function [purity, completeness] = check_cluster_purity_completeness( img_contigs, level, stats_required )
% checks only for 1 level and returns 2 numbers
% purity and completeness. 
% This function can only be called on contigs
% Does not include unclassified or unassigned contigs
%
% inputs:
% img_contigs: array of contigs
% stats_required: string used for computing diversity ie. 'mean' 'median'
%
% output:
% purity: number items represented in each cluster
% completeness: number of clusters each item is represented in.
%
% 2015.11.13 Brian Yu
% 2015.11.18 Changed to not do clustering but only check for the purity and
%            completeness of the clusters already generate

item_name = extract_field_profile(img_contigs,level,0);
item_name = item_name(:,2);
item_name(strcmpi(item_name,'Unassigned') | strcmpi(item_name,'unclassified')) = [];
contig_clusters = [img_contigs.groupNumber]'; % a column vector
tmpmat = false(length(item_name),max(contig_clusters));

% tabulating purity and completeness figures for one linkage cutoff
% and level
for m = 1:max(contig_clusters)
    p = extract_field_profile(extract_structarray_entries(img_contigs,'groupNumber',m),level,0);
    for n = 1:size(p,1)
        tmpmat(strcmp(item_name,p{n,2}), m) = true;
    end
end
fprintf('Calculation for level %s completed. \n',level);

switch stats_required
    case {'mean','average'}
        purity = nanmean(nansum(tmpmat,1));
        completeness = nanmean(nansum(tmpmat,2));
    case 'median'
        purity = nanmedian(nansum(tmpmat,1));
        completeness = nanmedian(nansum(tmpmat,2));
    otherwise
        error('Unknown stats_required argument %s',stats_required);
end


end
