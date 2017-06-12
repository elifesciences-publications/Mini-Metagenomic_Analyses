function diversity_index = contig_cluster_contig_diversity( contigs, groupID, phylo, metric, include_none )
% Calculates the diversity index for a group of contigs by extracting all
% the contigs, then look at the phylogenetic distribution at a particular
% level specified by phylo. Then, it calculates diversity based on a method
% specified by metric
%
% when calculating diversity, gene that are not annotated are not included
%
% contigs: array of structures
% groupID: double or string
%          if groupID <= 0 then the original contigs are used
% phylo: string or number
% metric: varies
% include_none: flag to include 'none'
%
% diversity_index: float
% 
% 2015.11.11 Brian Yu The difference is that this function uses contig
%                     phylo

tax_arr = {
    'domain'
    'phylum'
    'class'
    'order'
    'family'
    'genus'
    'species'
    };
    
if ~exist('include_none','var')
    include_none = false;
end

if groupID > 0
    C = extract_structarray_entries(contigs,'groupNumber',groupID);
else
    C = contigs;
end

% collect all gene functional abundances

% *************************
% I'm using 1 here as the last argument because I want to include the
% length of the contigs in the diversity calculations.
% **************************

if ischar(phylo)
    profile = extract_field_profile(C,phylo,1);
elseif isnumeric(phylo) && (length(phylo) == 1)
    profile = extract_filed_profile(C,tax_arr{phylo},1);
end

% profile
abundance = zeros(size(profile,1),1);
labelname = profile(:,2);
for i = 1:size(profile,1)
    abundance(i) = profile{i,1};
end

if ~include_none
    abundance = abundance(~(strcmp(labelname,'Unassigned') | strcmp(labelname,'unclassified')));
end

% keyboard;
% Compute diversity metric
if isempty(abundance)
    diversity_index = nan;
else
    switch metric
        case 'shannon' % larger means more diverse
            abundance = abundance ./ sum(abundance(:));
            diversity_index = 0 - (sum(abundance.*log2(abundance)));
        case 'simpson' % larger means more diverse
            diversity_index = 1 - sum(abundance.*(abundance-1))/(sum(abundance)*(sum(abundance)-1));
        case 'max' % larger means more divers
            abundance = abundance ./ sum(abundance(:));
            diversity_index = 1 - max(abundance);
        otherwise
            fprintf('Unknown diversity metric %s.\n',metric);
            diversity_index = nan;
    end
end

end

