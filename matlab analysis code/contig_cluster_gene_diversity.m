function diversity_index = contig_cluster_gene_diversity( contigs, groupID, phylo, metric, include_none )
% Calculates the diversity index for a group of contigs by extracting all
% the genes, then look at the phylogenetic distribution at a particular
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
% 2015.11.11 Brian Yu

tax_arr = {
    'domain'
    'phylum'
    'class'
    'order'
    'family'
    'genus'
    'species'
    'strain'
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
numcontigs = length(C);
genelabel = {};
for i = 1:numcontigs
    genenames = fieldnames(C(i).genes);
    for j = 1:length(genenames)
        geneID = genenames{j};
        if isempty(C(i).genes.(geneID).phylo)
            genelabel = [genelabel; 'none'];
        else
            phylo_str = C(i).genes.(geneID).phylo{1};
            tax = textscan(phylo_str,'%s','delimiter',';');
            tax = tax{1};
            if ischar(phylo)
                genelabel = [genelabel; tax{strcmpi(tax_arr,phylo)}];
            elseif isnumeric(phylo) && (length(phylo) == 1)
                genelabel = [genelabel; tax{phylo}];
            else
                fprintf('Unknown phylo type %s',phylo);
            end
        end
    end
end

% tabulate distriubtions
[labelname,~,ic] = unique(genelabel);
abundance = hist(ic, length(labelname));

% show profile
profile = cell(length(labelname),2);
for i = 1:length(labelname)
    profile{i,1} = abundance(i);
    profile{i,2} = labelname{i};
end
% profile

if ~include_none
    abundance = abundance(~strcmp(labelname,'none'));
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

