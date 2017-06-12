function prepare_supercontig_basic_structure(result_folder)
% Processing Annotated Contigs from IMG
%  2015.12.01 Brian Yu
% To run this function, must run prepare_supercontig_profiles() first.
% 2015.12.02 Changed to work without IMG annotations. This function will
%            only add basic contig information into the contig structure.

%% First load the supercontig profiles
%  This file contains variables like 'new_header'
load([result_folder '\super_contig_coverage_kmer_data.mat']);

%% Organize information, no img annotation required
%  2015.11.10
%  Consolidating IMG contig names 
%  This is only working on contigs >= 10kbp

% create a struct for each contig including fields like kingdom, phylum etc
numcontigs = length(new_header);
clear img_contigs;
s = struct('imgID','','imgName','','contigName','','domain','','phylum','',...
    'class','','order','','family','','genus','','species','','lineageIdentity','',...
    'contigLength','','geneCount','');
% clear img_contigs;
img_contigs(numcontigs) = s;

% adding in imgName, contigName, grouping information etc.
for i = 1:numcontigs
    
    s = i;
    if isempty(s)
        fprintf('IMG Contig with imgName %s is not found in new_header names.\n',img_contigs(i).imgName);
    else
        img_contigs(i).contigName = new_header{s};
    end
    
    % add kmer results, these should be precomputed and saved
    % kncount should have columns corresponding to new_header
    img_contigs(i).sequence = new_seq{s};
    img_contigs(i).kmer4 = k4count(:,s);
    img_contigs(i).kmer5 = k5count(:,s);
    % add coverage. contig_coverage corresponds to new_header
    img_contigs(i).coverage = contig_coverage(:, s);
    
end
fprintf('\n');

% save result
save([result_folder '\super_contig_properties.mat'],...
    'img_contigs');


end

