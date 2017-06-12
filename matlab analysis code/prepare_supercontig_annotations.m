function prepare_supercontig_annotations(result_folder, contig_annotation_file, ...
    namefile, phylodistfile, geneproductfile, geneseqfile)
% Processing Annotated Contigs from IMG
%  2015.12.01 Brian Yu
% To run this function, must run prepare_supercontig_profiles() first.

%% First load the supercontig profiles
%  This file contains variables like 'new_header'
load([result_folder '\super_contig_coverage_kmer_data.mat']);

%% Organize information from img annotations
%  2015.11.10
%  Consolidating IMG contig names 
%  This is only working on contigs >= 10kbp
[~,~,r] = xlsread([result_folder '\' contig_annotation_file]);

% create a struct for each contig including fields like kingdom, phylum etc
numcontigs = size(r,1)-1;
clear img_contigs;
s = struct('imgID','','imgName','','contigName','','domain','','phylum','',...
    'class','','order','','family','','genus','','species','','lineageIdentity','',...
    'contigLength','','geneCount','');
% clear img_contigs;
img_contigs(numcontigs) = s;
for i = 1:numcontigs
    tmp = r{i+1,1};
    img_contigs(i).imgID = tmp(22:end);
    img_contigs(i).imgName = '';
    img_contigs(i).contigName = '';
    if ~isnan(r{i+1,8}) img_contigs(i).domain = r{i+1,8}; else img_contigs(i).domain = 'Unassigned'; end;
    if ~isnan(r{i+1,9}) img_contigs(i).phylum = r{i+1,9}; else img_contigs(i).phylum = 'Unassigned'; end;
    if ~isnan(r{i+1,10}) img_contigs(i).class = r{i+1,10}; else img_contigs(i).class = 'Unassigned'; end;
    if ~isnan(r{i+1,11}) img_contigs(i).order = r{i+1,11}; else img_contigs(i).order = 'Unassigned'; end;
    if ~isnan(r{i+1,12}) img_contigs(i).family = r{i+1,12}; else img_contigs(i).family = 'Unassigned'; end;
    if ~isnan(r{i+1,13}) img_contigs(i).genus = r{i+1,13}; else img_contigs(i).genus = 'Unassigned'; end;
    if ~isnan(r{i+1,14}) img_contigs(i).species = r{i+1,14}; else img_contigs(i).species = 'Unassigned'; end;
    img_contigs(i).lineageIdentiy = r{i+1,15}; % already an array of doubles
    img_contigs(i).contigLength = r{i+1,5};
    img_contigs(i).geneCount = r{i+1,4};
end
fprintf('Supercontig lineage assignment completed.\n');

% import contig names from img annotation
img_name_map_filename = [result_folder '\' namefile];
fid = fopen(img_name_map_filename,'r'); d = textscan(fid,'%s%s'); 
fclose(fid);

% Import gene phylodist data from img annotation
img_gene_phylo_filename = [result_folder '\' phylodistfile];
fid = fopen(img_gene_phylo_filename,'r'); gene_phylo_dist = textscan(fid,'%s%s%s%s%s','delimiter','\t'); % gene_id in field 1, gene 
fclose(fid);

% Import gene function data from img annotation
img_gene_product_filename = [result_folder '\' geneproductfile];
fid = fopen(img_gene_product_filename,'r'); gene_product = textscan(fid,'%s%s%s','delimiter','\t'); % genename productName COG/Pfam/etc
fclose(fid);

% Import gene sequence data from img annotation
img_gene_sequence_filename = [result_folder '\' geneseqfile];
gene_sequence = cell(1,2); [gene_sequence{1},gene_sequence{2}] = fastaread(img_gene_sequence_filename); 
gene_sequence{1} = gene_sequence{1}'; gene_sequence{2} = gene_sequence{2}';

fprintf('File import completed.\n');
% keyboard;

% adding in imgName, contigName, grouping information etc.
for i = 1:numcontigs
    
    % find index using imgID then extract imgName
    if ismember(img_contigs(i).imgID, d{2})
        img_contigs(i).imgName = d{1}{strcmp(d{2},img_contigs(i).imgID)};
    else
        fprintf('IMG Contig with imgID %s is not found in contig names.\n',img_contigs(i).imgID);
    end
    
    % find index in new_header that contains, also add distmat entry
    s = strfind(new_header, img_contigs(i).imgName);
    s = find(~cellfun(@isempty, s));
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
    
    % Add geneIDs as fields in .genes for each gene include other info
    % gene fields are all cells of at most 1 string. It could also be empty
    % refer to the fields as genes.(geneID).seq{1}
    % This part takes ~1-5 minutes
    genelist1 = gene_phylo_dist{1}(~cellfun(@isempty,strfind(gene_phylo_dist{1},img_contigs(i).imgID)));
    genelist2 = gene_product{1}(~cellfun(@isempty,strfind(gene_product{1},img_contigs(i).imgID)));
    genelist3 = gene_sequence{1}(~cellfun(@isempty,strfind(gene_sequence{1},img_contigs(i).imgID)));
    genelist = union(union(genelist1,genelist2),genelist3);
    if length(genelist) == img_contigs(i).geneCount
        for k = 1:length(genelist)
            img_contigs(i).genes.(genelist{k}).seq = gene_sequence{2}(strcmp(gene_sequence{1},genelist{k}));
            img_contigs(i).genes.(genelist{k}).phylo = gene_phylo_dist{5}(strcmp(gene_phylo_dist{1},genelist{k}));
            img_contigs(i).genes.(genelist{k}).product = gene_product{2}(strcmp(gene_product{1},genelist{k}));
        end
    elseif length(genelist) < img_contigs(i).geneCount
        fprintf('Sequence of some genes for contig %s was not found\n', img_contigs(i).imgID);
    else
        fprintf('Something went terribly wrong. More genes found.\n');
    end
    
    % fprintf('i is %d, s is %d\n',i,s);
    % progress tracker
    if rem(i,50)==0 
        fprintf('.'); 
    end
    
end
fprintf('\n');

% save result
save([result_folder '\super_contig_properties.mat'],...
    'img_contigs');


end

