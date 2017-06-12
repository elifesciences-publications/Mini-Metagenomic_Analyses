%% 2015.10.02 Clustering SuperContigs
%  Consolidating data for clustering supercontigs based on different
%  metrics.
%  This file is used to analyze one sample
clear; clc; close all;

%% Block 1
%  2015.09.02 Evaluating clustering method and performance
clear; close all;
% we will be using the new_header and new_seq variables
folder_root = 'some_folder_location/';
datafile = [folder_root 'super_contig_coverage_kmer_data.mat'];
contigprops = [folder_root 'super_contig_properties.mat'];
% before this file, you need to comput contig coverage, pair-wise pvalue
% from coverage, kmer count with k=4,5, 
load(datafile);
load(contigprops);
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

new_length = zeros(length(new_seq),1);
for i=1:length(new_length)
    new_length(i) = length(new_seq{i});
end

%% Block 13
%  2015.11.17 Using tsne_d to plot contig clusters based on coverage or kmer
%  count. 
%
%  just need to load datafile and contig props
levels = {'phylum'}; %'class';'order';'family'};

r=corr(contig_coverage,'type','spearman');
distmat2=distmat;
for i=1:size(distmat2,1)
    for j=1:size(distmat2,2)
        if r(i,j)<0
            distmat2(i,j)=-distmat2(i,j);
        end
    end
end

% getting tsne coordinates
D = squareform(pdist(distmat2,'spearman'));
rng(2015);
x = tsne_d(D');
% x = tsne_d(distmatk4);
for i = 1:length(levels)
    % set level colors
    label = extract_field_profile(img_contigs,levels{i},0);
    label = label(:,2);
    groupName = cell(size(x,1),1);
    tmpcluster = zeros(size(x,1),1);
    for j = 1:size(x,1)
        % contigs are in the order of new_header
        c = extract_structarray_entries(img_contigs,'contigName',new_header{j});
        groupName{j} = c.(levels{i});
        tmpcluster(j) = find(strcmp(label,c.(levels{i})));
    end
    data2plot = dimensional_reduction_overlay([],tmpcluster,'tsne_d',new_length,'fid',8,'x',x,'label',groupName);
end

%% 2016.02.10 Plot verified clusters (with annotation) using block 13 output
%  used to unify colors across samples
%  could also be used to exclude certain contigs appearing in less than 3
%  wells.
load('phylum_colors.mat')
% contains color_arr and phyla, both lists of 28 entries
grey_color = ones(1,3)*160/256;
color_arr(ismember(phyla,'Unassigned'),:) = grey_color;
tmpgroup = unique(groupName,'stable');
carr = zeros(size(tmpgroup,1),3);
for i = 1:size(tmpgroup,1)
    ind = find(ismember(phyla,tmpgroup{i}) == 1);
    carr(i,:) = color_arr(ind,:);
end
figure(9); clf; 
% gscatter(y(:,1),y(:,2),groupName,carr,'.',13);
gscatter(data2plot(:,1),data2plot(:,2),groupName,carr,'.',13);
legend('location','eastoutside');
set(gca,'fontsize',14);


%% Block 15
%  2015.11.18 Clustering
% clustering tsne clusters by dbscan and then 
% tabulating clusters cell counts

% use 2 for 2 column kmer4 and 2.6 for 2 column coverage
eps = 2.6; % eps = 2 is not bad if x has 2 columns, 2.6 for 3 columns
minpts = 5;

rng(2015);
% y = tsne_d(distmatk4',[],3);
[a, ~] = dbscan(x,eps,minpts);
contig_clusters = a;

for i = 1:length(img_contigs)
    s = strfind(new_header, img_contigs(i).imgName);
    s = find(~cellfun(@isempty, s));
    if isempty(s)
        fprintf('IMG Contig with imgName %s is not found in new_header names.\n',img_contigs(i).imgName);
    else
        img_contigs(i).groupNumber = contig_clusters(s);
    end
end
[cellcount, labels] = tabulate_cellcount(img_contigs,0,'phylum',0);
[clustdata, color_arr] = dimensional_reduction_overlay([],contig_clusters,...
    'tsne_d',new_length,'fid',3,'x',x,'label',contig_clusters);

figure(2); clf;
imagesc(cellcount);  
set(gca,'xtick',1:size(cellcount,2),'xticklabel',labels);

% Plot cluster gram type graph
X = [];
for i = [1:max(labels) -1]
    a = extract_structarray_entries(img_contigs,'groupNumber',i);
    a = [a.coverage];
    X = [X a];
end
figure(4); imagesc(X); colorbar; colormap(redbluecmap);

figure(3);

%% Block 15
%  2016.02.04 Extract genes and gene profiles from certain contig clusters
%  make sure you run up to the previous block
[~,c,~] = xlsread('IMG_3300006065\IMG Data\cog2categories.xlsx');
cog2cat = c;
fid = fopen('IMG_3300006065\IMG Data\70201.assembled.faa.COG');
c=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s'); % 10 fields
fclose(fid);
cog_list = [c{1} c{2}];

%grouplist = {1,2,3,[5,17],[7,11],8,9,10,14};

for i = 1:length(grouplist)
    
    group_num = grouplist{i};
    
    % extracting contigs
    tmpcontigs = extract_structarray_entries(img_contigs, 'groupNumber', group_num(1));
    if length(group_num)>1
        for c = 2:length(group_num)
            tmpcontigs = [tmpcontigs extract_structarray_entries(img_contigs, 'groupNumber', group_num(c))];
        end
    end
    
    % finding cog profile
    [genelist,cogid,cogcat] = extract_COG_gene_profile(tmpcontigs,cog_list,cog2cat);
    
    % plotting
    figure,
    barh(cell2mat(cogcat(:,3)));
    axis([0 max(cell2mat(cogcat(:,3)))-100 0.5 size(cogcat,1)+0.5]);
    set(gca,'ytick',1:size(cogcat,1),'yticklabel',cogcat(:,2));
    xlabel('Number of Genes');
    title(['Group ' num2str(group_num)]);
    saveas(gca,['COG Profile Group ' num2str(group_num) '.tiff']),
    
    fprintf('Completed processing group number %d\n',group_num(1));
end

%% Block 16
%  2016.02.18  Run up to the previous block, must use full dismat (same
%  length as img_contigs)
%  Making distmat linkage plot

ind =[];

for i = [1:max(labels) -1]
    for j = 1:length(img_contigs)
        if img_contigs(j).groupNumber == i;
            ind = [ind find(ismember(new_header,img_contigs(j).contigName))];
        end
    end
end
distmat_temp = distmat;
distmat_temp(distmat_temp < 1e-14) = 1e-14; % putting a limit on this
figure(8); clf; imagesc(log10(distmat_temp(ind,ind))); colorbar; colormap('gray');
clear distmat_temp;
set(gca,'fontsize',14);



%% Block 17 (2016.04.26) separating each cluster of contigs
%  then extracting group based contig and gene information
%  need to be in the biosample's folder
%  also pulling out functions using KEGG

filelocation = 'grouped_data\';

% importing geneid to cog list reference
[~,c,~] = xlsread('IMG_3300006065\IMG Data\cog2categories.xlsx');
cog2cat = c;
fid = fopen('IMG_3300006065\IMG Data\70201.assembled.faa.COG');
c=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s'); % 10 fields
fclose(fid); 
cog_list = [c{1} c{2}];

% prepare arrays 
groups = [-1 1:max([img_contigs.groupNumber])];
group_genome_size = zeros(size(groups));
group_num_genes = zeros(size(groups));

for i = 1:length(groups)
    
    gr = groups(i);
    fprintf('Starting to process group %s\n',num2str(gr));
    tmp_contigs = extract_structarray_entries(img_contigs,'groupNumber',gr);
    
    % grouped contig sequences
    seqfilename = sprintf('group_%s_sequence.fasta',num2str(gr));
    fid_seq = fopen([filelocation seqfilename],'w');
    
    % grouped protein sequences
    genefilename = sprintf('group_%s_genes.fasta',num2str(gr));
    fid_gene = fopen([filelocation genefilename],'w');
    
    for c = 1:length(tmp_contigs)
        % write contig DNA sequence to file
        fprintf(fid_seq,'>%s\n%s\n',tmp_contigs(c).imgID,tmp_contigs(c).sequence);
        group_genome_size(i) = group_genome_size(i) + length(tmp_contigs(c).sequence);
        group_num_genes(i) = group_num_genes(i) + tmp_contigs(c).geneCount;
        % store protein sequences of genes
        tmp_gnames = fieldnames(tmp_contigs(c).genes);
        for j = 1:length(tmp_gnames)
            if isempty(tmp_contigs(c).genes.(tmp_gnames{j}).seq)
                fprintf('Gene %s has no amino acid sequence. The annotation is: %s\n',tmp_gnames{j},tmp_contigs(c).genes.(tmp_gnames{j}).product{1});
            else
                fprintf(fid_gene,'>%s\n%s\n',tmp_gnames{j},tmp_contigs(c).genes.(tmp_gnames{j}).seq{1});
            end
        end
    end
    
    % grouped COG ids
    cogfilename = sprintf('group_%s_cog.txt',num2str(gr));
    fid_cog = fopen([filelocation cogfilename],'w');
    % store genes and their cog ids, once per groups
    [tmp_gnames, cog_id, ~] = extract_COG_gene_profile(tmp_contigs, cog_list, cog2cat);
    [cog_id, ind] = sort(cog_id);
    tmp_gnames = tmp_gnames(ind);
    for j = 1:length(tmp_gnames)
        fprintf(fid_cog,'%s\t%s\n',tmp_gnames{j},cog_id{j});
    end
    
    % close files
    fclose(fid_seq); fclose(fid_gene); fclose(fid_cog);
    fprintf('Completed group %s\n', num2str(gr));
end
clear fid_seq fid_gene fid_cog;
[groups' group_genome_size' group_num_genes']

% adding data to an array of structures representing groups
clear group_data;
group_data(length(groups)) = struct();
for i=1:length(groups)
    group_data(i).groupNumber = groups(i);
    group_data(i).genomeSize = group_genome_size(i);
    group_data(i).geneCount = group_num_genes(i);
end

%% Block 18 Pulling out KEGG functions for groups (2016.04.28)

filelocation = 'grouped_data\';

% improting KEGG list reference
[~,c,~] = xlsread('IMG_3300006065\IMG Data\kmodlist25664.xls');
kegg2cat = c(2:end,:); % ignore headers
fid = fopen('IMG_3300006065\IMG Data\70201.assembled.faa.KO');
c=textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s'); % 10 fields
fclose(fid); clear fid;
ko_list = [c{1} c{3}]; % the second column includes KO:Kxxxxx

% prepare arrays 
groups = [-1 1:max([img_contigs.groupNumber])];
group_genome_size = zeros(size(groups));
group_num_genes = zeros(size(groups));

for i = 1:length(groups)
    gr = groups(i);
    fprintf('Starting to process group %s\n',num2str(gr));
    tmp_contigs = extract_structarray_entries(img_contigs,'groupNumber',gr);
    group_data(i).contigs = {tmp_contigs.imgID}';
    % grouped KEGG functions
%     keggfilename = sprintf('group_%s_kegg.txt',num2str(gr));
%     fid_kegg = fopen([filelocation keggfilename],'w');
    [group_data(i).geneID, group_data(i).keggID, group_data(i).keggCat, group_data(i).keggFunc] = ...
        extract_KO_gene_profile(tmp_contigs, ko_list, kegg2cat, 0);
    % close files
    fprintf('Completed group %s\n', num2str(gr));
end

%save('grouped_genes_and_functions.mat','img_contigs','group_data','groups');

%% 2016.12.07 Exploring Kegg modules and pathways
load('manually_curated_genomes\2016.09.17_lgb_manually_curated_genome_data.mat');
load('FunctionDatabase.mat');
% KeggPathway = 
%                     EnergyMetabolism: [1x1 struct]
%       CarbohydrateAndLipidMetabolism: [1x1 struct]
%     NucleotideAndAminoAcidMetabolism: [1x1 struct]
%                  SecondaryMetabolism: [1x1 struct]
%
% Use import_Kegg_modules function

genomeNames = {
    '#40 Ignavibacteria'
    '#7 Chlorobi'
    '#2 Bacteroidetes'
    '#15 Atribacteria'
    '#19 Nitrospirae'
    '#6 Bacteroidetes'
    '#27 Thermotogae'
    '#12 Atribacteria'
    '#4 Elusimicrobia'
    '#24 Chloroflexi'
    '#3 Ignavibacteria'
    '#25 Proteobacteria'
    '#8 Fervidibacteria'
    '#16 Euryarchaeota'
    '#14 Bathyarchaeota'
    '#11 Bacteroidetes'
    '#21 Thermodesulfobacteria'
    '#20 Bacteroidetes'
    };

groupNum = [40;7;2;15;19;6;27;12;4;24;3;25;8;16;14;11;21;20];
keggModules = import_Kegg_modules();
keggdata(length(groupNum),1) = keggModules;
genomeKOterms = zeros(size(groupNum));

for g = 1:length(groupNum)
    % Compute KO terms in modules for each genome
    keggdata(g) = extract_KO_module_profile(eval(['group_' num2str(groupNum(g)) '_contigs']), ...
        ko_list, keggModules, 1);
    [~, keggID, ~, ~] = ...
        extract_KO_gene_profile(eval(['group_' num2str(groupNum(g)) '_contigs']), ko_list, kegg2cat, 0);
    genomeKOterms(g) = length(keggID);
    fprintf('.');
end
fprintf('\n');

%% 2016.11.30 Making figure of abundant Kegg Functions
%  Run "Exploring Kegg modules and pathways" first

% [geneMat, moduleLabel] = tabulate_KeggModule_geneCount(keggdata,{'Pathway_module';'Energy_metabolism';'Methane_metabolism'});
[geneMat, moduleLabel] = tabulate_KeggModule_geneCount(keggdata,{'Pathway_module';'Carbohydrate_and_lipid_metabolism';'Central_carbohydrate_metabolism'});
% [geneMat, moduleLabel] = tabulate_KeggModule_geneCount(keggdata,{'Pathway_module'});
% geneMat contains geneCnt, moduleCnt, moduleRatio

% Plotting
fid = 3; figure(fid); clf; colormap(gray); imagesc(geneMat.moduleRatio);
% imagesc(geneNumber(selectedInd,:)./repmat(keggTermList,1,size(geneNumber,2)));
set(gca,'fontsize',11,'xtick',1:length(genomeNames),'xticklabel',genomeNames,...
    'ytick',1:length(moduleLabel),'yticklabel',moduleLabel);
rotateticklabel(gca,45); colorbar;

%% 2016.11.24 Saving module gene numbers to file

fid = fopen('2016.11.24_GroupBased_keggModuleGenes.csv','w');
for i=groupNum
fprintf(fid,'\t%d',i);
end
fprintf(fid,'\n');
for i=1:length(label)
fprintf(fid,'%s',kegg2cat{ismember(kegg2cat(:,1),label{i}),2});
for j = 1:length(groupNum)
fprintf(fid,'\t%d',geneNumber(i,j));
end
fprintf(fid,'\n');
end
fclose(fid);

%% 2016.11.24 Exploring location of SNPs in which genes

snpcat = unique(cog2cat(:,3));
snpcog = cell(2,length(groupNum));

synFuncCnt = zeros(length(groupNum),length(snpcat));
asynFuncCnt = synFuncCnt;

for g = 1:length(groupNum)
    genomeID = groupNum(g);
    genomeInd = find(ismember(selected_groups,genomeID));
    % synonymous snp
    tmplist = snp_genes{1,genomeInd}; tmpcog = {};
    for i=1:length(tmplist)
        ind = ismember(cog_list(:,1),tmplist{i});
        if sum(ind) > 0
            tmp = cog_list{ind,2};
            tmpcog = [tmpcog; cog2cat{ismember(cog2cat(:,1),tmp),3}];
        end
    end
    snpcog{1,g} = tmpcog;
    % tabulate counts
    for i = 1:length(tmpcog)
        ind = ismember(snpcat,tmpcog{i});
        synFuncCnt(g,ind) = synFuncCnt(g,ind) + 1;
    end
    % asynonymous snp
    tmplist = snp_genes{2,genomeInd}; tmpcog = {};
    for i=1:length(tmplist)
        ind = ismember(cog_list(:,1),tmplist{i});
        if sum(ind) > 0
            tmp = cog_list{ind,2};
            tmpcog = [tmpcog; cog2cat{ismember(cog2cat(:,1),tmp),3}];
        end
    end
    snpcog{2,g} = tmpcog;
    % tabulate counts
    for i = 1:length(tmpcog)
        ind = ismember(snpcat,tmpcog{i});
        asynFuncCnt(g,ind) = asynFuncCnt(g,ind) + 1;
    end
    fprintf('.');
end
fprintf('\n');


%% 2016.11.30 Making grouped genome snp plots

figure(4); clf;
imagesc(synFuncCnt); colorbar;
set(gca,'fontsize',12,'xtick',1:size(synFuncCnt,2),'xticklabel',snpcat',...
    'ytick',1:size(synFuncCnt,1),'yticklabel',genomeNames);
rotateticklabel(gca,90);

figure(5); clf;
imagesc(asynFuncCnt); colorbar;
set(gca,'fontsize',12,'xtick',1:size(asynFuncCnt,2),'xticklabel',snpcat',...
    'ytick',1:size(asynFuncCnt,1),'yticklabel',genomeNames);
rotateticklabel(gca,90);

%% 2016.11.30 Graph of snp distribution
genomeInd = [2 5 7:10 13 14 17 18];
xbin = linspace(1,50,7);
figure(3); clf;
for g = 1:length(genomeInd)
    % synonymous
    unique_genes = unique(snp_genes{1,genomeInd(g)});
    unique_gene_cnt = zeros(size(unique_genes));
    for i = 1:length(snp_genes{1,genomeInd(g)})
        k = find(ismember(unique_genes,snp_genes{1,genomeInd(g)}{i}));
        unique_gene_cnt(k) = unique_gene_cnt(k) + 1;
    end
    [n1,~] = hist(unique_gene_cnt,xbin);
    % asynonymous
    unique_genes = unique(snp_genes{2,genomeInd(g)});
    unique_gene_cnt = zeros(size(unique_genes));
    for i = 1:length(snp_genes{2,genomeInd(g)})
        k = find(ismember(unique_genes,snp_genes{2,genomeInd(g)}{i}));
        unique_gene_cnt(k) = unique_gene_cnt(k) + 1;
    end
    [n2,x] = hist(unique_gene_cnt,xbin);
    subplot(length(genomeInd),1,g);
    bar(x',[n1;n2]','grouped');
    axis([-5 max(xbin)+5 0 max([n1 n2])+5]);
end

%% 2016.08.24 Making grouped genome variation plots

genome_size = [group_data.genomeSize]';
genome_size = genome_size(selected_groups+1);

fid = fopen('grouped_variants/lgb_numberofvariants.txt','r');
c = textscan(fid,'%d%s%s');
num_variants = double(c{1});
labels = c{3};

figure(1); clf; bar_color = [120,35,180]/256;
bar(num_variants./genome_size*100,0.8,'facecolor',bar_color,'edgecolor','k');
set(gca,'xtick',1:length(labels),'xticklabel',labels);
rotateticklabel(gca,90); set(gca,'fontsize',14); axis([0 13 0 1.1]);


%% 2016.08.29 Tabulating Variant Results while counting coding reagion 
%  Edited 2016.11.22 added genelists as output for sym and asym genes

%  synonymous vs non-synonymous SNPs
%  only count SNPs
%  do not need distribution anymore.
%  needs to run variable selected_groups
% clear; 
close all;clc;
quality_thresh = 999;
%selected_groups = [2 3 4 6 7 8 12 19 21 24 25 27];
selected_groups = [2 3 4 6 7 8 11 39 12 14 15 16 19 20 21 24 25 27 40];
group_names = {
    'Bacteroidetes-2'
    'Unknown-Lineage-3'
    'Unknown-Lineage-4'
    'Bacteroidetes-6'
    'Unknown-Lineage-7'
    'Fervidibacteria-8'
    'Bacteroidetes-11'
    'Unknown-Lineage-39'
    'Atribacteria-12'
    'Unknown-Lineage-14'
    'Firmicutes-15'
    'Euryarchaeota-16'
    'Nitrospirae-19'
    'Bacteroidetes-20'
    'Thermodesulfobacteria-21'
    'Chloroflexi-24'
    'Proteobacteria-25'
    'Thermotogae-27'
    'Ignavibacteria-40'
    };
group_snp = []; % will have 3 columns, total, synonymous, asynonymous, rows are groups
snp_genes = {}; % 2 columns and rows equal to group number

% load img_contigs so that you have contig name and contigID correspondence
% load('grouped_genes_and_functions.mat');
contigNames = {img_contigs.contigName}';
[~,~,X] = xlsread('IMG_3300006065\lgb_gene_coordinates_combined.xlsx');
genecoord_header = X(1,:)';
locusTag = X(2:end, ismember(genecoord_header,'Locus Tag')); % cell
startCoord = cell2mat(X(2:end, ismember(genecoord_header,'Start Coord'))); % mat
endCoord = cell2mat(X(2:end, ismember(genecoord_header,'End Coord'))); % mat
strand = X(2:end, ismember(genecoord_header,'Strand'));% cell
scaffoldID = X(2:end, ismember(genecoord_header,'Scaffold ID')); % cell

% for each group ie. genome, read in contig names
tic;
for groupnum = selected_groups
    
    % Before was using grouped_variants folder
    fid = fopen(['manually_curated_genomes\group_' num2str(groupnum) '_variants.vcf'],'r');
    C = textscan(fid,'%s','delimiter','\n'); fclose(fid); C = C{1};
    numfield = textscan(C{1},'%s'); numfield = length(numfield{1});
    % matrix D contains all the vcf information for this group that passes
    % a certain quality threshold
    % D = cell(length(C),numfield);
    D = cell(1,numfield);
    for i = 1:length(C)
        t = textscan(C{i},'%s');
        if str2double(t{1}{6}) >= quality_thresh
            D(i,:) = t{1}'; 
        end
        % assert(str2double(t{1}{6}) > quality_thresh);
    end
    
    % prepare variables and arrays
    SNPnoncoding = 0;
    SNPrna = 0;
    synonymous = 0;
    asynonymous = 0;
    
    snp_asym_genes = {};
    snp_sym_genes = {};

    % for each location, first select if it's a SNP, is it over a certain
    % quality
    for varnum = 1:length(C)
        % All the tabulated SNPs have already been quality filtered
        % The variation is a SNP
        % Using 0.5* because assuming 2 cells but one genome
        if length(D{varnum,4})==1 && length(D{varnum,5})==1
            tmpout = extract_one_snp_type( D, img_contigs, contigNames, startCoord, endCoord, ...
                locusTag, strand, scaffoldID, varnum );
            switch tmpout.snpType
                case 'not_ORF'
                    SNPnoncoding = SNPnoncoding + (tmpout.alternate + 0.5*tmpout.heterozygous);
                case 'ORF_noncoding'
                    SNPrna = SNPrna + (tmpout.alternate + 0.5*tmpout.heterozygous);
                case 'ORF_synonymous'
                    synonymous = synonymous + (tmpout.alternate + 0.5*tmpout.heterozygous);
                    snp_sym_genes = [snp_sym_genes; tmpout.geneName];
                case 'ORF_asynonymous'
                    asynonymous = asynonymous + (tmpout.alternate + 0.5*tmpout.heterozygous);
                    snp_asym_genes = [snp_asym_genes; 0.5*tmpout.geneName];
                otherwise
                    fprintf('SNP type not recognized: %s\n',tmpout.snpType);
            end
        end
    end
    fprintf('Finished Processing Group %d in %f seconds\n', groupnum, toc);
    fprintf('Non-Coding SNP:\t%d\nrRNA SNP:\t%d\nSynonymous:\t%d\nNon-synonymous:\t%d\n',...
        SNPnoncoding,SNPrna,synonymous,asynonymous)
    group_snp = [group_snp; SNPnoncoding SNPrna synonymous asynonymous]; % rows are groups
    snp_genes = [snp_genes {snp_sym_genes; snp_asym_genes}];
end    


%% 2016.09.16 Plot Cell Presence
[c,~,~] = xlsread('manually_curated_genomes\lgb_genome_stats_manualCurate.xlsx','Abundance');
presence_mat = c(2:94,2:(length(group_names)+1)); % row ordered by well, column ordered by group
figure(3); clf; imagesc(presence_mat);
set(gca,'xtick',1:length(selected_groups),'xticklabel',group_names,'fontsize',14);
rotateticklabel(gca,90); colormap(gray);

%% 2016.09.20 Plot Clustered contigs and gray for other contigs
load('manually_curated_genomes\2016.09.17_lgb_manually_curated_genome_data.mat') % group_x_contigs data2plot, new_header
load('2016.09.18_grouped_SNP_total.mat') % selected_groups
% set colors
carr = distinguishable_colors(length(selected_groups));
grey_color = ones(1,3)*160/256;
dotsize = 20; 
% Save dot groups
Y = cell(length(selected_groups),1);
incorporated_contig_names = {};
for i = 1:length(selected_groups)
    tmpcontig = eval(['group_' num2str(selected_groups(i)) '_contigs']);
    Y{i} = data2plot(ismember(new_header,{tmpcontig.contigName}'),:);
    incorporated_contig_names = [incorporated_contig_names; {tmpcontig.contigName}'];
end
unincorporated2plot = data2plot(~ismember(new_header,incorporated_contig_names),:);
figure(2); clf;
scatter(unincorporated2plot(:,1),unincorporated2plot(:,2),dotsize,grey_color,'filled'); hold on;
for i = 1:length(Y)
    y = Y{i};
    scatter(y(:,1),y(:,2),dotsize,carr(i,:),'filled');
end
legend(['Unclustered-contigs'; group_names],'location','eastoutside'); hold off;
axis([-70 70 -60 80]); set(gca,'fontsize',14);

%% 2016.09.21 Plot abundance, genomesize, variation
[c,~,~] = xlsread('manually_curated_genomes\lgb_genome_stats_manualCurate.xlsx','GenomeStats');
abundance = c(:,6);
genome_size = c(:,2);
y = sum(group_snp,2)./group_coverage;
figure(5); clf; 
loglog(genome_size,y,'k.','markersize',20);
axis([3e5 4e6 1e-5 1e-2]); set(gca,'fontsize',14);
figure(6); clf;
loglog(abundance,y,'k.','markersize',20);
set(gca,'fontsize',14);

%% 2017.06.01 Sorting out shotgun coverage
load('manually_curated_genomes\2016.09.17_lgb_manually_curated_genome_data.mat')
load('2016.09.22_grouped_SNP_total.mat')
[a,b,~]=xlsread('manually_curated_genomes\Mound_shotgun_superContig_alignmentReport.xlsx');
shotgun_contig_names = b(1:end-1);
shotgun_contig_coverage = a(1:end-1,:);
% tabulate shotgun coverage by read number
group_shotgun_coverage = zeros(size(group_coverage));
grouped_genome_size = zeros(size(group_coverage));
for i = 1:length(selected_groups)
    tmpcontig = eval(['group_' num2str(selected_groups(i)) '_contigs']);
    grouped_genome_size(i) = sum([tmpcontig.contigLength]);
    for j = 1:length(tmpcontig)
        for k = 1:length(shotgun_contig_names)
            if strcmpi(tmpcontig(j).contigName, shotgun_contig_names{k})
                group_shotgun_coverage(i) = group_shotgun_coverage(i) ...
                    + shotgun_contig_coverage(k,2);
            end
        end
    end
end
% Plot
figure(7); clf;
subplot(1,2,1), barh(lgb_genome_data(:,8));
subplot(1,2,2), barh(group_shotgun_coverage./grouped_genome_size);

%% 2017.06.01 Checking genome taxonomy level
clc;
for j = 1:length(selected_groups)
    t = eval(['group_' num2str(selected_groups(j)) '_contigs']);
    num_tot_contigs = 0;
    fprintf('\n%s\n',['group_' num2str(selected_groups(j)) '_contigs']);
    for i=1:length(t)
        % if ~(strcmpi(t(i).family,'unclassified') || strcmpi(t(i).order,'Unassigned'))
        if ~strcmpi(t(i).class,'Unassigned')
            num_tot_contigs = num_tot_contigs + 1;
            fprintf('%s\t%s\t%s\t%s\n',t(i).phylum,t(i).class,t(i).family,t(i).species)
        end
    end
    fprintf('%d\t%d\n',length(t),num_tot_contigs);
end
