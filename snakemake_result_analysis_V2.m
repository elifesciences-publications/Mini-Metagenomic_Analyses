% Analyzing snakemake results
% Analyze script V2
% 2015.07.25 Brian Yu
clear; clc; close all;

%% all contig kmer count
%  edited to save memory
clear; 
l_th = 10000;
result_folder = '2015.07.27_Combined_Analysis\';
[header,seq] = fastaread('super_contigs.biosampleID.fasta');

% prefilter for long contigs
new_header = {};
new_seq = {};
for i = 1:length(seq)
    if length(seq{i}) >= l_th % make sure this is >=
        new_header = [new_header; header{i}];
        new_seq = [new_seq; seq{i}];
    end
end     

numseq = length(new_header)
clear header seq;
% kmer_cnt = cell(numseq,1);

% find kmer 4
k = 4; totmer = {};
for i = 1:numseq
    kmer_cnt = nmercount(new_seq{i},k);
    if length(totmer) ~= 4^k
        totmer = unique([totmer; kmer_cnt(:,1)]);
    else
        break;
    end
end
totmer = sort(totmer);
% collect kmer counts
fprintf('Starting to collect kmers\n');
kcount4 = zeros(length(totmer),numseq);
for s = 1:numseq
    temp_seq = nmercount(new_seq{s},k);
    for mer = 1:length(totmer)
        if sum(ismember(temp_seq(:,1),totmer{mer}))
            kcount4(mer,s) = temp_seq{ismember(temp_seq(:,1),totmer{mer}),2};
        end
    end
    if mod(s,1000)==0
        fprintf('.')
    end
end
fprintf('\n')

% find kmer 5
k = 5; totmer = {};
for i = 1:numseq
    kmer_cnt = nmercount(new_seq{i},k);
    if length(totmer) ~= 4^k
        totmer = unique([totmer; kmer_cnt(:,1)]);
    else
        break;
    end
end
totmer = sort(totmer);
% collect kmer counts
fprintf('Starting to collect kmers\n');
kcount5 = zeros(length(totmer),numseq);
for s = 1:numseq
    temp_seq = nmercount(new_seq{s},k);
    for mer = 1:length(totmer)
        if sum(ismember(temp_seq(:,1),totmer{mer}))
            kcount5(mer,s) = temp_seq{ismember(temp_seq(:,1),totmer{mer}),2};
        end
    end
    if mod(s,1000)==0
        fprintf('.')
    end
end
fprintf('\n')

%% Plot contig length distribution

contig_length = zeros(size(new_header));
for i = 1:length(new_header)
    contig_length(i) = length(new_seq{i});
end
figure(7); set(gca,'fontsize',18);
[n,x] = hist(contig_length,70); 
plot(x,n,'-','linewidth',2);
set(gca,'Yscale','log','Xscale','log'); grid on;
% axis([0 20000 0 3000]);
% line([5000 5000],[0 40000],'linestyle','--','color','r', 'linewidth',3);

%% 2015.07.12 Making contig coverage pca, Mostly importing coverage data

[num,label,~] = xlsread([result_folder 'super_contigs.biosampleID.alignment_report.xlsx']);
contig_names = label(2:end,1);
% in Y, the alignment values are normalized by length already
Y = log2(num+0.0001);
clear num label;
[loading,score] = pca(Y);

%% 2015.07.12 Making contig coverage pca, Mostly importing coverage data
%  This time, don't normalize to contig length

[num,label,~] = xlsread([result_folder 'super_contigs.biosampleID.alignment_report.xlsx']);
contig_names = label(2:end,1);
% in Y, the alignment values are normalized by length already
% Find contig lengths
contig_length = zeros(size(contig_names));
for i = 1:length(contig_names)
    c = textscan(contig_names{i},'%s','delimiter','_');
    contig_length(i) = str2double(c{1}{6});
end
% multiply by contig length to get absolute coverage profile
num = num .* repmat(contig_length,1,size(num,2));
Y = log2(num+1);
clear num label;
[loading,score] = pca(Y);

%% 2015.08.17 Find pairwise p value based on co-accorrance

% contig coverage matrix in in Y. each row is a contig, 
% each column is a well
% threshold is -2 for normalized coverage and 11 for absolute coverage
coverage_thresh = 11; % after log2 transform
Ytmp = Y(find_contig_ind(contig_names,new_header),:) > coverage_thresh;
tot = size(Ytmp,2);
distmat = zeros(size(Ytmp,1)); % distance matrix (can be based on anything)
for contig1 = 1:size(Ytmp,1)
    for contig2 = 1:(contig1-1)
        % create 2x2 count matrix
        tt = sum(Ytmp(contig1,:)&Ytmp(contig2,:));
        tn = sum(Ytmp(contig1,:)&~Ytmp(contig2,:));
        nt = sum(~Ytmp(contig1,:)&Ytmp(contig2,:));
        nn = sum(~(Ytmp(contig1,:)|Ytmp(contig2,:)));
        t1 = tt+tn;
        n1 = nt+nn;
        t2 = tt+nt;
        n2 = tn+nn;
        assert(tot == t1+n1);
        assert(tot == t2+n2);
        distmat(contig1,contig2) = factorial(n1)/factorial(tot)*factorial(n2)/factorial(nn)/...
            factorial(nt)*factorial(t1)/factorial(tn)*factorial(t2)/factorial(tt);
        %(nchoosek(t1,tt)/nchoosek(tot,t2))*(nchoosek(n1,nt)/nchoosek(tot,t2));
        distmat(contig2,contig1) = distmat(contig1,contig2); 
        %nchoosek(t2,tt)*nchoosek(n2,tn)/nchoosek(tot,t1);
    end
    distmat(contig1,contig1) = factorial(sum(~Ytmp(contig1,:)))/factorial(tot)*factorial(sum(Ytmp(contig1,:)));
    if rem(contig1,100) == 0
        fprintf('.');
    end
end
fprintf('\n');
clear Ytmp; 

%% 2015.08.17 Plotting pairwise contig pvalue using imagesc
% export all the clustered contig names k then
clustered_contig_names = k.ColumnNodeNames;
% the contig names is already changed to the order to new_header
t = find_contig_ind(new_header, clustered_contig_names);
figure(2); clf;
imagesc(distmat(t,t));

%%
% These still needs to be edited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2015.11.10 Scripts to process annotated contigs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2015.11.10 Scripts to process annotated contigs
% Requires contig annotation information from JGI IMG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,r] = xlsread('IMG_3300006065\2015.12.01_scaffold_cart25217_10kb.xlsx');

% create a struct for each contig including fields like kingdom, phylum etc
numcontigs = size(r,1)-1;
s = struct('imgID','','imgName','','contigName','','domain','','phylum','',...
    'class','','order','','family','','genus','','species','','lineageIdentity','',...
    'contigLength','','geneCount','');
clear img_contigs;
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

%% 2015.11.09 Plot contig gene count vs length

img_genecount = cell2mat(r(2:end,4)); % 4th column is gene count
img_contig_length = cell2mat(r(2:end,5)); % contig length is 5th column

% Plotting gene count histogram.
figure(8); clf; set(gca,'fontsize',18);
hist(img_genecount,70);
axis([0 205 0 180]); grid on;

% Plotting gene count vs length
figure(9); clf; set(gca,'fontsize',18);
loglog(img_contig_length,img_genecount,'r.','markersize',20);
axis([8000 300000 1 300]); grid on;

%% 2015.11.09 Plot contig lineages at phylum tax levels.
%  basically show pie graph of contig lineage by phylum and then show
%  family or species makeup. Also included mean and std for lineage
%  identity based on IMG annotation pipeline.

% find unique phylum names
phylum_profile = extract_field_profile(img_contigs,'phylum',1)
emptyarray = cell(size(phylum_profile,1),1);
for i=1:size(phylum_profile,1) emptyarray{i}=' '; end
figure(10);clf; set(gca,'fontsize',14,'linewidth',5);
h = pie(cell2mat(phylum_profile(:,1)));
%h = pie(cell2mat(phylum_profile(:,1)),phylum_profile(:,2));
for i=1:length(h)
    if rem(i,2)
        set(h(i),'linewidth',3);
    end
end

%% 2015.11.09 Plot stack bar graphs for different tax levels
%  find unique class names
%  need to change 
levels2plot = {'class'; 'order'; 'family'};

for taxID = 3%1:length(levels2plot);
    
    taxlevel = levels2plot{taxID};
    tempProfiles = cell(size(phylum_profile,1)-1,1);
    tempContigNumber = cell(size(phylum_profile,1)-1,1);
    max_profile_len = 0;
    max_contig_cnt = 0;
    % populating matrices
    for phylumID = 1:(size(phylum_profile,1)-1)
        phylum_name = phylum_profile{phylumID,2};
        tempProfiles{phylumID} = extract_field_profile(extract_structarray_entries(...
            img_contigs,'phylum',phylum_name),taxlevel,1);
        tmp = extract_field_profile(extract_structarray_entries(...
            img_contigs,'phylum',phylum_name),taxlevel,0);
        tempContigNumber{phylumID} = tmp;
        max_profile_len = max(max_profile_len, size(tempProfiles{phylumID},1));
        max_contig_cnt = max(max_contig_cnt, max(cell2mat(tmp(:,1))));
    end
    % making plots
    figure(taxID); clf; set(gca,'fontsize',14); hold on;
    color_reference = colormap(jet(20));
    for phylumID = 1:(size(phylum_profile,1)-1)
        if phylumID < size(phylum_profile,1)-1
            x = [phylumID phylumID+1];
        else 
            x = [phylumID phylumID-1];
        end
        tmp1 = tempProfiles{phylumID};
        tmp2 = tempContigNumber{phylumID};
        y = [cell2mat(tmp1(:,1))'; zeros(1,length(tmp1(:,1)))]./1e3;
        z = find_color([cell2mat(tmp2(:,1)); 1; max_contig_cnt],color_reference);
        z = z(1:(end-2),:); % added 2 numbers here to get the correct range
        h = bar(x,y,0.9,'stacked');
        assert(length(h) == size(z,1));
        for i = 1:length(h)
            if strcmpi(tmp1{i,2},'unclassified')
                set(h(i),'linewidth',1,'facecolor',0.7*ones(1,3));
            elseif strcmpi(tmp1{i,2},'unassigned')
                set(h(i),'linewidth',1,'facecolor',0.4*ones(1,3));
            else
                set(h(i),'linewidth',1,'facecolor',z(i,:));
            end
        end
        [~,ind] = max(cell2mat(tmp1(:,1)));
        k = text(phylumID, sum(cell2mat(tmp1(:,1)))/1000+70, tmp1{ind,2},'fontsize',12);
        set(k,'rotation',90);
    end
    hold off;
    axis([0.5 phylumID+0.5 0 2.3e3]);
    set(gca,'xtick',1:length(phylum_profile(1:(end-1),2)),...
        'xticklabel',phylum_profile(1:(end-1),2));
    k = colorbar; % largest value is 5
    
    % This is a hack to just make this plot work. Handles colorbar label
    %k.TickLabels = {'1','13','25','37','49','61','73','85','97'};
    k.TickLabels = {'1','10','19','28','37','46','55','64','73','82','91'};
    
    ylabel('Combined Contig Length (kbp)');
    rotateticklabel(gca,90);
    
end

