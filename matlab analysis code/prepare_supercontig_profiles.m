function prepare_supercontig_profiles(result_folder, super_contig_fasta, coverage_file, ...
    l_th, coverage_thresh)
% Analyzing snakemake supercontig results
% Adapted from analyze script V2
% 2015.07.25 Brian Yu
% 2015.12.01 Brian Yu Made into a function and only handles snakemake
%                     results, not annotations required from JGI IMG

%% all contig kmer count
%  needs variables:
% result_folder
% super_contig_fasta
% l_th

[header,seq] = fastaread([result_folder '\' super_contig_fasta]);

% prefilter for long contigs
new_header = {};
new_seq = {};
for i = 1:length(seq)
    if length(seq{i}) > l_th
        new_header = [new_header; header{i}];
        new_seq = [new_seq; seq{i}];
    end
end     

numseq = length(new_header);
fprintf('Number of sequences greater than %d bp long is %d\n',l_th,numseq);
kmer_cnt = cell(numseq,1);

% find kmer for k = 4
k = 4;
totmer = {};
for i = 1:numseq
    kmer_cnt{i} = nmercount(new_seq{i},k);
    if length(totmer) ~= 4^k
        totmer = unique([totmer; kmer_cnt{i}(:,1)]);
    end
end
totmer = sort(totmer);
% collect kmer counts
fprintf('Starting to collect kmers for k = %d\n',k);
kcount = zeros(length(totmer),numseq);
for s = 1:numseq
    temp_seq = kmer_cnt{s};
    for mer = 1:length(totmer)
        if sum(ismember(temp_seq(:,1),totmer{mer}))
            kcount(mer,s) = temp_seq{ismember(temp_seq(:,1),totmer{mer}),2};
        end
    end
    if mod(s,100)==0
        fprintf('.')
    end
end
k4count = kcount;
fprintf('\n')

% find kmer for k = 5
k = 5;
totmer = {};
for i = 1:numseq
    kmer_cnt{i} = nmercount(new_seq{i},k);
    if length(totmer) ~= 4^k
        totmer = unique([totmer; kmer_cnt{i}(:,1)]);
    end
end
totmer = sort(totmer);
% collect kmer counts
fprintf('Starting to collect kmers for k = %d\n',k);
kcount = zeros(length(totmer),numseq);
for s = 1:numseq
    temp_seq = kmer_cnt{s};
    for mer = 1:length(totmer)
        if sum(ismember(temp_seq(:,1),totmer{mer}))
            kcount(mer,s) = temp_seq{ismember(temp_seq(:,1),totmer{mer}),2};
        end
    end
    if mod(s,100)==0
        fprintf('.')
    end
end
k5count = kcount;
fprintf('\n')


% %% Plot contig length distribution
% 
% contig_length = zeros(size(header));
% for i = 1:length(header)
%     contig_length(i) = length(seq{i});
% end
% figure; set(gca,'fontsize',18);
% hist(contig_length,200); grid on;
% % axis([0 20000 0 3000]);
% % line([5000 5000],[0 40000],'linestyle','--','color','r', 'linewidth',3);

% keyboard;

%% 2015.07.12 Making contig coverage pca, Mostly importing coverage data
%  This time, don't normalize to contig length
%  Need variables
% result_folder
% coverage_file

[num,label,~] = xlsread([result_folder '\' coverage_file]);
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
contig_coverage = log2(num+1);
clear num label;

keyboard;

%% 2015.08.17 Find pairwise p value based on co-accorrance

% contig coverage matrix in in Y. each row is a contig, 
% each column is a well
% threshold is -2 for normalized coverage and 11 for absolute coverage
% cannot use the distmetric_binomial_probability.m function
% Needs variable: 
% coverage_thresh;

fprintf('Starting to calculate pairwise pvalues.\n');
Ytmp = contig_coverage(find_contig_ind(contig_names,new_header),:) > coverage_thresh;
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

%% saving variables into file called super_contig_coverage_kmer_data.mat
contig_coverage = contig_coverage(find_contig_ind(contig_names,new_header),:)';
save([result_folder '\super_contig_coverage_kmer_data.mat'], ...
    'contig_coverage', ...
    'contig_length', ...
    'coverage_thresh', ...
    'distmat', ...
    'k4count','k5count', ...
    'l_th','new_header','new_seq','numseq');


end



