function [ absolute_coverage, shrinked_coverage, distmat, removed_contig_ind ] = ...
    make_minimeta_distmat( coverage_file, new_header, coverage_thresh, exclude_thresh )
% Givin the file path containing coverage, produces absolute # of bp
% coverage per contig and also the distmat using fisher's exact test.
%
% 2016.02.17 Brian Yu

[num,label,~] = xlsread(coverage_file);
contig_names = label(2:end,1);
% in Y, the alignment values are normalized by length already
% Find contig lengths
contig_length = zeros(size(contig_names));
for i = 1:length(contig_names)
    c = textscan(contig_names{i},'%s','delimiter','_');
    contig_length(i) = str2double(c{1}{6});
end
% multiply by contig length to get absolute coverage profile (total number
% of basepairs
num = num .* repmat(contig_length,1,size(num,2));
Y = log2(num+1);
clear num label;

% 2015.08.17 Find pairwise p value based on co-accorrance

% contig coverage matrix in in Y. each row is a contig, 
% each column is a well
% threshold is -2 for normalized coverage and 11 for absolute coverage
% cannot use the distmetric_binomial_probability.m function

Ytmp = Y(find_contig_ind(contig_names,new_header),:) > coverage_thresh;

% exclude contigs that are not present in more than exclude_thresh chambers
removed_contig_ind = sum(Ytmp,2) < exclude_thresh;
Ytmp(removed_contig_ind, :) = [];

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

absolute_coverage = Y(find_contig_ind(contig_names,new_header),:)'; % should be number of wells by number of contigs
shrinked_coverage = absolute_coverage(:,~removed_contig_ind);

end

