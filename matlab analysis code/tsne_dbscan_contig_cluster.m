function img_contigs = tsne_dbscan_contig_cluster( contig_coverage, distmat, img_contigs, new_header, new_length, keep_ind )
% generates contig clusters 
% make_minimeta_distmat returns removed ind and not kept ind.
% 2016.02.17 Brian Yu

% new_header and new_length needs to be adjusted to match contig_coverage
% since some of the contigs may be excluded due to low presence.
removed_contigs = new_header(~keep_ind);
new_header = new_header(keep_ind);
new_length = new_length(keep_ind);
% need to remove excluded contigs from img_contigs but the contigs are in
% different orders
removed_imgcontigs_ind = false(size(img_contigs));
for i = 1:length(removed_contigs)
    for j = 1:length(img_contigs)
        if strcmpi(removed_contigs{i}, img_contigs(j).contigName);
            removed_imgcontigs_ind(j) = true;
        end
    end
end
img_contigs(removed_imgcontigs_ind) = [];


%% Block 13 !!!!!
%  2015.11.17 Using tsne_d to plot contig clusters based on coverage or kmer
%  count. Then color the contigs based on phylogenetic group annotation in
%  order to assess the goodness of clusters. Or just run Block 12
% 
%  use zscore for only kmer counts. Use pvalue distance for coverage based.
%  Need D from distmat, distmatk4 and distmatk5; also new_header
%  has nothing to do with linkage cutoffs
%
%  just need to load datafile and contig props
levels = {'phylum';'class';'order';'family'};

r=corr(contig_coverage,'type','spearman');
distmat2=distmat;
for i=1:size(distmat,1)
    for j=1:size(distmat,2)
        if r(i,j)<0
            distmat2(i,j)=-distmat2(i,j);
        end
    end
end

% getting tsne coordinates
D = squareform(pdist(distmat2,'spearman'));
rng(2015);
x = tsne_d(D');

for i = 1%:length(levels)
    % set level colors
    label = extract_field_profile(img_contigs,levels{i},0);
    label = label(:,2);
    groupName = cell(size(x,1),1);
    tmpcluster = zeros(size(x,1),1);
    for j = 1:size(x,1)
        c = extract_structarray_entries(img_contigs,'contigName',new_header{j});
%         if isempty(c)
%             keyboard
%         end
        groupName{j} = c.(levels{i});
        tmpcluster(j) = find(strcmp(label,c.(levels{i})));
    end
%     dimensional_reduction_overlay([],tmpcluster,'tsne_d',new_length,'fid',4,'x',x,'label',groupName);
    dimensional_reduction_overlay([],tmpcluster,'tsne_d',new_length,'x',x,'label',groupName);
end

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
% dimensional_reduction_overlay([],contig_clusters,...
%     'tsne_d',new_length,'fid',3,'x',x,'label',contig_clusters);
dimensional_reduction_overlay([],contig_clusters,...
    'tsne_d',new_length,'x',x,'label',contig_clusters);

% figure(2); clf;
% imagesc(cellcount); 
% set(gca,'xtick',1:size(cellcount,2),'xticklabel',labels);

% Plot cluster gram type graph
X = [];
for i = [1:max(labels) -1]
    a = extract_structarray_entries(img_contigs,'groupNumber',i);
    a = [a.coverage];
    X = [X a];
end
% figure(1); imagesc(X); colorbar;

end

