function [cellcount, groups] = tabulate_cellcount( img_contigs, ...
    true_inferred, level, observe )
% Based on the clustering of contigs, tabulate which clusters appeared in
% which chambers. This data can then be translated back into a cell countin
% the entire chip.
%
% img_contigs: array of structures, needs to be in the new_header order
%              Need to add group number and coverage before hand
% true_inferred: flag if true use level of contigs annotated. If inferred
%                each cluster is seen as a separate sequence specific group.
% level: just a string, (not a cell of strings)
% 
% cellcount: number of chambers x number of clusters/groups of binaries
% groups: cluster names or strings 
%
% 2015.11.17 Brian Yu uses dbscan

cov_th = 6; % This is log2 coverage threshold
noise_th = 0.4; % this is a ratio

contig_coverage = [img_contigs.coverage]; % ordered same as img_contigs
clusters = [img_contigs.groupNumber]'; % a column vector
contig_length = [img_contigs.contigLength]';

if true_inferred
    % create empty matrix
    groups = extract_field_profile(img_contigs,level,0);
    groups = groups(:,2);
    % groups(strcmpi(groups,'Unassigned') | strcmpi(groups,'unclassified')) = [];
    cellcount = false(size(contig_coverage,1),length(groups));
    num_contigs = zeros(length(groups),1);
    % loop over all clusters
    for i = 1:max(clusters)
        coverage = contig_coverage(:,clusters==i);
        if observe
            figure(15); clf; imagesc(coverage);
            title(['cluster ' num2str(i)]); 
            pause(4);
        end
        tmpcount = sum(coverage > cov_th, 2) > ...
            (noise_th*size(coverage,2));
        contigs = extract_structarray_entries(img_contigs,'groupNumber',i);
        p = extract_field_profile(contigs,level,1); % including length
        [~,ind] = max(cell2mat(p(:,1)));
        cellcount(:,strcmpi(groups,p{ind,2})) = ...
            cellcount(:,strcmpi(groups,p{ind,2})) | tmpcount;
        num_contigs(strcmpi(groups,p{ind,2})) = length(contigs);
    end
else
    cellcount = false(size(contig_coverage,1),max(clusters));
    group_length = zeros(max(clusters),1);
    for i = 1:max(clusters)
        group_length(i) = sum(contig_length(clusters==i));
    end
%     groups = hist(clusters(clusters>0), 1:max(clusters));
    [num_contigs, cluster_order] = sort(group_length, 'descend');
    num_contigs = num_contigs';
    for i = 1:max(clusters)
        coverage = contig_coverage(:,clusters==i);
        if observe 
            figure(15); clf; imagesc(coverage);
            title(['cluster ' num2str(i)]); 
            pause(4);
        end
        cellcount(:,cluster_order==i) = sum(coverage > cov_th, 2) > ...
            (noise_th*size(coverage,2));
    end
    groups = cluster_order;
end

fprintf('Contig Number of Each Group is:\n');
display(num_contigs);

% use Poisson stats to compute number of cells
% lambda = -log(n/tot)
num_cells = -log(sum(~cellcount)/size(contig_coverage,1))*size(contig_coverage,1);
fprintf('\nTotal Number of Cells is: %5.0f\n',sum(num_cells));

end

