function [ydata, color_arr, carr] = dimensional_reduction_overlay( original_matrix, clusters, func_name, contig_length, varargin )
% Plots original matrix entries as dots in a dimensionally reduced format.
% Then colors the points according to clusters. This function can be used
% to visualize effectiveness of clusters generated.
%
% original_matrix: rows are variables and columns are observations. ie.
% rows are genes or kmer profile or chamber, columns are contig or cell.
% clusters: an array the same length as the number of columns of
% original_matrix. The entries specify which cluster each item belongs to.
% func_name: function used to plot the matrix. ie. tsne, pca etc. a string
% contig_length: length of the contigs in the order of clusters
%
% optional arguments:
% 'func_arguments', func_arguments: other arguments used for the function (shoud be a string)
% 'fid', fid: figure number
% 'x', x: Mx2 coordinates to plot
% 'label', label: names of the groups to put into legend
%
% output:
% ydata: dimensionally reduced coordinates for plotting
% 
% 2015.11.16 Brian Yu 
% 2015.11.19 added contig length as an arugment so that dots can have
%            different sizes based on contig length.
%            dot size does not work right now
% 2016.02.18 Added color_arr and carr as an addition output to plot tsne independently

if ~isempty(varargin)
    if rem(length(varargin),2)
        excep_err = MException('UserValueException:ValueError','Variable varargin is not even');
        throw(excep_err);
    else
        args = reshape(varargin,2,length(varargin)/2);
        for i = 1:size(args,2)
            %keyboard;
            eval([args{1,i} '=args{2,i};']);
        end
    end
end

grey_color = ones(1,3)*160/256;
dotsize = 10*(log2(contig_length) - min(log2(contig_length))) + 15;

if ~exist('x','var')
    original_matrix = original_matrix';
    if exist('func_arguments','var')
        x = eval([func_name '(original_matrix,' func_arguments ');']);
    else
        x = eval([func_name '(original_matrix);']);
    end
end

color_arr = distinguishable_colors(max(clusters));
% populating colorarr carr
carr = zeros(length(clusters),3);
for i = 1:length(clusters)
    if clusters(i) == -1
        carr(i,:) = grey_color;
    else
        carr(i,:) = color_arr(clusters(i),:);
    end
end

if sum(clusters<=0)
    color_arr = [grey_color; color_arr]; % put grey in front because -1 goes in front
end

if exist('fid','var')
    figure(fid);
else
    figure; 
end

clf; set(gca,'fontsize',18);
if exist('label','var')
    gscatter(x(:,1),x(:,2),label,color_arr,'.',7); % default is 15
    legend('location','eastoutside');
else
    scatter(x(:,1),x(:,2),10,carr,'filled'); % default is 20
end
grid on;

ydata = x;

end

