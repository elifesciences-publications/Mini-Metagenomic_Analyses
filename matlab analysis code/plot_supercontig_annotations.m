function phylum_profile = plot_supercontig_annotations(result_folder, contig_annotation_file, img_contigs)
% Processing Annotated Contigs from IMG
%  2015.12.01 Brian Yu
% Must run prepare_supercontig_annotations() before running this function

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2015.11.10 Scripts to process annotated contigs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,~,r] = xlsread([result_folder '\' contig_annotation_file]);

% create a struct for each contig including fields like kingdom, phylum etc
numcontigs = size(r,1)-1;


%% 2015.11.09 Plot contig gene count vs length

img_genecount = cell2mat(r(2:end,4)); % 4th column is gene count
img_contig_length = cell2mat(r(2:end,5)); % contig length is 5th column

% Plotting gene count histogram.
figure; clf; set(gca,'fontsize',18);
hist(img_genecount,70);
% axis([0 205 0 180]); 
grid on;

% Plotting gene count vs length
figure; clf; set(gca,'fontsize',18);
loglog(img_contig_length,img_genecount,'r.','markersize',20);
% axis([8000 300000 1 300]); 
grid on;

%% 2015.11.09 Plot contig lineages at phylum tax levels.
%  basically show pie graph of contig lineage by phylum and then show
%  family or species makeup. Also included mean and std for lineage
%  identity based on IMG annotation pipeline.

% find unique phylum names
phylum_profile = extract_field_profile(img_contigs,'phylum',1);
emptyarray = cell(size(phylum_profile,1),1);
for i=1:size(phylum_profile,1) emptyarray{i}=' '; end
figure; clf; set(gca,'fontsize',14,'linewidth',5);
h = pie(cell2mat(phylum_profile(:,1)),emptyarray);
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

% remove 'Unassigned' contigs from phylum profile
adj_phylum_profile = phylum_profile;
adj_phylum_profile(strcmp(adj_phylum_profile(:,2),'Unassigned'),:) = [];

for taxID = 1:length(levels2plot);
    
    taxlevel = levels2plot{taxID};
    tempProfiles = cell(size(adj_phylum_profile,1),1);
    tempContigNumber = cell(size(adj_phylum_profile,1),1);
    max_profile_len = 0;
    max_contig_cnt = 0;
    % populating matrices
    for phylumID = 1:size(adj_phylum_profile,1)
        phylum_name = adj_phylum_profile{phylumID,2};
        tempProfiles{phylumID} = extract_field_profile(extract_structarray_entries(...
            img_contigs,'phylum',phylum_name),taxlevel,1);
        tmp = extract_field_profile(extract_structarray_entries(...
            img_contigs,'phylum',phylum_name),taxlevel,0);
        tempContigNumber{phylumID} = tmp;
        max_profile_len = max(max_profile_len, size(tempProfiles{phylumID},1));
        max_contig_cnt = max(max_contig_cnt, max(cell2mat(tmp(:,1))));
    end
    % making plots
    figure; clf; set(gca,'fontsize',14); hold on;
    color_reference = colormap(summer(20));
    for phylumID = 1:size(adj_phylum_profile,1)
        % this part of the code is to make a fake bar so that 
        % when plotting one bar it can be stacked.
        if phylumID < size(adj_phylum_profile,1)
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
    axis([0.5 phylumID+0.5 0 5e3]);
    set(gca,'xtick',1:length(adj_phylum_profile(:,2)),...
        'xticklabel',adj_phylum_profile(:,2));
    k = colorbar; % largest value is 5
    
    % This is a hack to just make this plot work. Handles colorbar label
    %k.TickLabels = {'1','13','25','37','49','61','73','85','97'};
    %k.TickLabels = {'1','10','19','28','37','46','55','64','73','82','91'};
    
    ylabel('Combined Contig Length (kbp)');
    rotateticklabel(gca,90);
    
end


end

