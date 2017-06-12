function keggModule = import_Kegg_modules( filename )
% Takes a text file containing all Kegg modules in A-E levels
% output keggModule is a structure containing all kegg modules all the way to K terms
% 2016.12.07 Brian Yu Created

% import KeggModule list
% This location never changes
if exist('filename','var')
    fid = fopen(filename);
else
    fid = fopen('B:\Stanford\Research-Quake\YNP Metagenome\matlab analysis code\KeggTable.keg','r');
end
keggList = textscan(fid,'%s','delimiter','\n'); keggList = keggList{1};
fclose(fid);

% Compute number of modules and create a structure
for i = 1:length(keggList)
    switch keggList{i}(1)
        case 'A'
            if length(keggList{i}) > 1
                Aleveltext = strrep(strrep(keggList{i},'A<b>',''),'</b>','');
                Aleveltext = strrep(Aleveltext,' ','_');
            end
        case 'B'
            % If the line is not empty
            if length(keggList{i}) > 1
                Bleveltext = strrep(strrep(keggList{i},'B  <b>',''),'</b>','');
                Bleveltext = strrep(Bleveltext,' ','_');
            end
        case 'C'
            tmp = textscan(keggList{i},'%s'); tmp = tmp{1};
            Cleveltext = tmp{2};
            for j = 3:length(tmp)
                Cleveltext = strcat(Cleveltext,'_',tmp{j});
            end
            for spchar = ',()'
                Cleveltext = strrep(Cleveltext,spchar,'');
            end
            Cleveltext = strrep(strrep(Cleveltext,'-','_'),'/','_');
        case 'D'
            tmp = textscan(keggList{i},'%s'); tmp = tmp{1};
            Dleveltext = tmp{2};
            Dlevel_description = tmp{3};
            for j = 4:length(tmp)
                if isempty(strfind(tmp{j},'['))
                    Dlevel_description = [Dlevel_description ' ' tmp{j}];
                else
                    break;
                end
            end
            keggModule.(Aleveltext).(Bleveltext).(Cleveltext).(Dleveltext).description = Dlevel_description;
        case 'E'
            tmp = textscan(keggList{i},'%s'); tmp = tmp{1};
            Eleveltext = tmp{2};
            Elevel_description = tmp{3};
            for j = 4:length(tmp)
                if isempty(strfind(tmp{j},'['))
                    Elevel_description = [Elevel_description ' ' tmp{j}];
                else
                    break;
                end
            end
            keggModule.(Aleveltext).(Bleveltext).(Cleveltext).(Dleveltext).(Eleveltext).description = Elevel_description;
        otherwise
            fprintf('%s\n',keggList{i});
    end
end

end

