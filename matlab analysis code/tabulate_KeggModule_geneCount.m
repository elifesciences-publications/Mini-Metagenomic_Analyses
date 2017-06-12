function [ output_mat, module_label ] = tabulate_KeggModule_geneCount( keggData, modulePath )
% Takes in kegg module data containing gene names of each genome bin.
% providing the level to tabulate, returns number of KO terms found, ratio,
% or number of gene found based on option_string
%
% Input:
% keggData is an array of keggModules representing different genomes
% modulePath is a cell array containing up to 3 levels in Kegg Module
% 
% Output:
% output_mat is 3 matrices with size number of M terms by number 
% of genomes (ie. length of keggData)
% module_label is the y axislabel
%
% 2016.12.07 Brian Yu Created

keggModules = keggData(1);

% Find all labels based on length of modulePath
Mlabel = {}; module_label = {}; moduleCnt = 0;
switch length(modulePath)
    case 1
        Blevel = fieldnames(keggModules.(modulePath{1}));
        for b = 1:length(Blevel)
            Clevel = fieldnames(keggModules.(modulePath{1}).(Blevel{b}));
            for c = 1:length(Clevel)
                Dlevel = fieldnames(keggModules.(modulePath{1}).(Blevel{b}).(Clevel{c}));
                Mlabel = [Mlabel; Dlevel];
                moduleCnt = moduleCnt + length(Dlevel);
                for d = 1:length(Dlevel)
                    module_label = [module_label; ...
                        keggModules.(modulePath{1}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).description];
                end
            end
        end
    case 2
        Clevel = fieldnames(keggModules.(modulePath{1}).(modulePath{2}));
        for c = 1:length(Clevel)
            Dlevel = fieldnames(keggModules.(modulePath{1}).(modulePath{2}).(Clevel{c}));
            Mlabel = [Mlabel; Dlevel];
            moduleCnt = moduleCnt + length(Dlevel);
            for d = 1:length(Dlevel)
                module_label = [module_label; ...
                    keggModules.(modulePath{1}).(modulePath{2}).(Clevel{c}).(Dlevel{d}).description];
            end
        end
    case 3
        Dlevel = fieldnames(keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}));
        Mlabel = Dlevel;
        moduleCnt = moduleCnt + length(Mlabel);
        for d = 1:length(Dlevel)
            module_label = [module_label; ...
                keggModules.(modulePath{1}).(modulePath{2}).(modulePath{3}).(Dlevel{d}).description];
        end
    otherwise
        error('modulePath variable length is wrong.')
end
assert(~isempty(Mlabel));
assert(length(Mlabel) == length(module_label));
% keyboard;

% Tabulate K number counts, could use moduleCnt as an index
geneCnt = zeros(length(Mlabel),length(keggData)); kCnt = geneCnt; ktot = geneCnt;
Alevel = fieldnames(keggModules);
for a = 1:length(Alevel)
    Blevel = fieldnames(keggModules.(Alevel{a}));
    for b = 1:length(Blevel)
        Clevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}));
        for c = 1:length(Clevel)
            Dlevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}));
            for d = 1:length(Dlevel)
                if ((length(modulePath)==1 && strcmpi(Alevel{a},modulePath{1})) || ...
                        (length(modulePath)==2 && strcmpi(Alevel{a},modulePath{1}) ...
                        && strcmpi(Blevel{b},modulePath{2})) || ...
                        (length(modulePath)==3 && strcmpi(Alevel{a},modulePath{1}) ...
                        && strcmpi(Blevel{b},modulePath{2}) && strcmpi(Clevel{c},modulePath{3}))) && ...
                        sum(ismember(Mlabel,Dlevel{d}))>=1
                    moduleInd = find(ismember(Mlabel,Dlevel{d}));
                    Elevel = fieldnames(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}));
                    for g = 1:length(keggData)
                        keggModules = keggData(g); % This is a place that might generate bugs
                        geneCnt(moduleInd,g) = 0;
                        kCnt(moduleInd,g) = 0;
                        ktot(moduleInd,g) = 0;
                        % go through each E level and count K
                        for e = 1:length(Elevel)
                            if ~strcmpi(Elevel{e},'description')
                                ktot(moduleInd,g) = ktot(moduleInd,g)+1;
                                if iscell(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename)
                                    kCnt(moduleInd,g) = kCnt(moduleInd,g)+1;
                                    geneCnt(moduleInd,g) = geneCnt(moduleInd,g) +...
                                        length(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename);
                                else
                                    assert(isnan(keggModules.(Alevel{a}).(Blevel{b}).(Clevel{c}).(Dlevel{d}).(Elevel{e}).genename));
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% preparing output matrices
output_mat.geneCnt = geneCnt;
output_mat.moduleCnt = kCnt;
output_mat.moduleRatio = kCnt./ktot;
output_mat.ktot = ktot;

end

