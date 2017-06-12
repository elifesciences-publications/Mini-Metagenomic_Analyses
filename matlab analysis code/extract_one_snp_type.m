function output = extract_one_snp_type( D, img_contigs, contigNames, startCoord, endCoord, ...
    locusTag, strand, scaffoldID, varnum )
% This function takes in one SNP identified via samtools and vcftools. then
% it goes throught the vcf profile and annotated contigs to find the type
% of the SNP and how many major, heterozygous and minor alleles are
% present. The heterozygous alleles represent chambers with more than 1
% cell of the same species.
%
% Input Variables:
% D - snp profile array, otherwise known as the output of the vcf file
% group_x_variants.vcf
% img_contigs - contains all contig information from IMG
% contigNames - contig name in "SuperContig_xxxx_NODE_xxx ..." format
% belonging to one group
% startCoord - array containing all starting positions in a gene on contigs
% endCoord - array containing all end positions for genes on contigs
% locusTag - cell array of gene names
% strand - cell array of + or -
% scaffoldID - contig names extracted from SNP .vcf file
% varnum - the index into the 
%
% Output Variables:
% output - a structure including a field describing SNP type, then fields
% describing how many homozygous, heterozygous, or alternate alleles were
% seen
% 
% Revision History
% 2016.09.14 Brian Yu Created
% 2016.12.20 Brian Yu each snp must be seen in one chamber as the dominant
%            allele

% User adjustable constants
readcov_thresh = 5;
data_start = 10; % from inspecting the .vcf file
% geneticCode = 11; % NCBI bacterial codon table is 11 but 1 is standard
geneticCode = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23];

% we want to know if it's in the coding region, if it
% causes a synonymous or non-synonymous mutation in the protein sequence.
snpLocation = str2double(D{varnum,2});
imgInd = ismember(contigNames,D{varnum,1});
imgContigName = img_contigs(imgInd).imgID;
contigEndORF = endCoord(ismember(scaffoldID,imgContigName)); % shorter array
contigStartORF = startCoord(ismember(scaffoldID,imgContigName)); % shorter array
tmpLocusTag = locusTag(ismember(scaffoldID,imgContigName));
tmpStrand = strand(ismember(scaffoldID,imgContigName));
geneInd = find(((contigStartORF-snpLocation).*(contigEndORF-snpLocation)) < 0);
allele_dist = D(varnum,data_start:end);

% Format 0/0:0,0,0:20 --> This is per SNP location
% This for loop handles all situations
% Tabulate SNP allele type
[output.homozygous, output.heterozygous, output.alternate] = ...
    tabulate_snp(allele_dist, readcov_thresh);

% if alternate allele exists in at least one chamber
if output.alternate > 0

% if the SNP is in an open reading frame
if sum(((contigStartORF-snpLocation).*(contigEndORF-snpLocation)) < 0)
    
    if isempty(geneInd)
        fprintf('SNP in zero genes\n'); % should never execute.
    elseif length(geneInd) > 1
        % Overlapping genes
        geneName = tmpLocusTag(geneInd);
        assert(~strcmp(geneName{1},geneName{2}),'Two genes are the same');
        fprintf('Overlapping Genes Detected at Contig: %s SNP Position: %d\n',imgContigName,snpLocation);
        % use the longer gene
        refProtSeqLengths = zeros(size(geneName));
        for i = 1:length(refProtSeqLengths)
            refProtSeqLengths(i) = length(img_contigs(imgInd).genes.(geneName{i}).seq{1});
        end
        geneInd = geneInd(refProtSeqLengths == max(refProtSeqLengths));
    end
    
    % for each SNP in coding reagions, get alleles and
    % construct two dna sequences for the gene.
    allele1 = D{varnum,4}; allele2 = D{varnum,5};
    geneName = tmpLocusTag{geneInd};
    output.geneName = geneName;
    dnaSeq1 = img_contigs(imgInd).sequence(contigStartORF(geneInd):contigEndORF(geneInd));
    dnaSeq2 = img_contigs(imgInd).sequence;
    assert(strcmp(dnaSeq2(snpLocation),allele1),'Sequence is wrong or allele position is off');
    dnaSeq2(snpLocation) = allele2; dnaSeq2 = dnaSeq2(contigStartORF(geneInd):contigEndORF(geneInd));
    assert(strcmp(dnaSeq1,dnaSeq2)==0,'SNP is outside of the DNA sequence')
    
    % get the protein sequences
    if isempty(img_contigs(imgInd).genes.(geneName).seq) % SNP in non-protein coding sequence
        
        fprintf('Detected contig: %s gene: %s product: %s with no protein sequence.\n',...
            imgContigName,geneName,img_contigs(imgInd).genes.(geneName).product{1})
        output.snpType = 'ORF_noncoding';
        
    else % the SNP is in a protein coding sequence
        
        refProtSeq = img_contigs(imgInd).genes.(geneName).seq{1};
        if strcmp(tmpStrand{geneInd},'+')
            % fprintf('Sequence is on positive strand\n')
            % run all the translation tables from NCBI
            protSeqScore = 100*ones(length(geneticCode),1);
            for i = 1:length(geneticCode)
                protSeq1 = nt2aa(dnaSeq1,'geneticcode',geneticCode(i)); % 11 for using NCBI bacteria codon usage
                if length(protSeq1) == length(refProtSeq)
                    protSeqScore(i) = sum(~(protSeq1==refProtSeq));
                end
            end
            trans_table_ind = find(protSeqScore == min(protSeqScore));
            trans_table_ind = trans_table_ind(1);
            % fprintf('Using Translation Table %d\n',trans_table_ind)
            protSeq1 = nt2aa(dnaSeq1,'geneticcode',geneticCode(trans_table_ind));
            protSeq2 = nt2aa(dnaSeq2,'geneticcode',geneticCode(trans_table_ind));
        else
            % fprintf('Sequence is on reverse strand\n')
            % run all the translation tables from NCBI
            protSeqScore = 100*ones(length(geneticCode),1);
            for i = 1:length(geneticCode)
                protSeq1 = nt2aa(seqrcomplement(dnaSeq1),'geneticcode',geneticCode(i)); % 11 for using NCBI bacteria codon usage
                if length(protSeq1) == length(refProtSeq)
                    protSeqScore(i) = sum(~(protSeq1==refProtSeq));
                end
            end
            trans_table_ind = find(protSeqScore == min(protSeqScore));
            trans_table_ind = trans_table_ind(1);
            % fprintf('Using Translation Table %d\n',trans_table_ind)
            protSeq1 = nt2aa(seqrcomplement(dnaSeq1),'geneticcode',geneticCode(trans_table_ind));
            protSeq2 = nt2aa(seqrcomplement(dnaSeq2),'geneticcode',geneticCode(trans_table_ind));
        end
        
        % This asserts that protSeq1 and refProtSeq does not
        % differ by more than 2 amino acids (because using
        % different translation tables.
        % 2016.09.13 removed the requirement protSeq1 be the
        % same as the reference refProtSeq
        if length(protSeq1) == length(protSeq2)% && sum(~(protSeq1==refProtSeq)) <= 2
            
            % We are going to normalize by total amount of assembled contigs over
            % coverage of readcov_thresh
            
            % do something
            
            % check if SNP is synnonymous
            if strcmp(refProtSeq,protSeq2)
                output.snpType = 'ORF_synonymous';
            else
                output.snpType = 'ORF_asynonymous';
            end
            
            % Debugging
            %[~,m] = nwalign(dnaSeq1,dnaSeq2,'alphabet','NT');
            %showalignment(m);
            %fprintf('%s\t%s\t%d\n%s\n%s\n%s\n%s\n%s\n\n',imgContigName,geneName,snpLocation,dnaSeq1,dnaSeq2,refProtSeq,protSeq1,protSeq2)
            %keyboard;
        else
            % I'm handling protein translation differently so this should
            % not be executed now.
            fprintf('Detected contig: %s gene: %s that could not be translated into JGI protein sequence.\n',imgContigName, geneName)
        end
    end
    
else % protein sequence is not in a coding region
    
    output.snpType = 'not_ORF';
    
end

else % if alternate allele is not seen in at least one chamber
    
    output.snpType = 'not confident';
    
end

end


function [major_allele, both_alleles, minor_allele] = tabulate_snp(allele_dist, readcov_thresh)
% Given a cell array of strings describing allele distribution, tabulate
% number of major, minor, or both alleles among different chambers.
%

major_allele = 0; 
both_alleles = 0;
minor_allele = 0;

% Format 0/0:0,0,0:20 --> This is per SNP location
% Tabulate SNP allele type
for subsample = 1:length(allele_dist)
    t = textscan(allele_dist{subsample},'%s%s%s','delimiter',':');
    t = [t{1};t{2};t{3}]; % These are still all strings
    % last value is the number reads supporting this position
    if str2double(t{3}) >= readcov_thresh
        prob_score = textscan(t{2},'%s%s%s','delimiter',',');
        prob_score = [str2double(prob_score{1});str2double(prob_score{2});str2double(prob_score{3})];
        % 0/0 homozygous, corresponding to quality of first number
        if prob_score(1) == 0
            major_allele = major_allele + 1;
            % 0/1 heterozygous, corresponding to quality of 2nd number
        elseif prob_score(2) == 0
            both_alleles = both_alleles + 1;
            % 1/1 homozygous for alternate, corresponding to quality of 3rd number
        elseif prob_score(3) == 0
            minor_allele = minor_allele + 1;
        else
            if strcmp(t{1},'0/0')
                major_allele = major_allele + 1;
            elseif strcmp(t{1},'0/1')
                both_alleles = both_alleles + 1;
            elseif strcmp(t{1},'1/1')
                minor_allele = minor_allele + 1;
            end
        end
    end
end


end




