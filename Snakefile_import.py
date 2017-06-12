##########################################################
# Rules for importing, trimming, and clustering fastq
##########################################################

#####################################
# rules
#####################################



# use an input file to specify and cat the fastq files into local10G
rule concatenate:
  # for each subsample, find fastq files from all sequencing runs and combine them.
  # need to be able to handle certain subsamples appearing on in some sequencing runs but not others
  # input: No Input
  # temp removes the files once it has been used
  output: temp("{subsample}/read1.{subsample}.fastq"), temp("{subsample}/read2.{subsample}.fastq")
  params: 
    name="combine_fastq", 
    partition="general", 
    mem="5000", # Don't change this
  threads: 1
  version: "1.0"
  run: 
    scratch = get_scratch(False)
    output_on_scratch = names_on_scratch(output, scratch)
    print(scratch)
    s = sample_table.ix[wildcards.subsample,:] 
    # if only one row, automatically become a column so need to handle that
    if len(s.shape) == 1:
      print('Sample '+wildcards.subsample+' exists in only one sequencing run.')
      fastq_dir = s.ix["biosample"] + "/" + s.ix["sequencingrun"] + "/" + s.ix["subsamplename"]
      shell("bash {code_dir}/snakehelper_combine_fastq.sh {fastq_dir} {scratch} {output_on_scratch[0]} {output_on_scratch[1]}")
    else:
      numseqruns = len(s.index)
      print('Sample '+wildcards.subsample+' exists in '+str(numseqruns)+' sequencing runs.')
      for i in range(numseqruns):
        fastq_dir = s.ix[i,"biosample"] + "/" + s.ix[i,"sequencingrun"] + "/" + s.ix[i,"subsamplename"]
        shell("bash {code_dir}/snakehelper_combine_fastq.sh {fastq_dir} {scratch} {output_on_scratch[0]} {output_on_scratch[1]}")
    assert(file_empty(output_on_scratch)),"Output fastq files are empty."
    assert(check_lines(output_on_scratch[0],output_on_scratch[1])),"Output fastq files have different number of lines."
    assert(check_fastq_ids(output_on_scratch[0],output_on_scratch[1])),"Output fastq ids are not in the same order."
    cp_from_scratch(output, scratch)


# You can't put temp on inputs need the temp for inputs
rule quality_trim:
  input: 
    "{subsample}/read1.{subsample}.fastq", 
    "{subsample}/read2.{subsample}.fastq"
  output: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/S1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq", 
    "{subsample}/S2.{subsample}.fastq"
  params: 
    name="quality_trim", 
    partition="long", 
    mem="30000" # Don't change this
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Performing trimmomatic trimming; note there are 4 outputs
    shell("bash {code_dir}/snakehelper_trimmingWithTrimmomatic.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir}")
    assert(check_lines(output_on_scratch[0],output_on_scratch[2])),"Paired output files have different number of lines."
    assert(check_fastq_ids(output_on_scratch[0],output_on_scratch[2])),"Paired output files have different read id numbers/orders."
    cp_from_scratch(output, scratch)


rule overrep_seq:
  # use "wild card" to specify for multiple parallel processes
  input: "{subsample}/{type}.{subsample}.fastq"
  output: "{subsample}/{type}.{subsample}.fastqc_results.txt"
  params: 
    name="fastqc", 
    partition="general", 
    mem="3000" # Don't change this
  threads: 1
  version: "1.0"
  run: 
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Quality checks
    shell("bash {code_dir}/snakehelper_fastqc.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir}")
    cp_from_scratch(output, scratch)


rule subsample_splitNcluster_reads:
  input: 
    "{subsample}/P1.{subsample}.fastq", 
    "{subsample}/P2.{subsample}.fastq"
  output: 
    temp("{subsample}/clusterR1_{split}.fastq"), 
    temp("{subsample}/clusterR2_{split}.fastq")
  params: 
    name="subsample_splitNcluster_reads", 
    partition=parameters.ix['subsample_clust_partition','entry'], 
    mem=parameters.ix["subsample_clust_memory",'entry'],
    similarity=parameters.ix['subsample_similarity','entry'],
    kmer_len=parameters.ix['kmer_len','entry']
  threads: 1
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = get_scratch(False)
    output_on_scratch = names_on_scratch(output, scratch)
    # Splitting off the particular segment of P1 and P2 files for clustering. split starts with 0
    num_files = int(depth_table.ix[wildcards.subsample,'subsample_clust_file_numbers'])
    print('Number of files for this subsample '+wildcards.subsample+' is: '+str(num_files))
    print('This file is '+wildcards.split+'th segment of the fastq file.')
    intermediate_files = [scratch+'/R1_'+wildcards.split+'.fastq', scratch+'/R2_'+wildcards.split+'.fastq']
    extract_subsample_fastq_segment(input, intermediate_files, num_files, int(wildcards.split))
    # Perform read clustering with dnaclust
    shell("bash {code_dir}/snakehelper_clusterFastqPair.sh {scratch} {intermediate_files} {output_on_scratch} {code_dir} {tool_dir} {params.similarity} {params.kmer_len}")
    assert(check_lines(output_on_scratch[0],output_on_scratch[1])),"Output files after clustering have different number of lines."
    assert(file_empty(output_on_scratch)),"At least one of the output files is empty."
    cp_from_scratch(output, scratch)


rule subsample_merge_fastq:
  input: get_split_filenames
  output: 
    temp("{subsample}/ClustPair1.{subsample}.fastq"), 
    temp("{subsample}/ClustPair2.{subsample}.fastq")
  params: 
    name="subsample_merge_clustered_fastq", 
    partition="general", 
    mem="5000" # don't change this
  threads: 1
  version: "1.0"
  run: 
    merge_fastq_files(input, output)


"""
rule fastq_to_fasta:
  # This function is not used in this file
  input: "{f}.fastq"
  output: temp("{f}.fasta")
  params: 
    name="fastq2fasta", 
    partition="general", 
    mem="1000" # Don't change this
  threads: 1
  version: "1.0"
  shell: "{tool_dir}/fastx/fastq_to_fasta -Q33 -i {input} -o {output}"
"""

