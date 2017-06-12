# Mini-Metagenomic_Analyses

This repository contains the original code used to produce the results in the following publication:
http://biorxiv.org/content/early/2017/03/07/114496

## Disclaimer:

All scripts are provided for reference purposes only. There may be file paths hard coded into the scripts such that they will not work on computing clusters set up differently. The bioinformatic pipeline has undergone significant changes as a result of
1. including newer versions of bioinformatic tools,
2. incorporating new functionalities, and
3. taking advantage of more powerful compution infrastructure at Stanford.

As a result, the older version of the bioinformatics pipeline is no longer maintained and is provided as is in order for the reader to get a general idea of the bioinformatic processes. 

## Contig Assembly Fram Raw Sequencing Reads

### Required tools

Trimmomatic-0.30
Fastqc
SPAdes-3.5.0
Quast-2.3
Dnaclust
Fastx
Blast-2.2.30
Bowtie2-2.1.0
Samtools-0.1.19
Snakemake

### Process Flow

Prepare a seed files including location of all sub-sample fastq files

```
python generate_snakemake_seed.py biosample_directory output_file_name read_threshold
```

Fill in the parameters file. A sample is provided

```
parameters.txt
```

Run Snakefile_toplevel1.py to produce sub-sample contigs and combined corrected reads. The following command submits at most 20 jobs at a time to a cluster managed by Slurm. root_folder is specified as the biosample directory.

```
snakemake -j 20 -w 600 -k --config location=$root_folder --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}_%j.log" --rerun-incomplete -s Snakefile_toplevel1.py
```

Perform joint assembly using SPAdes on the large memory node separately. Then, continue with the second part of the analysis including aligning sub-sample reads back to combined contigs

'''
snakemake -j 20 -w 600 -k --config location=$root_folder --cluster "sbatch --job-name={params.name} --ntasks=1 --cpus-per-task={threads} --partition={params.partition} --mem={params.mem} -o {params.name}_%j.log" --rerun-incomplete -s Snakefile_toplevel2.py
```

The final output files are

```
super_contigs.[biosample_id].fasta
super_contigs.[biosample_id].alignment_report.txt
```

## Contig Analysis and Plotting

After assembly and alignment, all contig analysis and plotting are done using MATLAB. First, alignment based contig occurrence and p values extracted from Fisher's Exact Test are computed using

```
snakemake_result_analysis_V2.m
```

tSNE plots, functional annotation, SNP, and abundance analyses are carried out with MATLAB scripts in following file

```
cluster_supercontigs.m
```

A folder containing helper functions to explore contigs and genomes is also included

## Questions and Comments

For questions or comments, please contact [Brian Yu](http:/brianyu.org)
