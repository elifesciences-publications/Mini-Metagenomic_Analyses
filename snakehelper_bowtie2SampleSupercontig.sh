#!/bin/bash

###################################
# Using Bowtie2 to align unclustered
# reads from one subsample to the 
# final super contigs obtain through
# a variety of methods. 
# The output is a vector containing
# the coverage depth normalized by 
# contig length.
###################################

scratch=$1
read1=$2
read2=$3
contigs=$4
report=$5
tool_dir=$6
code_dir=$7
length_thresh=$8   # this is an integer number 
thread=$9

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
cp $code_dir/threshold_scaffolds.py $scratch

# filter out scaffolds shorter than contig_thresh                                                                                     
python threshold_scaffolds.py $length_thresh $contigs temp.fasta

# Build Bowtie2 index
samfile="alignResults"
$tool_dir"bowtie2-2.1.0/bowtie2-build" -f temp.fasta spadeContigs

# Align back with Bowtie2
$tool_dir"bowtie2-2.1.0/bowtie2" --very-sensitive-local -I 0 -X 1000 -p $thread -t -x spadeContigs \
-1 $read1 -2 $read2 -S $samfile.sam 

# Do pileup for each contig
$tool_dir"samtools-0.1.19/samtools" view -bSf 0x003 $samfile".sam" | $tool_dir"samtools-0.1.19/samtools" sort - $samfile"_sorted"
$tool_dir"samtools-0.1.19/samtools" index $samfile"_sorted.bam"
# in future samtool 1.1 version, -D is replaced with -t DB
$tool_dir"samtools-0.1.19/samtools" mpileup -f temp.fasta $samfile"_sorted.bam" > $samfile".pile"

# Tabulate contig coverage normalized by depth (could use python)
# this statement allows me to sum all the numbers in column 4 corresponding to each id
# in column 1. The output is contig id \t total depth of coverage (not normalized by contig length)
sort -k1 $samfile".pile" > $samfile"_sorted.pile"
awk 'a1==$1 {a2+=$4; next} {print a1,"\t", a2; a1=$1; a2=$4} END {print a1,"\t",a2}' $samfile"_sorted.pile" > $report

echo "Alignment Completed"

date
exit


