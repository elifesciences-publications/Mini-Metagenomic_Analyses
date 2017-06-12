#!/bin/bash

###################################
# Using Bowtie2 to align all 
# trimmed reads back to contigs
###################################

scratch=$1
read1=$2
read2=$3
contigs=$4
report=$5
tool_dir=$6
thread=$7

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/
samfile="alignResults"

# Build Bowtie2 index
$tool_dir"bowtie2-2.1.0/bowtie2-build" -f $contigs spadeContigs

# Align back with Bowtie2
$tool_dir"bowtie2-2.1.0/bowtie2" --very-sensitive-local -I 0 -X 1000 -p $threads -t -x spadeContigs \
-1 $read1 -2 $read2 -S $samfile.sam 
echo -e "New_Sample\t"$sample_id > $report

# Get all coordinates
echo -e "Total_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x001 $samfile".sam" | cut -d: -f 4-7 \
| cut -f 1 | sort | uniq | wc -l ) >> $report

# Uniquely mapped coordinates
echo -e "Uniquely_Mapped_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x003 $samfile".sam" | cut -d: -f 4-7 | cut -f 1-9 \
| sort | uniq | wc -l ) >> $report

# All mapped coordinates (even only one of the paired reads
echo -e "All_Mapped_Reads\t"$( $tool_dir"samtools-0.1.19/samtools" view -Sf 0x001 -F 0x00c $samfile".sam" | cut -d: -f 4-7 | cut -f 1-9 \
| sort | uniq | wc -l ) >> $report

echo "Alignment Completed"

date
exit


