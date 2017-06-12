#!/bin/bash

###################################
# Quantifying Assembly Statistics
###################################

scratch=$1
contigs=$2
report=$3
tool_dir=$4
threads=$5

echo $scratch
cd $scratch
date

source /local10G/brianyu/anaconda/bin/activate /local10G/brianyu/anaconda/

# Process scaffold stats with Quast
# commandline output goes to log file
python $tool_dir"/quast-2.3/quast.py" -T $threads $contigs -o quast_output 

mv quast_output/report.txt $report

echo "Quast Completed"
date
exit


