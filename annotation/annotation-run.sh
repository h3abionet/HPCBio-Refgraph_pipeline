#!/bin/bash

#SBATCH --mem 18G
#SBATCH --job-name annotation
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out

### This Runs Nextflow Annotation UIUC pipeline 
## Date File Created: Dec 5, 2021


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021

# Load nextflow ------
module load nextflow/21.04.1-Java-1.8.0_152

# Run nextflow UIUC workflow -----
nextflow run HPCBio-Refgraph_pipeline/annotation.nf \
-c HPCBio-Refgraph_pipeline/annotation-config.conf \
-qs 3 -resume \
-with-report nextflow_reports/annotation_nf_report.html \
-with-timeline nextflow_reports/annotation_nf_timeline.html \
-with-trace nextflow_reports/annotation_nf_trace.txt

# -log custom.log  #add this for log not hidden
# -q  # Disable the printing of information to the terminal.

# -with-report nf_exec_report_annotation.html \
# -with-timeline nf_timeline_annotation.html \
# -with-trace > nf_trace_annotation.txt \ # this is the same as slurm output, if you use this, slurm output will be empty
# -with-dag nf_flowchart_annotation.pdf

#if [ echo "wc -l ${keep}" == echo "grep -E ">" ${id}_kn_filtered.fasta | wc -l" ]
 #   then 
  #      echo "The filtering has not been done correctly. Please check your blastncontam script"
   # fi


   