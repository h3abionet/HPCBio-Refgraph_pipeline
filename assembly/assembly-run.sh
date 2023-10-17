#!/bin/bash

#SBATCH --mem 50G
#SBATCH --job-name assembly
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -p hpcbio
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out


### This Runs Nextflow assembly UIUC pipeline 
## Date File Created: July 28th, 2021


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021

# Load nextflow ------
module load nextflow/21.04.1-Java-1.8.0_152

# Run nextflow UIUC workflow -----
nextflow run HPCBio-Refgraph_pipeline/assembly/assemble_bowtie2.nf \
-c HPCBio-Refgraph_pipeline/assembly/main.conf \
-qs 1 -resume \
-with-report nextflow_reports/assembly_nf_report.html \
-with-timeline nextflow_reports/assembly_nf_timeline.html \
-with-trace nextflow_reports/assembly_nf_trace.txt

# -log custom.log  #add this for log not hidden

# -q  # Disable the printing of information to the terminal.

# Put stderr and stdout into a log ------
#2>&1 | tee ./assembly.log




