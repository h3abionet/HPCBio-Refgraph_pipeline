#!/bin/bash

#SBATCH --mem 18G
#SBATCH --job-name assembly
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH --output slurm-%j.out
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -A h3abionet


### This Runs Nextflow assembly UIUC pipeline 
## Date File Created: July 28th, 2021


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021

# Load nextflow ------
module load nextflow/21.04.1-Java-1.8.0_152

# Run nextflow UIUC workflow -----
nextflow run HPCBio-Refgraph_pipeline/assemble.nf -c HPCBio-Refgraph_pipeline/test-config.conf -qs 1 -resume 

# -log custom.log  #add this for log not hidden

# -q  # Disable the printing of information to the terminal.

# Put stderr and stdout into a log ------
2>&1 | tee ./assembly.log




