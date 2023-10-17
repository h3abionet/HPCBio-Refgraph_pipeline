#!/bin/bash

#SBATCH --mem 30G
#SBATCH --job-name compare-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -p hpcbio
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out

# -n 8 #core, cpu, thread, or processor
#-N 1 #node or computer

### HPCBio: This will give seqkit stats for output of each assembly step; 
# Created by Negin Valizadegan May 24, 2023; valizad2@illinois.edu


##############################################################################
##																			##
##				                 STEP 1: SET UP                             ##
##																			##	
##############################################################################


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/
echo "Working directory is set to" && pwd

assembly="masurca"
alignment="bowtie2"


##############################################################################
##																			##
##				         STEP 2: SeqKit Read-Prep/Bowtie                    ## 
##																		    ##	
##############################################################################


# Load seqkit -----
module purge
module load seqkit/0.12.1

# Variables -----
#folder="Unmapped"
folder="cram-to-fastq"

# Create seqkit folder -----
mkdir -p ${alignment}/assembly/Read_testing/Read-Prep/

# Run seqkit -----
echo "Started SeqKit"

seqkit stats ${alignment}/assembly/Read-Prep/Bowtie/${folder}/*.fastq* \
--all --basename --tabular --threads 8 \
--out-file ${alignment}/assembly/Read_testing/Read-Prep/seqkit_Read-Prep_${folder}.tsv

echo "SeqKit was completed"