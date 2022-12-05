#!/bin/bash

#SBATCH --mem 18G
#SBATCH --job-name assembly
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out


### This Runs Nextflow assembly UIUC pipeline 
## Date File Created: July 28th, 2021


# Load modules -----
module load Bowtie2/2.4.2-IGB-gcc-8.2.0
module load SAMtools/1.12-IGB-gcc-8.2.0

# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38_bowtie2/


# # Create a bowtie2 index --------
# bowtie2-build --threads ${SLURM_NTASKS} \
# GRCh38_bowtie2/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# GRCH38.decoy.hla



# To align fastq sequences using bowtie2 -------
   
    # PE
    bowtie2 -p ${SLURM_NTASKS} -x GRCH38_decoy_hla \
    -1 /home/n-z/valizad2/NeginV_Test_Summer2021/results/assembly/Read-Prep/Bowtie/cram-to-fastq/NA19238.cramconverted.R1.fastq.gz \
    -2 /home/n-z/valizad2/NeginV_Test_Summer2021/results/assembly/Read-Prep/Bowtie/cram-to-fastq/NA19238.cramconverted.R2.fastq.gz \
    | samtools view -bS - > /home/n-z/valizad2/NeginV_Test_Summer2021/results/assembly/Read-Prep/Bowtie/NA19238.bowtie.pe.bam

