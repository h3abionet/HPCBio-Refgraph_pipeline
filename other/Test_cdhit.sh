#!/bin/bash

#SBATCH --mem 18G
#SBATCH --job-name cd-hit-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out

### This tests various parameters for cdhit at the annotation step UIUC pipeline 
## Date File Created: April 23, 2022


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/annotation

# Load nextflow ------
module load CD-HIT/4.8.1-IGB-gcc-8.2.0


# Use cd-hit to cluster and remove redundancy ------

# GRCH38_decoys -----

# megahit -----

# 0.9
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_decoys_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS

# 0.92
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_decoys_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS

# 0.94
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_decoys_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS

# 0.96
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_decoys_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS


# masurca -----

# 0.9
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_decoys_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS

# 0.92
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_decoys_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS

# 0.94
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_decoys_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS

# 0.96
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_decoys_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS


# GRCH38_p0 -----

# megahit -----

# 0.9
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_p0_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS  

# 0.92
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_p0_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS

# 0.94
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_p0_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS  

# 0.96
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/megahit/clustered_GRCH38_p0_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS  


# masurca -----

# 0.9
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_p0_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS 

# 0.92
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_p0_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS 

# 0.94
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_p0_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS 

# 0.96
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/masurca/clustered_GRCH38_p0_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS 


# CHM13 -----

# megahit

# 0.9
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/megahit/clustered_CHM13_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS 

# 0.92
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/megahit/clustered_CHM13_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS 

# 0.94
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/megahit/clustered_CHM13_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS 

# 0.96
cd-hit-est \
-i Merged_Reads/megahit/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/megahit/clustered_CHM13_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS 


# masurca -----

# 0.9
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/masurca/clustered_CHM13_n5_i0.9.fasta \
-c 0.9 \
-n 5 \
-T $SLURM_NPROCS 

# 0.92
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/masurca/clustered_CHM13_n5_i0.92.fasta \
-c 0.92 \
-n 5 \
-T $SLURM_NPROCS 

# 0.94
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/masurca/clustered_CHM13_n5_i0.94.fasta \
-c 0.94 \
-n 5 \
-T $SLURM_NPROCS 

# 0.96
cd-hit-est \
-i Merged_Reads/masurca/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/masurca/clustered_CHM13_n5_i0.96.fasta \
-c 0.96 \
-n 5 \
-T $SLURM_NPROCS 