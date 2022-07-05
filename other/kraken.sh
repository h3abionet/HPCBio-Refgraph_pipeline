#!/bin/bash

#SBATCH --mem 70G
#SBATCH --job-name kraken
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH --output slurm-%j.out
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



# HPCBio UIUC Kraken; Created by Negin Valizadegan Oct 6, 2021; valizad2@illinois.edu

# ----- Load necessary modules ------ 
module load Kraken2/2.0.8-beta-IGB-gcc-4.9.4

# Set working directory ------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/

echo "Input directory is set to" | tr '\n' ' ' && pwd

# Create a directory for outputs ----
mkdir -p annotation/kraken/masurca
mkdir -p annotation/kraken/megahit

echo "Start of kraken process masurca"

# Masurca -----

# Run Kraken ------

# kraken2 --use-names --threads 6 --quick   \
# --db /home/groups/h3abionet/RefGraph/data/kraken2/human \
# annotation/seqkit/masurca/HG03563.masurca.filtered.fasta

kraken2 --use-names --threads 6 --quick   \
--report annotation/kraken/masurca/HG02465_kraken2_report.txt \
--classified-out annotation/kraken/masurca/HG02465_masurca_kraken2_classified.fasta \
--unclassified-out annotation/kraken/masurca/HG02465_masurca_kraken2_unclassified.fasta \
--db /home/groups/h3abionet/RefGraph/data/kraken2/pluspf_20200919 \
annotation/seqkit/masurca/HG02465.masurca.filtered.fasta > annotation/kraken/masurca/HG02465_masurca_kraken2_output.txt

echo "End of kraken process masurca"