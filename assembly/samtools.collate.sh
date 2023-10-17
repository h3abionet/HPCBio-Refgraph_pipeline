#!/bin/bash

#SBATCH --mem 50G
#SBATCH --job-name  collate_convert
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -p hpcbio
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out


### This will sort the cram file based on read names
## Date File Created: April 24, 2023

# Set working directory -----
cd /home/n-z/valizad2/NeginV_Test_Summer2021/crams/AFR/test/test/test/

# Load the modules -----
module load SAMtools/1.12-IGB-gcc-8.2.0


# Create prefix directory for the temporary files ----
# If the prefix option not used, will run out of memory trying to write to tmp/

# Sort based on read names -----

tmpdir=$( mktemp -d collate.XXXXXXXXX )

    samtools collate \
    -@ 8 \
    --reference /home/groups/h3abionet/RefGraph/data/genomes/human/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -O ../HG02011.final.cram \
    $tmpdir/HG02011 | samtools fastq \
     -1 HG02011.cramconverted.R1.fastq.gz \
     -2 HG02011.cramconverted.R2.fastq.gz

# This one worked
# samtools collate -u \
# --threads 8 \
# --reference /home/groups/h3abionet/RefGraph/data/genomes/human/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# -O HG02011.final.cram \
# tmp.fq.gz \
# | samtools fastq \
# -1 HG02011.cramconverted.R1.collated.fastq.gz \
# -2 HG02011.cramconverted.R2.collated.fastq.gz 

