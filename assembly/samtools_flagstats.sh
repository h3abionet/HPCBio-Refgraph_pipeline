#!/bin/bash

#SBATCH --mem 30G
#SBATCH --job-name assembly_qc
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -p hpcbio
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%A.out


### This runs QC for assembled files 
## Date File Created: Jan 31, 2023


##############################################################################
##																			                                    ##
##				                    STEP 1: SET UP                                ##
##																			                                    ##	
##############################################################################

# Set working directory and parameters -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/
alignment="bowtie2"

mkdir -p results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/stats

# Load modules ------
module load SAMtools/1.12-IGB-gcc-8.2.0


##############################################################################
##																			                                    ##
##				               STEP 2: Run samtools                               ## 
##																		                                      ##	
##############################################################################

for i in results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/*.pe.bam #need to change Bowtie/ here if we use it for bwa-mem
do
newi=$(echo $(basename ${i}) | sed 's/.pe.bam//g')
start=`date +%s` # capture start time 
echo "Calculating stats for $(basename ${i})"

samtools flagstat ${i} \
--output-fmt tsv \
--threads 8 \
> results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/stats/${newi}0.tsv #need to change Bowtie/ here if we use it for bwa-mem

# Create headers for the blast file -----
echo -e "QC-passed,QC-failed,Category" | tr ',' '\t'| cat - results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/stats/${newi}0.tsv \
> results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/stats/${newi}.tsv

rm results/${alignment}/assembly/Read-Prep/Bowtie/Aligned/stats/${newi}0.tsv

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run samtools stats for $(basename ${i})"

done 