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
module load picard/2.27.5-Java-1.8.0_201

# check to make sure picard version is right
echo "Picard version" $(picard MarkDuplicates --version) 


##############################################################################
##																			                                    ##
##				               STEP 2: Run picard tools                           ## 
##																		                                      ##	
##############################################################################


# Run picard -----
for j in ./crams/AFR/*.final.cram #need to change Bowtie/ here if we use it for bwa-mem
do
newj=$(echo $(basename ${j}) | sed 's/.final.cram//g')
start=`date +%s` # capture start time 
echo "Calculating stats for $(basename ${j})"

picard CollectInsertSizeMetrics \
I=${j} \
O=results/bowtie2/assembly/Read-Prep/Bowtie/Aligned/stats/${newj}_insertSize_metrics.txt \
H=results/bowtie2/assembly/Read-Prep/Bowtie/Aligned/stats/${newj}_insertSize_hist.pdf \
R=/home/groups/h3abionet/RefGraph/data/genomes/human/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
M=0.5 

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run picard stats for $(basename ${j})"

done