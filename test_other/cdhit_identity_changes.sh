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
## Edited: Aug 19, 2022



# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/annotation

# Load nextflow ------
module load CD-HIT/4.8.1-IGB-gcc-8.2.0

# Variables ----
identity_array=(90 92 94 96 99)
assembly="masurca"


## Use cd-hit to cluster and remove redundancy ------

##############################################################################
##																			##
##				                 GRCH38_decoys                              ##
##																			##	
##############################################################################

GRCH38_decoys ()

{
for i in ${identity_array[@]}
do
if (( ${i} >= 90 && ${i} < 95 )); then
wordsize='8'
echo "${i} is between 90 - 95 so word size is ${wordsize}"
elif (( ${i} >= 95 && ${i} < 100 )); then
wordsize='10'
echo "${i} is between 95 - 100 so word size is ${wordsize}"
else
    echo "[[${i} is out of defined range in this script. Please edit the code to define. See below]]"
    echo "
#Choose of word size: 
#For DNAs:
 # * Word size 10-11 is for thresholds 0.95 ~ 1.0
 # * Word size 8,9 is for thresholds 0.90 ~ 0.95
 # * Word size 7 is for thresholds 0.88 ~ 0.9
 # * Word size 6 is for thresholds 0.85 ~ 0.88
 # * Word size 5 is for thresholds 0.80 ~ 0.85
 # * Word size 4 is for thresholds 0.75 ~ 0.8
 "
fi
newi=$(echo "0.$i")
cd-hit-est \
-i Merged_Reads/${assembly}/merged_sequences_GRCH38_decoys.fasta \
-o Cluster_CDHIT/${assembly}/clustered_GRCH38_decoys_i${newi}_n${wordsize}.fasta \
-c ${newi} \
-n ${wordsize} \
-T $SLURM_NPROCS
done
}

##############################################################################
##																			##
##				                 GRCH38_p0                                  ##
##																			##	
##############################################################################

GRCH38_p0 ()

{
for i in ${identity_array[@]}
do
if (( ${i} >= 90 && ${i} < 95 )); then
wordsize='8'
echo "${i} is between 90 - 95 so word size is ${wordsize}"
elif (( ${i} >= 95 && ${i} < 100 )); then
wordsize='10'
echo "${i} is between 95 - 100 so word size is ${wordsize}"
else
    echo "[[${i} is out of defined range in this script. Please edit the code to define. See below]]"
    echo "
#Choose of word size: 
#For DNAs:
 # * Word size 10-11 is for thresholds 0.95 ~ 1.0
 # * Word size 8,9 is for thresholds 0.90 ~ 0.95
 # * Word size 7 is for thresholds 0.88 ~ 0.9
 # * Word size 6 is for thresholds 0.85 ~ 0.88
 # * Word size 5 is for thresholds 0.80 ~ 0.85
 # * Word size 4 is for thresholds 0.75 ~ 0.8
 "
fi
newi=$(echo "0.$i")
cd-hit-est \
-i Merged_Reads/${assembly}/merged_sequences_GRCH38_p0.fasta \
-o Cluster_CDHIT/${assembly}/clustered_GRCH38_p0_i${newi}_n${wordsize}.fasta \
-c ${newi} \
-n ${wordsize} \
-T $SLURM_NPROCS
done
}

##############################################################################
##																			##
##				                   CHM13                                    ##
##																			##	
##############################################################################

CHM13 ()

{
for i in ${identity_array[@]}
do
if (( ${i} >= 90 && ${i} < 95 )); then
wordsize='8'
echo "${i} is between 90 - 95 so word size is ${wordsize}"
elif (( ${i} >= 95 && ${i} < 100 )); then
wordsize='10'
echo "${i} is between 95 - 100 so word size is ${wordsize}"
else
    echo "[[${i} is out of defined range in this script. Please edit the code to define. See below]]"
    echo "
#Choose of word size: 
#For DNAs:
 # * Word size 10-11 is for thresholds 0.95 ~ 1.0
 # * Word size 8,9 is for thresholds 0.90 ~ 0.95
 # * Word size 7 is for thresholds 0.88 ~ 0.9
 # * Word size 6 is for thresholds 0.85 ~ 0.88
 # * Word size 5 is for thresholds 0.80 ~ 0.85
 # * Word size 4 is for thresholds 0.75 ~ 0.8
 "
fi
newi=$(echo "0.$i")
cd-hit-est \
-i Merged_Reads/${assembly}/merged_sequences_CHM13.fasta \
-o Cluster_CDHIT/${assembly}/clustered_CHM13_i${newi}_n${wordsize}.fasta \
-c ${newi} \
-n ${wordsize} \
-T $SLURM_NPROCS
done
}

##############################################################################
##																			##
##						           MAIN                                     ##
##																		    ##
##############################################################################

# Main function runs each step/function of the pipeline separately so that
# user can choose to run steps one at a time.

main ()
{
	# Determine whether running full pipeline or single step
	runtype="PARTIAL"
    #runtype="FULL"
    echo ""
	echo "*** RUNNING ${runtype} BLAST COMPARE PIPELINE ***"
  
GRCH38_decoys
GRCH38_p0
CHM13

}

# Run main function
main