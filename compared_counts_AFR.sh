#!/bin/bash

#SBATCH --mem 10G
#SBATCH --job-name compare-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



### HPCBio RefGraph Pipeline Testing Count: This compares test data set's long read insertions with our files; Created by Negin Valizadegan Feb 19, 2022; valizad2@illinois.edu


##############################################################################
##																			##
##				                 STEP 1: SET UP                             ##
##																			##	
##############################################################################

# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/pipeline_testing

echo "output directory is set"


##############################################################################
##																			##
##				             STEP 2: Count number of matches                ## 
##																		    ##	
##############################################################################


# Remove if the file exists -----
rm shared_sequences_pident100_qcovs100*.txt

# Count number of pident of 100 and qcov 100 for each file and put into a table ------
for i in *.txt
do 
echo ${i} >> shared_sequences_pident100_qcovs100.0.txt
awk '($4 == 100)' ${i} | cut -f1 | uniq | wc -l >> shared_sequences_pident100_qcovs100.1.txt
awk '($4 == 100 && $7 == 100)' ${i} | cut -f1 | uniq | wc -l >> shared_sequences_pident100_qcovs100.2.txt
done 

echo "numbers with pident == 100 and qcovs == 100 added"

# Put the two above into columns of a table -----
paste shared_sequences_pident100_qcovs100.0.txt shared_sequences_pident100_qcovs100.1.txt shared_sequences_pident100_qcovs100.2.txt | sort -r -k2,2 > shared_sequences_pident100_qcovs100.3.txt

# Create headers -----
echo -e "qseqid,pident100,pident100_qcovs100" | tr ',' '\t'| cat - shared_sequences_pident100_qcovs100.3.txt > shared_sequences_pident100_qcovs100.txt

# Remove extra files and view final file -----
rm shared_sequences_pident100_qcovs100.0.txt shared_sequences_pident100_qcovs100.1.txt shared_sequences_pident100_qcovs100.2.txt shared_sequences_pident100_qcovs100.3.txt
cat shared_sequences_pident100_qcovs100.txt





