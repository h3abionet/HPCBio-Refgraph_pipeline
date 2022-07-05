#!/bin/bash

#SBATCH --mem 100G
#SBATCH --job-name compare-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/

echo "Working directory is set to" && pwd


assembly="megahit"



#module load seqkit/0.12.1

#echo "started seqkit"
#seqkit stats pipeline_testing/merged_AFR/${assembly}/*.fasta \
#--all --basename --tabular --threads 8 \
#--out-file pipeline_testing/SeqKit//${assembly}/All_seqkit_${assembly}.txt

#echo "finished seqkit"




# Load BLAST ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# # blast to compare two sequences --- 
# for j in pipeline_testing/merged_AFR/${assembly}/*seqkit_kn*.fasta
# do 
# start=`date +%s` # capture start time
# newj=$(echo $(basename ${j}) | sed 's/.fasta//g')
# echo "BLAST to reference genome is starting for ${j}"
# blastn -db ../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
# -query ${j} \
#  -outfmt "6 qseqid sseqid stitle pident \
#     length evalue qcovs bitscore mismatch \
#     gapopen qstart qend qlen sstart send slen" \
# -perc_identity 100 \
# -max_target_seqs 5 \
# -max_hsps 10 \
# -evalue 1e-5 \
# -num_threads $SLURM_NPROCS \
# -out pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv

# # Create headers for the blast file -----
# echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
# qlen,sstart,send,slen" | tr ',' '\t'| cat - pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv \
# > pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.tsv

# rm pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv

# end=`date +%s`
# runtime=$((end-start))
# runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
# echo "It took $runtime minutes to run blast ref genome for ${j}"

# done


# blast to compare two sequences ---
for i in pipeline_testing/merged_AFR/${assembly}/3.seqkit_kn_filtered_All_AFR_megahit.fasta
do 
start=`date +%s` # capture start time 
newi=$(echo $(basename ${i}) | sed 's/.fasta//g')
echo "Subject-Query BLAST is starting for ${i}"
blastn \
-query ${i} \
-subject /home/groups/h3abionet/RefGraph/results/2022-01-30-HGSVCv2-truthset/data/insertions.fasta \
 -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
-perc_identity 100 \
-max_target_seqs 5 \
-max_hsps 10 \
-evalue 1e-5 \
-num_threads $SLURM_NPROCS \
-out pipeline_testing/BLAST/${assembly}/${newi}_blast_compare.0.tsv

# Create headers for the blast file -----
echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
qlen,sstart,send,slen" | tr ',' '\t'| cat - pipeline_testing/BLAST/${assembly}/${newi}_blast_compare.0.tsv \
> pipeline_testing/BLAST/${assembly}/${newi}_blast_compare.tsv

rm pipeline_testing/BLAST/${assembly}/${newi}_blast_compare.0.tsv

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run blast subject-query for ${i}"

done


#for i in pipeline_testing/BLAST/${assembly}/GRCH38_p0_filter_All_AFR_megahit_blast_ref.tsv
#do 
#echo $(basename ${i}) >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.00.tsv
#awk -F"\t" '($4 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.11.tsv
#awk -F"\t" '($4 == 100 && $7 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.22.tsv
#done 

#paste pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.00.tsv \
#pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.11.tsv \
#pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.22.tsv \
#| sort -r -k1,1 > pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref0.tsv

# Create headers -----
#echo -e "seqid,pident100_ref,pident100_qcovs100_ref" \
#| tr ',' '\t'| cat - pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref0.tsv \
#> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_${assembly}.tsv

#echo "table created"