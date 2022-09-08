#!/bin/bash

#SBATCH --mem 100G
#SBATCH --job-name compare-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out

# -n 8 #core, cpu, thread, or processor
#-N 1 #node or computer

### HPCBio: This merges AFR files at each step of pipeline and give blast and seqkit stats; 
# Created by Negin Valizadegan March 5, 2022; valizad2@illinois.edu



##############################################################################
##																			##
##				                 STEP 1: SET UP                             ##
##																			##	
##############################################################################


# Set working directory -------
cd /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/
echo "Working directory is set to" && pwd

assembly="megahit"


##############################################################################
##																			##
##				            STEP 2: Merge AFR sequences                         ## 
##																		    ##	
##############################################################################

# Create a directory for merged files ------
mkdir -p pipeline_testing/merged_AFR/${assembly}/

echo "Sequence files from various RefGraph pipeline steps are merging"

# Final assembly -----
rm pipeline_testing/merged_AFR/${assembly}/1.Final_Assembly_All_AFR_${assembly}.fasta # you can add -f to get rid of error about file not existing.
cat assembly/Final-Assembly/${assembly}/AFR/*.${assembly}.final.fasta \
>> pipeline_testing/merged_AFR/${assembly}/1.Final_Assembly_All_AFR_${assembly}.fasta

# Filted to >500 length sequence by seqkit -----
rm pipeline_testing/merged_AFR/${assembly}/2.seqkit_500_All_AFR_${assembly}.fasta
cat filter/SeqKit/${assembly}/*_minLength_500.fasta \
>> pipeline_testing/merged_AFR/${assembly}/2.seqkit_500_All_AFR_${assembly}.fasta

# Filtered by Kraken -----
rm pipeline_testing/merged_AFR/${assembly}/3.seqkit_kn_filtered_All_AFR_${assembly}.fasta
cat filter/Kraken2/${assembly}/*_seqkit_kn_filtered.fasta \
>> pipeline_testing/merged_AFR/${assembly}/3.seqkit_kn_filtered_All_AFR_${assembly}.fasta

# Filtered by blast NT -----
rm pipeline_testing/merged_AFR/${assembly}/4.seqkit_kn_blast_filtered_All_AFR_${assembly}.fasta
cat filter/BLAST_Contam/${assembly}/*_seqkit_kn_blast_filtered.fasta \
>> pipeline_testing/merged_AFR/${assembly}/4.seqkit_kn_blast_filtered_All_AFR_${assembly}.fasta

# Clustered by CD-HIT -----
rm pipeline_testing/merged_AFR/${assembly}/5.seqkit_kn_blast_cdhit_All_AFR_${assembly}.fasta
cat filter/CD-HIT/${assembly}/*_seqkit_kn_blast_cdhit.fasta \
>> pipeline_testing/merged_AFR/${assembly}/5.seqkit_kn_blast_cdhit_All_AFR_${assembly}.fasta

# Filtered by blast GRCH38.p0 -----
rm pipeline_testing/merged_AFR/${assembly}/6.GRCH38_p0_filter_All_AFR_${assembly}.fasta
cat filter/Final-Filtered/${assembly}/*_GRCH38_p0_filter.final.fasta \
>> pipeline_testing/merged_AFR/${assembly}/6.GRCH38_p0_filter_All_AFR_${assembly}.fasta

# Filtered by blast GRCH38.decoys -----
rm pipeline_testing/merged_AFR/${assembly}/6.GRCH38_decoys_filter_All_AFR_${assembly}.fasta
cat filter/Final-Filtered/${assembly}/*_GRCH38_decoys_filter.final.fasta \
>> pipeline_testing/merged_AFR/${assembly}/6.GRCH38_decoys_filter_All_AFR_${assembly}.fasta

# Filtered by blast CHM13 -----
rm pipeline_testing/merged_AFR/${assembly}/6.CHM13_filter_All_AFR_${assembly}.fasta
cat filter/Final-Filtered/${assembly}/*_CHM13_filter.final.fasta \
>> pipeline_testing/merged_AFR/${assembly}/6.CHM13_filter_All_AFR_${assembly}.fasta

# CD-HIT in annotation by blast GRCH38.p0 -----
rm pipeline_testing/merged_AFR/${assembly}/7.clustered_*.fasta
for file in annotation/Cluster_CDHIT/${assembly}/clustered_*.fasta 
do
newfile=$(echo $(basename ${file}) | sed 's/.fasta//g')
cp -a ${file} pipeline_testing/merged_AFR/${assembly}/7.${newfile}_${assembly}.fasta
done

echo "All sequencing files are merged"


##############################################################################
##																			                                    ##
##				         STEP 3: Counts of Contigs: SeqKit                        ## 
##																		                                      ##	
##############################################################################


# Load seqkit
module purge
module load seqkit/0.12.1

# Create seqkit folder ------
mkdir -p pipeline_testing/SeqKit/${assembly}/

# Run seqkit -------
echo "Started SeqKit"

seqkit stats pipeline_testing/merged_AFR/${assembly}/*.fasta \
--all --basename --tabular --threads 8 \
--out-file pipeline_testing/SeqKit/${assembly}/All_seqkit_${assembly}.tsv

echo "SeqKit was completed"


##############################################################################
##																			##
##				         STEP 3: Counts of Contigs: QUAST                ## 
##																		    ##	
##############################################################################

quast ()
{
# Load quast
module purge
module load quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1

# Run QUAST ------
mkdir -p pipeline_testing/QUAST/${assembly}/

for k in pipeline_testing/merged_AFR/${assembly}/*.fasta
do 
newk=$(echo $(basename ${k}) | sed 's/.fasta//g')
echo "QUAST is starting for ${k}"
quast.py ${k} -o pipeline_testing/QUAST/${assembly}/${newk}_quast
echo "QUAST was completed for ${k}"
done
}

##############################################################################
##																			##
##				       STEP 4: BLAST AGAINST TEST DATA                      ## 
##																		    ##	
##############################################################################

blast_test() 
{

# Load BLAST ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

mkdir -p pipeline_testing/BLAST/${assembly}/

# blast to compare two sequences ---
for i in pipeline_testing/merged_AFR/${assembly}/*.fasta
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

}


##############################################################################
##																			##
##				     STEP 5: BLAST AGAINST HUMAN REFERENCE GENOME           ## 
##																		    ##	
##############################################################################

blast_ref () 
{

# Load BLAST ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# blast to compare two sequences --- 
for j in pipeline_testing/merged_AFR/${assembly}/*.fasta
do 
start=`date +%s` # capture start time
newj=$(echo $(basename ${j}) | sed 's/.fasta//g')
echo "BLAST to reference genome is starting for ${j}"
blastn -db ../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-query ${j} \
 -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
-perc_identity 100 \
-max_target_seqs 5 \
-max_hsps 10 \
-evalue 1e-5 \
-num_threads $SLURM_NPROCS \
-out pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv

# Create headers for the blast file -----
echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
qlen,sstart,send,slen" | tr ',' '\t'| cat - pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv \
> pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.tsv

rm pipeline_testing/BLAST/${assembly}/${newj}_blast_ref.0.tsv

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run blast ref genome for ${j}"

done

}

##############################################################################
##																			##
##				            STEP 6: Count number of matches                 ## 
##																		    ##	
##############################################################################

count ()

{
echo "count process started"

# Remove if the file exists -----
rm pipeline_testing/Count/${assembly}/shared_sequences*.tsv
mkdir -p pipeline_testing/Count/${assembly}/

# Get number of reads from seqkit results -----
head -n 1 pipeline_testing/SeqKit/${assembly}/All_seqkit_${assembly}.tsv \
&& tail -n +2 pipeline_testing/SeqKit/${assembly}/All_seqkit_${assembly}.tsv \
| sort -k1,1 > pipeline_testing/Count/${assembly}/shared_sequences_qstats0.tsv

# Create headers -----
echo -e "seqid,format,type,num_seqs,sum_len,min_len,avg_len,max_len,Q1,Q2,Q3,sum_gap,N50,Q20(%),Q30(%)" \
| tr ',' '\t' | cat - pipeline_testing/Count/${assembly}/shared_sequences_qstats0.tsv \
> pipeline_testing/Count/${assembly}/shared_sequences_qstats_${assembly}.tsv

# Count number of pident of 100 and qcov 100 for each file and put into a table ------
for i in pipeline_testing/BLAST/${assembly}/*blast_compare.tsv
do 
echo $(basename ${i}) >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.0.tsv
awk '($4 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.1.tsv
awk '($4 == 100 && $7 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.2.tsv
done 

# Put the two above into columns of a table -----
paste pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.0.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.1.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.2.tsv \
| sort -k1,1 > pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare0.tsv

# Create headers -----
echo -e "seqid,pident100_insertions,pident100_qcovs100_insertions" \
| tr ',' '\t'| cat - pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare0.tsv \
> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_${assembly}.tsv


# Count number of pident of 100 and qcov 100 for each file and put into a table ------
for i in pipeline_testing/BLAST/${assembly}/*blast_ref.tsv
do 
echo $(basename ${i}) >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.00.tsv
awk -F"\t" '($4 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.11.tsv
awk -F"\t" '($4 == 100 && $7 == 100)' ${i} | cut -f1 | uniq | wc -l >> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.22.tsv
done 

paste pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.00.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.11.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.22.tsv \
| sort -k1,1 > pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref0.tsv

# Create headers -----
echo -e "seqid,pident100_ref,pident100_qcovs100_ref" \
| tr ',' '\t'| cat - pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref0.tsv \
> pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_${assembly}.tsv

cat pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_${assembly}.tsv \
| cut -f2- > pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_merge.tsv

cat pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_${assembly}.tsv \
| cut -f2- > pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_merge.tsv

#echo "Assembly: Final assembly file,Filter 1: After limiting to > 500 length, \
#Filter 2: Filtered based on kraken,Filter 3: Filtered based on kraken and blast nt, \
#Filter 4: Filtered based on kraken and blast nt and clustered, \
#Filter 5: Final filtered file (GRCH38_p0),Filter 5: Final filtered file (GRCH38_decoys_hla), \
#Filter 5: Final filtered file (CHM13), Annotation: File clustered after filtering (GRCH38_p0), \
#Annotation: File clustered after filtering (GRCH38_decoys_hla),Annotation: File clustered after filtering (CHM13)" \
#| tr ',' '\t' > pipeline_testing/Count/${assembly}/descriptions.tsv

#ls -v pipeline_testing/Count/${assembly}/shared_sequences_qstats_${assembly}.tsv \
#pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_megahit.tsv \
#pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_megahit.tsv \
#| xargs paste > data.txt

paste -d '\t' pipeline_testing/Count/${assembly}/shared_sequences_qstats_${assembly}.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_merge.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_merge.tsv \
> pipeline_testing/Count/${assembly}/Final_stats_${assembly}.tsv

# Remove extra files and view final file -----
rm pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.0.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.1.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.2.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.00.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.11.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.22.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare0.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref0.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_qstats0.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_ref_merge.tsv \
pipeline_testing/Count/${assembly}/shared_sequences_pident100_qcovs100.blast_compare_merge.tsv 

echo "Created a table with pident=100 and qcovs=100 and number of sequences from seqkit"
}


##############################################################################
##																			                                    ##
##						                      MAIN                                    ##
##																		                                      ##
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
  
quast
blast_test
blast_ref
count

}

# Run main function
main

##############################################################################
##																		                                    	##
##				                     STEP 7: MultiQC                              ##
##																		                                      ##	
##############################################################################


# Load MultiQC ------  # multiqc has issues with processing () so do not do that
  module purge
  module load MultiQC/1.11-IGB-gcc-8.2.0-Python-3.7.2
  mkdir -p pipeline_testing/MultiQC/

  # Run multiqc

  echo "started multiqc"

  multiqc -f pipeline_testing/ -o pipeline_testing/MultiQC/ # -f will force it to overwrite previous runs

  echo "completed multiqc"
  
