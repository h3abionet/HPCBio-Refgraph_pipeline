#!/bin/bash

#SBATCH --mem 70G
#SBATCH --job-name annotation
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH --output slurm-%j.out
#SBATCH -n 6
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



# HPCBio UIUC Annotation pipeline (Filtering); Created by Negin Valizadegan Sep 2, 2021; valizad2@illinois.edu

##############################################################################
##																		                                    	##
##			           GENERAL WRAPPER RELATED SCRIPTS	                        ##
##																			                                    ##	
##############################################################################

# Set fancy fonts for the help message ------
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# Help ------
function HELP {
  echo ""
  echo "${BOLD}Help Documentation for the HPCBio UIUC Annotation (Filtering) Pipeline${NORM}"
  echo ""
  echo "The Following Options Must Be Specified:"
  echo "${REV}-d${NORM}   The full path to the main results directory${NORM} (Required)"
  echo "${REV}-s${NORM}   The name of the input sequence (Required)"
  echo "${REV}-h${NORM}   Displays this help message without complaints (Optional)"
  echo ""
  echo "[ ${NORM}${BOLD}Example:${NORM} sbatch filter.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/ -s HG02461.masurca.final.fasta ]"
  echo ""
  exit 1
}


# Check the number of arguments. If none are passed, print message and exit ------
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  echo ""
  echo "You Did Not Pass Any Arguments. Please Specify the Arguments Below:" 
  echo ""
  HELP
fi


# Parse the inputs
while getopts :d:s:h FLAG; do
  case $FLAG in
    d)  #set option "d"
      OPT_d=$OPTARG
      ;;
    s)  #set option "s"
      OPT_s=$OPTARG
      ;;
     h)  #set option "h"
      OPT_h=$OPTARG
      HELP
      ;;
    \?) #unrecognized option - show help
      echo "Option -${BOLD}$OPTARG${NORM} not allowed."
      exit 1
      ;;
  esac
done


# Exit if necessary options are not passed ------
if [[ -z "$OPT_d" ]]; then
        echo "No project directory specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_s" ]]; then
        echo "No input sequence name is specified, aborting script"
        exit 1
fi



startall=`date +%s`  # record start time 

##############################################################################
##																			                                    ##
##			           STEP 1: LOAD THE NECESSARY MODULES                       ##
##																			                                    ##	
##############################################################################

setup ()
{

# Create a variable to be used later -----
news=$(echo ${OPT_s} | sed 's/.masurca.final.fasta//g')
echo "Running the pipeline for $news"

}

##############################################################################
##																			                                    ##
##				            STEP 2: FILTER THE ASSEMBLY FILES                     ##
##																			                                    ##	
##############################################################################

filter ()
{

# ----- Set working directory -----
cd ${OPT_d}

echo "Input directory is set to" | tr '\n' ' ' && pwd

# Create a directory for outputs ----
mkdir -p annotation/seqkit/masurca
mkdir -p annotation/seqkit/megahit

# Load seqkit ------
module purge
module load seqkit/0.12.1

# ----- Filter sequences < 500bp length -----
start=`date +%s` # capture start time 

echo "Start of filtering process masurca"

# Masurca -----
echo "[Removing reads with length < 500]"

seqkit seq --min-len 500 --remove-gaps assembly/Final-Assembly/masurca/${OPT_s} > annotation/seqkit/masurca/$news.masurca.filtered.fasta

echo "[Creating stats on fasta files]"

seqkit stats assembly/Final-Assembly/masurca/${OPT_s} > annotation/seqkit/masurca/$news.masurca.filter-stats.txt
seqkit stats annotation/seqkit/masurca/$news.masurca.filtered.fasta | sed -e '1d' >> annotation/seqkit/masurca/$news.masurca.filter-stats.txt
echo "End of Filtering Process"

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to filter [$news.masurca.final.fasta]"

}

##############################################################################
##																			                                    ##
##	               STEP 3: CREATE BLAST DATABASE (GRCh38)                   ##
##																			                                    ##	
##############################################################################

blastdbGRCh38 ()
{

# ----- Set working directory -----
cd ${OPT_d}

echo "Input directory is set to" | tr '\n' ' ' && pwd

# Load blast ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# ----- Creat Blast database from GRCh38.decoy.hla [?] -----
start=`date +%s`
makeblastdb -in ../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

# Building a new DB, current time: 09/13/2021 18:24:45
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
# New DB title:  GRCh38.decoy.hla
# Sequence type: Nucleotide
# Deleted existing Nucleotide BLAST database named /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 3366 sequences in 29.3228 seconds.

end=`date +%s` # record finish time 
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to create blast database from [GRCh38.decoy.hla]"


# ----- Creat Blast database from GRCh38.no.decoy.hla [GRCh38.p0; GCA_000001405.15] -----
start=`date +%s`
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
makeblastdb -in ../GRCh38.p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -parse_seqids -title "GRCh38.p0.no.decoy.hla" -dbtype nucl

# Building a new DB, current time: 09/13/2021 18:25:15
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38.p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# New DB title:  GRCh38.p0.no.decoy.hla
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 195 sequences in 28.1475 seconds.

end=`date +%s` # record finish time 
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to create blast database from [GRCh38.p0.no.decoy.hla]"

}


##############################################################################
##																			                                    ##
##			        STEP 4: RUN BLAST ON FILTERED READS (GRCh38)                ##
##																			                                    ##	
##############################################################################

blastnGRCh38 ()
{

# ----- Run blastn on filtered sequences from previous step ------

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to" | tr '\n' ' ' && pwd

# Make a directory for blast ------
mkdir -p annotation/blast/masurca
mkdir -p annotation/blast/megahit

echo "[Start of BLASTN:GRCh38 Process]"

# Load blast ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# Masurca -----

# Set directory to GRCh38 ------
cd ../GRCh38/


# BLAST against GRCh38.decoy.hla ------
echo "[Running BLAST+ on $news.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run blast ------
   blastn -db GRCh38_full_analysis_set_plus_decoy_hla.fa \
   -query ../results/annotation/seqkit/masurca/$news.masurca.filtered.fasta \
   -out ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.decoy.hla.asn \
   -outfmt 11 \
   -max_target_seqs 5 \
   -max_hsps 10 \
   -evalue 1e-5 \
   -perc_identity 90

echo "[Convert BLAST Archive (.asn) to Tabular Format]"

blast_formatter -archive ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.decoy.hla.asn \
-outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.decoy.hla.txt

end=`date +%s` # record finish time 
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [$news.masurca.filtered.fasta] against [GRCh38.decoy.hla]"

echo "Blasted $news.masurca.filtered.fasta against GRCh38.decoy.hla"


# BLAST against GRCh38.p0.no.decoy.hla ------
start=`date +%s`

# Download the ref sequence ------
# mkdir GRCh38.p0
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Set directory to GRCh38.p0 ------
cd ../GRCh38.p0/

# Run blast ------
    blastn -db GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -query ../results/annotation/seqkit/masurca/$news.masurca.filtered.fasta \
    -out ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.p0.no.decoy.hla.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.p0.no.decoy.hla.asn \
-outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> ../results/annotation/blast/masurca/blast_$news.masurca.GRCh38.p0.no.decoy.hla.txt

end=`date +%s`
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [$news.filtered.fasta] against [GRCh38.p0.no.decoy.hla]"

echo "Blasted $news.masurca.filtered.fasta against GRCh38.p0.no.decoy.hla"

echo "End of BLASTN GRCh38 process"

}


##############################################################################
##																			                                    ##
##	                 STEP 5: CREATE BLAST DATABASE (CHM13)                  ##
##																		                                    	##	
##############################################################################

blastdbCHM13()
{

# Load blast ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# ----- Creat Blast database from CHM13 [GCA_009914755.3_CHM13_T2T_v1.1] + chrY of GRCh38.p13 [GCA_000001405.28_GRCh38.p13] -----
start=`date +%s`
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_CHM13_T2T_v1.1/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chrY.fna.gz

# Combine CHM13 (female) + chrY of GRCh38.p13 ------
# cat GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna chrY.fna  > CHM13.v1.1_GRCh38.p13.chrY.fna

# blast ------
makeblastdb -in ../CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna -parse_seqids -title "CHM13.v1.1_GRCh38.p13.chrY" -dbtype nucl

# Building a new DB, current time: 09/15/2021 17:21:25
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna
# New DB title:  CHM13.v1.1_GRCh38.p13.chrY
# Sequence type: Nucleotide
# Deleted existing Nucleotide BLAST database named /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 25 sequences in 25.1712 seconds.

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to create blast database for [CHM13.v1.1_GRCh38.p13.chrY]"

}

##############################################################################
##																			                                    ##
##			        STEP 6: RUN BLAST ON FILTERED READS (CHM13)                 ##
##																			                                    ##	
##############################################################################

blastnCHM13 ()
{

# ----- Run blastn on filtered sequences from previous step ------

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to:"  |  tr '\n' ' ' && pwd 

# Make a directory for blast ------
mkdir -p annotation/blast/masurca
mkdir -p annotation/blast/megahit

echo "[Start of BLASTN:CHM13 Process]"

# Load blast ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# Masurca -----

# Set directory to CHM13.v1.1_GRCh38.p13.chrY ------
cd ../CHM13.v1.1_GRCh38.p13.chrY/


# BLAST against GRCh38.decoy.hla ------

echo "[Running BLAST+ on $news.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run blast ------
    blastn -db CHM13.v1.1_GRCh38.p13.chrY.fna \
    -query ../results/annotation/seqkit/masurca/$news.masurca.filtered.fasta \
    -out ../results/annotation/blast/masurca/blast_$news.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../results/annotation/blast/masurca/blast_$news.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
-outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> ../results/annotation/blast/masurca/blast_$news.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.txt

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [$news.masurca.filtered.fasta] against [CHM13.v1.1_GRCh38.p13.chrY]"

echo "Blasted $news.masurca.filtered.fasta against CHM13.v1.1_GRCh38.p13.chrY.fna"

echo "End of BLASTN CHM13 process"

}


##############################################################################
##																			                                    ##
##			        STEP 7: RUN REPEAT MASKER ON FILTERED READS                 ##
##																			                                    ##	
##############################################################################

repeatmasker ()
{

# ----- Run RepeatMasker on filtered sequences from previous step ------

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to" | tr '\n' ' ' && pwd

# Make a directory for blast ------
mkdir -p annotation/RepeatMasker/masurca
mkdir -p annotation/RepeatMasker/megahit

echo "[Start of RepeatMasker Process]"


# Load what we need -----
module purge
module load RepeatMasker/4.1.2-p1-IGB-gcc-8.2.0-Perl-5.28.1


# Masurca -----
echo "[Running RepeatMasker on $news.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run Repeat Masker ------
RepeatMasker -species human \
annotation/seqkit/masurca/$news.masurca.filtered.fasta \
-dir annotation/RepeatMasker/masurca

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes for RepeatMasker to run [$news.masurca.filtered.fasta]"


echo "End of RepeatMasker process"

}


##############################################################################
##																			                                    ##
##			       STEP 8: RUN QUAST ON FILTERED & MASKED READS                ##
##																			                                    ##	
##############################################################################

quast ()
{

# ----- Run RepeatMasker on filtered sequences from previous step ------

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to" | tr '\n' ' ' && pwd

# Make a directory for blast ------
mkdir -p annotation/QUAST/$news

echo "[Start of RepeatMasker Process]"

# Load quast ------
module purge
module load quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1

# Masurca -----
echo "[Running QUAST on $news.masurca.filtered.fasta"

start=`date +%s`  # record start time 

# Run QUAST ------
quast.py annotation/seqkit/masurca/$news.masurca.filtered.fasta \
--output-dir annotation/QUAST/$news \
--threads 2 \

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes for QUAST to run [$news.masurca.filtered.fasta]"


echo "End of QUAST process"

}

 
##############################################################################
##																			##
##			      STEP 9: RUN BLAST ON FILTERED READS              ##
##																			##	
##############################################################################

blastncontam ()
{

# ----- Run blastn on filtered sequences from previous step ------

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to:"  |  tr '\n' ' ' && pwd 

# Make a directory for blast ------
mkdir -p annotation/blast-contam/masurca
mkdir -p annotation/blast-contam/megahit

echo "[Start of BLASTN for contamination detetcion Process]"

# Load blast database ------
module purge
module load BLAST+/2.10.1-IGB-gcc-8.2.0
module load ncbi-blastdb/20201212

# Masurca -----

# BLAST against built in blast database ------
echo "[Running BLASTN on $news.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run blast ------
    blastn -db nt \
    -query annotation/seqkit/masurca/$news.masurca.filtered.fasta \
    -out annotation/blast-contam/masurca/blast_$news.masurca.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

# Convert blast archive to tabular -----
blast_formatter -archive annotation/blast-contam/masurca/blast_$news.masurca.asn \
-outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> annotation/blast-contam/masurca/blast_$news.masurca.txt

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [$news.masurca.filtered.fasta] against [nt]"

echo "Blasted $news.filtered.fasta.masked against nt"

echo "End of blastn.contam process"

}


##############################################################################
##																			                                    ##
##			           STEP 10: RUN KRAKEN ON FILTERED READS                    ##
##																			                                    ##	
##############################################################################

kraken ()

{

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to" | tr '\n' ' ' && pwd

# Create a directory for outputs ----
mkdir -p annotation/kraken/masurca
mkdir -p annotation/kraken/megahit

echo "Start of kraken process masurca"

# Load kraken ------
module purge
module load Kraken2/2.0.8-beta-IGB-gcc-4.9.4

# Masurca -----
start=`date +%s`  # record start time 

# Run Kraken ------
kraken2 --use-names --threads 6 --quick   \
--report annotation/kraken/masurca/${news}.masurca_kraken2_report.txt \
--classified-out annotation/kraken/masurca/${news}.masurca_kraken2_classified.fasta \
--unclassified-out annotation/kraken/masurca/${news}.masurca_kraken2_unclassified.fasta \
--db /home/groups/h3abionet/RefGraph/data/kraken2/pluspf_20200919 \
annotation/seqkit/masurca/${news}.masurca.filtered.fasta > annotation/kraken/masurca/${news}.masurca_kraken2_output.txt

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run kraken on [$news.masurca.filtered.fasta]"

echo "End of kraken process masurca"

}

##############################################################################
##																			                                    ##
##								                 MAIN                                     ##
##																		                                      ##
##############################################################################

# Main function runs each step/function of the pipeline separately so that
# user can choose to run steps one at a time.

main()
{
	# Determine whether running full pipeline or single step
	runtype="PARTIAL"
  # runtype="FULL"
  echo ""
	echo "*** RUNNING ${runtype} ANNOTATION PIPELINE ***"
  
  setup
  filter
  #blastdbGRCh38 
  blastnGRCh38
  #blastdbCHM13
  blastnCHM13
  repeatmasker
  quast
  blastncontam
  kraken
  
  }

# Run main function
main

endall=`date +%s`
runtimeall=$( echo "scale=2;$((endall-startall)) / 60" | bc )
echo "It took $runtimeall minutes to complete the annotation pipeline for [$news]"