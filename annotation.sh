#!/bin/bash

#SBATCH --mem 70G
#SBATCH --job-name annotation
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH --output slurm-%j.out
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



# HPCBio UIUC Annotation pipeline; Created by Negin Valizadegan Sep 2, 2021; valizad2@illinois.edu

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
  echo "${BOLD}Help Documentation for the HPCBio UIUC Annotation Pipeline${NORM}"
  echo ""
  echo "The Following Options Must Be Specified:"
  echo "${REV}-d${NORM}   The full path to the main results directory${NORM} (Required)"  
  #echo "${REV}-s${NORM}   TRUE for single-end and FALSE for paired-end${NORM}"
  #echo "${REV}-f${NORM}   Adapter/primer sequence to be trimmed from R1${NORM}"
  #echo "${REV}-r${NORM}   Adapter/primer sequence to be trimmed from R2${NORM}"
  #echo "${REV}-l${NORM}   DADA2 truncation length for R1${NORM}"
  #echo "${REV}-n${NORM}   DADA2 truncation length for R2${NORM}"
  echo "${REV}-h${NORM}   Displays this help message without complaints (Optional)"
  echo ""
  echo "[ ${NORM}${BOLD}Example:${NORM} sbatch annotation.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/ ]"
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
while getopts :d:s:f:r:l:n:h FLAG; do
  case $FLAG in
    d)  #set option "d"
      OPT_d=$OPTARG
      ;;
    #s)  #set option "s"
     # OPT_s=$OPTARG
      #;;
    #f)  #set option "f"
     # OPT_f=$OPTARG
      #;;
   # r)  #set option "r"
    #  OPT_r=$OPTARG
     # ;;
    #l)  #set option "l"
     # OPT_l=$OPTARG
      #;;
    #n)  #set option "n"
     # OPT_n=$OPTARG
      #;;
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

#if [[ -z "$OPT_s" ]]; then
#	echo "Single-end TRUE or FALSE not specified, aborting script"
#	exit 1
#fi

#if [[ -z "$OPT_f" ]]; then
#	echo "No adapter sequence for read 1 specified, aborting script"
#	exit 1
#fi

#if [[ -z "$OPT_l" ]]; then
#	echo "No truncation length for read 1 specified, aborting script"
#	exit 1
#fi


##############################################################################
##																			                                    ##
##			           STEP 1: LOAD THE NECESSARY MODULES                       ##
##																			                                    ##	
##############################################################################
setup ()
{

# ----- Load necessary modules ------ 

# Unload loaded modules -----
module purge

# Load what we need -----
module load seqkit/0.12.1
module load BLAST+/2.10.1-IGB-gcc-8.2.0

echo "- Necessary modules are loaded"

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

# ----- Filter sequences < 500bp length -----
echo "Start of filtering process masurca"

# Masurca -----
cd masurca

echo "[Removing reads with length < 500]"

seqkit seq --min-len 500 --remove-gaps assembly/Final-Assembly/HG03563.masurca.final.fasta > annotation/seqkit/masurca/HG03563.filtered.fasta

echo "[Creating stats on fasta files]"

seqkit stats assembly/Final-Assembly/HG03563.masurca.final.fasta >> annotation/seqkit/masurca/HG03563.stats_before_filter.masurca.txt
seqkit stats annotation/seqkit/masurca/HG03563.filtered.fasta >> annotation/seqkit/masurca/HG03563.stats_after_filter.masurca.txt

echo "End of Filtering Process"

}


##############################################################################
##																			                                    ##
##	               STEP 3: CREATE BLAST DATABASE (GRCh38)                   ##
##																			                                    ##	
##############################################################################

blastdb.GRCh38 ()
{

# ----- Set working directory -----
cd ${OPT_d}

echo "Input directory is set to" | tr '\n' ' ' && pwd


# ----- Creat Blast database from GRCh38.decoy.hla [?] -----

makeblastdb -in ../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

# Building a new DB, current time: 09/13/2021 18:24:45
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
# New DB title:  GRCh38.decoy.hla
# Sequence type: Nucleotide
# Deleted existing Nucleotide BLAST database named /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 3366 sequences in 29.3228 seconds.


# ----- Creat Blast database from GRCh38.no.decoy.hla [GRCh38.p0; GCA_000001405.15] -----

# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
makeblastdb -in ../GRCh38.p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -parse_seqids -title "GRCh38.p0.no.decoy.hla" -dbtype nucl

# Building a new DB, current time: 09/13/2021 18:25:15
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38.p0/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
# New DB title:  GRCh38.p0.no.decoy.hla
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 195 sequences in 28.1475 seconds.

}


##############################################################################
##																			                                    ##
##			        STEP 4: RUN BLAST ON FILTERED READS (GRCh38)                ##
##																			                                    ##	
##############################################################################

blastn.GRCh38 ()
{

# ----- Run blastn on filtered sequences from previous step ------

# ----- Load necessary modules ------ 
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to" | tr '\n' ' ' && pwd

# Make a directory for blast ------
mkdir -p annotation/blast/masurca
mkdir -p annotation/blast/megahit

echo "[Start of BLASTN:GRCh38 Process]"

# Masurca -----

# Set directory to GRCh38 ------
cd ../GRCh38/


# BLAST against GRCh38.decoy.hla ------

echo "[Running BLAST+ on HG03563.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run blast ------
   blastn -db GRCh38_full_analysis_set_plus_decoy_hla.fa \
   -query ../results/annotation/seqkit/masurca/HG03563.masurca.filtered.fasta \
   -out ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.decoy.hla.asn \
   -outfmt 11 \
   -max_target_seqs 5 \
   -max_hsps 10 \
   -evalue 1e-5 \
   -perc_identity 90

echo "[Convert BLAST Archive (.asn) to Tabular Format]"

#blast_formatter -archive ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.decoy.hla.asn \
#-outfmt "7 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
#> ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.decoy.hla.tsv

end=`date +%s` # record finish time 
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [HG03563.filtered.fasta] against [GRCh38.decoy.hla]"

echo "Blasted HG03563.filtered.fasta against GRCh38.decoy.hla"


# BLAST against GRCh38.p0.no.decoy.hla ------
start=`date +%s`

# Download the ref sequence ------
# mkdir GRCh38.p0
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Set directory to GRCh38.p0 ------
cd ../GRCh38.p0/

# Run blast ------
    blastn -db GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -query ../results/annotation/seqkit/masurca/HG03563.masurca.filtered.fasta \
    -out ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.p0.no.decoy.hla.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.p0.no.decoy.hla.asn \
-outfmt "7 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> ../results/annotation/blast/masurca/blast_HG03563.masurca.GRCh38.p0.no.decoy.hla.tsv

end=`date +%s`
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [HG03563.filtered.fasta] against [GRCh38.p0.no.decoy.hla]"

echo "Blasted HG03563.filtered.fasta against GRCh38.p0.no.decoy.hla"

echo "End of BLASTN GRCh38 process"

}


##############################################################################
##																			                                    ##
##	                 STEP 5: CREATE BLAST DATABASE (CHM13)                  ##
##																		                                    	##	
##############################################################################

blastdb.CHM13()
{

# ----- Creat Blast database from CHM13 [GCA_009914755.3_CHM13_T2T_v1.1] + chrY of GRCh38.p13 [GCA_000001405.28_GRCh38.p13] -----

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.3_CHM13_T2T_v1.1/GCA_009914755.3_CHM13_T2T_v1.1_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chrY.fna.gz
makeblastdb -in ../CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna -parse_seqids -title "CHM13.v1.1_GRCh38.p13.chrY" -dbtype nucl

# Building a new DB, current time: 09/15/2021 17:21:25
# New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna
# New DB title:  CHM13.v1.1_GRCh38.p13.chrY
# Sequence type: Nucleotide
# Deleted existing Nucleotide BLAST database named /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/CHM13.v1.1_GRCh38.p13.chrY/CHM13.v1.1_GRCh38.p13.chrY.fna
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 25 sequences in 25.1712 seconds.

}

##############################################################################
##																			                                    ##
##			        STEP 6: RUN BLAST ON FILTERED READS (CHM13)                 ##
##																			                                    ##	
##############################################################################

blastn.CHM13 ()
{

# ----- Run blastn on filtered sequences from previous step ------

# ----- Load necessary modules ------ 
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# Set working directory -----
cd ${OPT_d}
echo "Directory is set to:"  |  tr '\n' ' ' && pwd 

# Make a directory for blast ------
mkdir -p annotation/blast/masurca
mkdir -p annotation/blast/megahit

echo "[Start of BLASTN:CHM13 Process]"

# Masurca -----

# Set directory to CHM13.v1.1_GRCh38.p13.chrY ------
cd ../CHM13.v1.1_GRCh38.p13.chrY/


# BLAST against GRCh38.decoy.hla ------

echo "[Running BLAST+ on HG03563.masurca.filtered.fasta]"

start=`date +%s`  # record start time 

# Run blast ------
    blastn -db CHM13.v1.1_GRCh38.p13.chrY.fna\
    -query ../results/annotation/seqkit/masurca/HG03563.masurca.filtered.fasta \
    -out ../results/annotation/blast/masurca/blast_HG03563.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../results/annotation/blast/masurca/blast_HG03563.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
-outfmt "7 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
> ../results/annotation/blast/masurca/blast_HG03563.masurca.CHM13.v1.1_GRCh38.p13.chrY.fna.tsv

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to blast [HG03563.filtered.fasta] against [CHM13.v1.1_GRCh38.p13.chrY]"

echo "Blasted HG03563.filtered.fasta against CHM13.v1.1_GRCh38.p13.chrY.fna"

echo "End of BLASTN CHM13 process"

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
  #filter
  # blastdb.GRCh38
  # blastn.GRCh38
  # blastdb.CHM13
   blastn.CHM13
}

# Run main function
main