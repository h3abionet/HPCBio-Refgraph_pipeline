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
##																			##
##			           GENERAL WRAPPER RELATED SCRIPTS	                    ##
##																			##	
##############################################################################

# Set fancy fonts for the help message ------
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`

# Help ------
function HELP {
  echo "${BOLD}Help documentation for the HPCBio UIUC Annotation pipeline${NORM}"
  echo ""
  echo "The following options must be specified:"
  echo "${REV}-d${NORM}   The full path to the main results directory${NORM}"  
  #echo "${REV}-s${NORM}   TRUE for single-end and FALSE for paired-end${NORM}"
  #echo "${REV}-f${NORM}   Adapter/primer sequence to be trimmed from R1${NORM}"
  #echo "${REV}-r${NORM}   Adapter/primer sequence to be trimmed from R2${NORM}"
  #echo "${REV}-l${NORM}   DADA2 truncation length for R1${NORM}"
  #echo "${REV}-n${NORM}   DADA2 truncation length for R2${NORM}"
  echo "${REV}-h${NORM}   Displays this help message without complaining about lack of arguments"
  echo ""
  echo "Example: sbatch annotation.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/"
  exit 1
}


# Check the number of arguments. If none are passed, print message and exit ------
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  echo "Did not pass any arguments, you probably need help"
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
##																			##
##			        STEP 1: LOAD THE NECESSARY MODULES                      ##
##																			##	
##############################################################################
setup ()
{

# ----- Load necessary modules ------ 

# Unload loaded modules -----
module purge

# Load what we need -----
module load seqkit/0.12.1
module load BLAST+/2.10.1-IGB-gcc-8.2.0

echo "Necessary modules are loaded"

}

##############################################################################
##																			##
##				    STEP 2: FILTER THE ASSEMBLY FILES                       ##
##																			##	
##############################################################################

filter ()
{

# ----- Set working directory -----
cd ${OPT_d}/assembly/Final-Assembly/

echo "Input directory is set to" | tr '\n' ' ' && pwd

# Create a directory for outputs ----
mkdir -p ../../annotation/seqkit/masurca
mkdir -p ../../annotation/seqkit/megahit

# ----- Filter sequences < 500bp length -----
echo "Start of filtering process masurca"

# Masurca -----
cd masurca
echo "Now directory is set to" | tr '\n' ' ' && pwd

seqkit seq --min-len 500 --remove-gaps HG03563.masurca.final.fasta > ../../../annotation/seqkit/masurca/HG03563.filtered.fasta
echo "[Now creating stats on filtered fasta files]"
seqkit stats HG03563.masurca.final.fasta >> ../../../annotation/seqkit/masurca/HG03563.stats_before_filter.masurca.txt
seqkit stats ../../../annotation/seqkit/masurca/HG03563.filtered.fasta >> ../../../annotation/seqkit/masurca/HG03563.stats_after_filter.masurca.txt


echo "End of filtering process"

}


##############################################################################
##																			##
##	             STEP 3: CREATE BLAST DATABASE (GRCh38)                     ##
##																			##	
##############################################################################

blastdbGRCh38 ()
{

# ----- Set working directory -----
cd ${OPT_d}

echo "Input directory is set to" | tr '\n' ' ' && pwd

# ----- Creat Blast database from GRC38.decoy.hla -----
makeblastdb -in ../GRCh38 -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

#Building a new DB, current time: 09/03/2021 18:47:07
#New DB name:   /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa
#New DB title:  GRCh38.decoy.hla
#Sequence type: Nucleotide
#Keep MBits: T
#Maximum file size: 1000000000B
#Adding sequences from FASTA; added 3366 sequences in 30.3125 seconds.

}


##############################################################################
##																			                                    ##
##			   STEP 4: RUN BLAST ON FILTERED READS (GRCh38)                     ##
##																			                                    ##	
##############################################################################

blastGRCh38 ()
{

# ----- Run blastn on filtered sequences from previous step ------

# ----- Load necessary modules ------ 
module load BLAST+/2.10.1-IGB-gcc-8.2.0

# Set working directory -----
cd ${OPT_d}/annotation
echo "Directory is set to" | tr '\n' ' ' && pwd

# Make a directory for blast ------
mkdir -p blast/masurca
mkdir -p blast/megahit

echo "Start of BLASTN process"

# Masurca -----
cd seqkit/masurca

# In order for the following codes to run properly without warning, you need to create a hidden file named .ncbirc and
# put it inside seqkit/masurca and seqkit/megahit
# Contents of this file (remove the # from each line):

#; Start the section for BLAST configuration
#[BLAST]
#; Specifies the path where BLAST databases are installed
#BLASTDB=/home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/GRCh38/

# Run blast ------

#  srun -A h3abionet -n 1 -N 1 --mem 20g -p lowmem --pty /bin/bash

echo "RUNNING HG03563.masurca.filtered.fasta"

start=`date +%s`
# Restrict number of hits -----
    blastn -db ../../../../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -query HG03563.masurca.filtered.fasta\
    -out ../../blast/masurca/blast_HG03563.masurca.filtered-limit.asn \
    -outfmt 11 \
    -max_target_seqs 50 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../../blast/masurca/blast_HG03563.masurca.filtered-limit.asn \
-outfmt "7 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" > ../../blast/masurca/blast_HG03563.masurca.filtered-limit.tsv

end=`date +%s`
runtime=$((end-start))
echo "$runtime"


start=`date +%s`
# Without limiting -----
    blastn -db ../../../../GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -query HG03563.masurca.filtered.fasta\
    -out ../../blast/masurca/blast_HG03563.masurca.filtered.asn \
    -outfmt 11 \
    -evalue 1e-5 \
    -perc_identity 90

echo "Convert blast archive .asn to tabular format"

blast_formatter -archive ../../blast/masurca/blast_HG03563.masurca.filtered.asn \
-outfmt "7 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" > ../../blast/masurca/blast_HG03563.masurca.filtered.tsv

end=`date +%s`
runtime=$((end-start))
echo "$runtime"

#   echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,sblastnames,mismatch,gapopen,qstart,qend,qlen sstart send slen" | tr ',' '\t'| cat - ../blast/masurca/blast_HG03563.filtered.0.asn > ../blast/masurca/blast_HG03563.filtered.tsv
  # rm ../../blast/masurca/blast_HG03563.filtered.0.tsv
   echo "Blasted HG03563.filtered.fasta"
   
echo "End of BLASTN process"
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
	runtype="SINGLE STEP"

	echo "*** RUNNING ${runtype} ANNOTATION PIPELINE ***"
    setup
    #filter
    #blastdbGRCh38
    blastGRCh38
}

# Run main function
main
