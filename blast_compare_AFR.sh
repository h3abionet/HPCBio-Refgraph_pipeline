#!/bin/bash

#SBATCH --mem 50G
#SBATCH --job-name compare-test
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



### HPCBio: This compares test long read insertions with our files -----; Created by Negin Valizadegan Feb 19, 20222; valizad2@illinois.edu

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
  echo "${BOLD}Help Documentation for the HPCBio UIUC RefGraph (Testing) Pipeline${NORM}"
  echo ""
  echo "The Following Options Must Be Specified:"
  echo "${REV}-d${NORM}     The full path to the results directory${NORM} (Required)"
  echo "${REV}-s${NORM}     Path to the subject fasta file [file to be tested against] (Required)"
  echo "${REV}-q${NORM}     Path to the query fasta file [file to be tested] (Required)"
  echo "${REV}-i${NORM}     percentage of identity for blast (Required)"
  echo "${REV}-t${NORM}   number of aligned sequences to keep in blast [max_target_seqs] (Required)"
  echo "${REV}-p${NORM}  maximum number of HSPs (alignments) to keep for any single query-subject pair in blast [max_hsps] (Required)"
  echo "${REV}-e${NORM}  e-value for blast (Required)"
  echo "${REV}-h${NORM}     Displays this help message without complaints (Optional)"
  echo ""
  echo "[ ${NORM}${BOLD}Example:${NORM} sbatch blast_compare_AFR.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/ \
  -s /home/groups/h3abionet/RefGraph/results/2022-01-30-HGSVCv2-truthset/data/insertions.fasta \
  -q assembly/Final-Assembly/masurca/AFR/merged_sequences_AFR_Final_Assembly.fasta \
  -i 100 \
  -t 5 \
  -p 10 \
  -e 1e-5 ]"
  echo ""
  echo "${BOLD}Other Example Queries:${NORM}"
  echo "filter/SeqKit/masurca/AFR/merged_sequences_AFR_Seqkit_500.fasta;"
  echo "filter/BLAST_Contam/masurca/AFR/merged_masurca_kn_filtered.fasta;"
  echo "filter/BLAST_Contam/masurca/AFR/merged_masurca_blst_kn_filtered.fasta;"
  echo "filter/CD-HIT/masurca/AFR/merged_masurca_blst_kn_cdhit.fasta;"
  echo "filter/Final-Filtered/masurca/AFR/merged_masurca_GRCH38_p0_filter.final.fasta;"
  echo "filter/Final-Filtered/masurca/AFR/merged_masurca_GRCH38_decoys_hla_filter.final.fasta"
  echo "filter/Final-Filtered/masurca/AFR/merged_masurca_CHM13_filter.final.fasta;"
  echo "annotation/Merged_Reads/masurca/merged_sequences_GRCH38_p0.fasta;"m
  echo "annotation/Merged_Reads/masurca/merged_sequences_GRCH38_decoys.fasta"
  echo "annotation/Merged_Reads/masurca/merged_sequences_CHM13.fasta;"
  echo "annotation/Cluster_CDHIT/masurca/clustered_GRCH38_p0.fasta;"
  echo "annotation/Cluster_CDHIT/masurca/clustered_GRCH38_decoys.fasta"
  echo "annotation/Cluster_CDHIT/masurca/clustered_CHM13.fasta"
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
while getopts :d:s:q:i:c:t:p:e:h FLAG; do
  case $FLAG in
    d)  #set option "d"
      OPT_d=$OPTARG
      ;;
    s)  #set option "s"
      OPT_s=$OPTARG
      ;; 
    q)  #set option "q"
      OPT_q=$OPTARG
      ;;
    i)  #set option "i"
      OPT_i=$OPTARG
      ;; 
     t)  #set option "t"
      OPT_t=$OPTARG
      ;;
    p)  #set option "p"
      OPT_p=$OPTARG
      ;; 
    e)  #set option "e"
      OPT_e=$OPTARG
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
        echo "No subject sequence path is specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_q" ]]; then
        echo "No  query sequence path is specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_i" ]]; then
        echo "No percent identity is specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_t" ]]; then
        echo "No max_target_seqs is specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_p" ]]; then
        echo "No maximum number of HSPs is specified, aborting script"
        exit 1
fi

if [[ -z "$OPT_e" ]]; then
        echo "No e-value is specified, aborting script"
        exit 1
fi

##############################################################################
##																			##
##				                 STEP 1: SET UP                             ##
##																			##	
##############################################################################

# Load nextflow ------
module load BLAST+/2.10.1-IGB-gcc-8.2.0

echo "blast module loaded"


# Set working directory -------
cd ${OPT_d}
mkdir -p pipeline_testing

echo "output directory is set"


##############################################################################
##																			##
##				           STEP 2: Blast run                                ## 
##																			##	
##############################################################################

# blast to compare two sequences ---
echo "blast is starting"
start=`date +%s` # capture start time 

news=$(echo $(basename ${OPT_q}) | sed 's/.fasta//g')
blastn \
-query ${OPT_q} \
-subject ${OPT_s} \
 -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
-perc_identity ${OPT_i} \
-max_target_seqs ${OPT_t} \
-max_hsps ${OPT_p} \
-evalue ${OPT_e} \
-out pipeline_testing/${news}_blast_testdata.0.txt

echo "blast has ended"

# Create headers for the blast file -----
echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
qlen,sstart,send,slen" | tr ',' '\t'| cat - pipeline_testing/${news}_blast_testdata.0.txt > pipeline_testing/${news}_blast_testdata.txt

rm pipeline_testing/${news}_blast_testdata.0.txt

echo "headers are added to blast output"

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run maker on ${OPT_s}"


