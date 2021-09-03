#!/bin/bash

#SBATCH --mem 18G
#SBATCH --job-name annotation
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH --output slurm-%j.out
#SBATCH -n 2
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
  echo "${REV}-d${NORM}   The full path to the Final Assembly directory${NORM}"  
  #echo "${REV}-s${NORM}   TRUE for single-end and FALSE for paired-end${NORM}"
  #echo "${REV}-f${NORM}   Adapter/primer sequence to be trimmed from R1${NORM}"
  #echo "${REV}-r${NORM}   Adapter/primer sequence to be trimmed from R2${NORM}"
  #echo "${REV}-l${NORM}   DADA2 truncation length for R1${NORM}"
  #echo "${REV}-n${NORM}   DADA2 truncation length for R2${NORM}"
  echo "${REV}-h${NORM}   Displays this help message without complaining about lack of arguments"
  echo ""
  echo "Example: sh annotation.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/assembly/Final-Assembly/"
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
##			 STEP 1: SET UP DIRECTORY AND LOAD MODULES                      ##
##																			##	
##############################################################################

setup()
{
# Set working directory -----
# command + / for commenting
cd ${OPT_d}

echo "Input directory is set to" | tr '\n' ' ' && pwd

# Create a directory for outputs ----
mkdir -p ../../annotation/seqkit/masurca
mkdir -p ../../annotation/seqkit/megahit

# Load necessary modules ------ 
module purge
module load seqkit/0.12.1

echo "Necessary modules are loaded"
}


##############################################################################
##																			##
##				    STEP 2: FILTER THE ASSEMBLY FILES                       ##
##																			##	
##############################################################################

filter ()
{

# Filter sequences < 500bp length -----

# Masurca -----
cd masurca
for i in *.masurca.final.fasta
do id=$(echo $i | sed 's/.final.fasta//1')
    echo "Now filtering ${i}"
    seqkit seq --min-len 500 ${i} > ../../../annotation/seqkit/masurca/$id.filtered.fasta
done


# Megahit -----
cd ../megahit
for i in *.megahit.final.fasta
do id=$(echo $i | sed 's/.final.fasta//1')
    echo "Now reading ${i}"
    seqkit seq --min-len 500 ${i} > ../../../annotation/seqkit/megahit/$id.filtered.fasta
done

# Reset to assembly directory
cd ..
pwd

}

##############################################################################
##																			##
##								 MAIN                                       ##
##																		    ##	
##############################################################################

# Main function runs each step/function of the pipeline separately so that
# user can choose to run steps one at a time.

main()
{
	# Determine whether running full pipeline or single step
	runtype="FULL"

	echo "*** RUNNING ${runtype} Annotation Pipeline ***"
	setup
    filter
}

# Run main function
main
