#!/bin/bash

#SBATCH --mem 50G
#SBATCH --job-name maker
#SBATCH --mail-user valizad2@illinois.edu ## CHANGE THIS TO YOUR EMAIL
#SBATCH --mail-type ALL
#SBATCH -n 24
#SBATCH -N 1
#SBATCH -A h3abionet
#SBATCH -o /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/slurm_output/slurm-%j.out



# HPCBio UIUC Gene Annotation pipeline (MAKER + EVidenceModeler); Created by Negin Valizadegan Jan 18, 2022; valizad2@illinois.edu

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
  echo "[ ${NORM}${BOLD}Example:${NORM} sbatch MAKER.sh -d /home/groups/h3abionet/RefGraph/results/NeginV_Test_Summer2021/results/ -s clustered_GRCH38_p0.fasta ]"
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

##############################################################################
##																			                                    ##
##				                   STEP 0: LOAD MODULES                           ##
##																			                                    ##	
##############################################################################

setup ()
{

# Load modules ------
module load MAKER/3.01.03-IGB-gcc-4.9.4-Perl-5.26.1-unthreaded

# Create 3 control files needed for maker ----- (do not run if present; control files should be edited manually)
# cd ${OPT_d}/../HPCBio-Refgraph_pipeline/
# maker -CITL

echo "Control ctl files are created if not already exist. They are usually needed to be manually modified."

}


##############################################################################
##																			                                    ##
##				                    STEP 1: RUN MAKER                             ##
##																			                                    ##	
##############################################################################

maker ()
{

# Set working directory -----
cd ${OPT_d}/annotation

# Create output directory
mkdir -p MAKER
cd MAKER

echo "Working directory is set to" | tr '\n' ' ' && pwd

# Create a temp directory ------
mkdir -p /scratch/valizad2/maker # change valizad2 to your username

start=`date +%s` # capture start time 
echo "Start of maker annotation"

#export AUGUSTUS_CONFIG_PATH=/home/n-z/valizad2/NeginV_Test_Summer2021/augustus/3.2.3-IGB-gcc-4.9.4/config export PATH=$PATH:=/home/n-z/valizad2/NeginV_Test_Summer2021/augustus/3.2.3-IGB-gcc-4.9.4/bin

# Run maker -----
mpiexec -n $SLURM_NPROCS maker \
${OPT_d}/../HPCBio-Refgraph_pipeline/maker_opts.ctl \
${OPT_d}/../HPCBio-Refgraph_pipeline/maker_bopts.ctl \
${OPT_d}/../HPCBio-Refgraph_pipeline/maker_exe.ctl \
-genome ${OPT_d}/annotation/Cluster_CDHIT/masurca/${OPT_s} \
-fix_nucleotides # This will change Ys to Ns

echo "Maker gene prediction is completed for ${OPT_s}"

end=`date +%s`
runtime=$((end-start))
runtime=$( echo "scale=2;$((end-start)) / 60" | bc )
echo "It took $runtime minutes to run maker on ${OPT_s}"

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
	#runtype="PARTIAL"
  runtype="FULL"
  echo ""
	echo "*** RUNNING ${runtype} ANNOTATION PIPELINE ***"
  
  setup
  maker
  
  }


# Run main function
main