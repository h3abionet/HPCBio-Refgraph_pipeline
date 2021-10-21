#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome1                = false          /*genome fasta file GRCh38, must specify complete path. Required parameters*/
params.genome2                = false          /*genome fasta file GRCh38.p0, must specify complete path. Required parameters*/
params.genome3                = false          /*genome fasta file, CHM13, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = "./results"    /*output folder, must specify complete path. Required parameters*/
params.assembler             = 'masurca'      /*options: megahit|masurca. Default is masurca*/

/* parameters for readprep = seqkit */
params.min_read_length       = '500'          /*minimum length of read to be kept after trimming for downstream analysis. Default is 500*/

/*Stage*/
stage = "annotation"

resultsPath = "${params.outputDir}/${stage}"


/*cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbio'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '6'
assemblerMemory              = '100'
params.clusterAcct           = " -A h3bionet "

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()

// Sanity checks
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
//if( params.assembler != "megahit" || params.assembler != "masurca" ) exit 1, "Unknown assembler: ${params.assembler}"

CRAM_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)

/*
  prepare_genome 
  This process is executed only once
*/

process prepare_genome{
    tag                    { "PREP:${genome}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    storeDir               genomeStore
    
    input:
    file genome from genome_file

    output:
    file "*.fai" into genome_index_ch
    
    script:
    """
    samtools faidx ${genome}
    """
}