#!/usr/bin/env nextflow

/*
  LIST OF TOOLS USED IN THIS PIPLINE:
  1. repeatmasker  -->  RepeatMasker/4.1.2
  2. quast         -->  quast/5.0.0
  3. cd-hit        -->  CD-HIT/4.8.1
*/

/*Parameters that are specified at the command line or via config file*/
params.genome1                = false          /*genome fasta file GRCh38, must specify complete path. Required parameter*/
params.genome2                = false          /*genome fasta file GRCh38.p0, must specify complete path. Required parameter*/
params.genome3                = false          /*genome fasta file, CHM13, must specify complete path. Required parameter*/
params.samplePath             = false          /*input folder, must specify complete path. Required parameter*/

/*Parameters to be used inside the pipeline */
params.outputDir              = "./results"    /*output folder, must specify path from current directory. Required parameter*/
params.assembler              = 'masurca'      /*options: megahit|masurca. Default is masurca*/

/*Parameters for cdhit */
params.cdhit_identity         = '0.9'          /*proportion of idenitity for clustering using cdhit. Default is 0.9*/
params.cdhit_wordsize         = '7'            /*word size for cdhit. Default is 7*/

/*Stage*/
stage = "annotation"

/*Results path*/
resultsPath = "${params.outputDir}/${stage}"

/*Cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbio'
defaultCPU                   = '9'
defaultMemory                = '120'
params.clusterAcct           = " -A h3bionet "

/*Prepare input*/
genome_file1                 = file(params.genome1)
genome_file2                 = file(params.genome2)
genome_file3                 = file(params.genome3)
genomeStore1                 = genome_file1.getParent()
genomeStore2                 = genome_file2.getParent()
genomeStore3                 = genome_file3.getParent()

// Sanity checks
if( !genome_file1.exists() ) exit 1, "Missing reference genome file: ${genome_file1}"
if( !genome_file2.exists() ) exit 1, "Missing reference genome file: ${genome_file2}"
if( !genome_file3.exists() ) exit 1, "Missing reference genome file: ${genome_file3}"
//if( params.assembler != "megahit" || params.assembler != "masurca" ) exit 1, "Unknown assembler: ${params.assembler}"

/* Create channcel for input files */
filtered_fasta_Ch = Channel.fromFilePairs("${params.samplePath}", size: 1)
filtered_fasta_Ch2 = Channel.fromFilePairs("${params.samplePath}", size: 1)
filtered_fasta_Ch3 = Channel.fromFilePairs("${params.samplePath}", size: 1)

/*
  STEP 1: RUN REPEAT MASKER ON FILTERED READS
*/
process repeatmasker {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "RepeatMasker/4.1.2-p1-IGB-gcc-8.2.0-Perl-5.28.1"
    publishDir             "${resultsPath}/RepeatMasker/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(fasta) from filtered_fasta_Ch 

    output:
    tuple val(id), file('*.fasta.cat')
    tuple val(id), file('*.fasta.masked') 
    tuple val(id), file('*.fasta.out') 
    tuple val(id), file('*.fasta.tbl')

    script:
    """
    # Run Repeat Masker ------
    RepeatMasker -species human ${fasta} -nolow

    """
}

/*
  STEP 2: RUN QUAST ON FILTERED READS
*/
process quast {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             "${resultsPath}/QUAST/${id}/",mode:"copy"

    input:
    tuple val(id), file(fasta2) from filtered_fasta_Ch2

    output:
    file '*'

    script:
    """
    # Run QUAST ------
    quast.py ${fasta2} \
    --threads ${task.cpus}

    """
}

/*
  STEP 3: MERGE ALL SEQUENCES FROM ALL SAMPLES
*/
process merge_reads {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Merged_Reads/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(fasta3) from filtered_fasta_Ch3

    output:
    file('merged_sequences.fasta') into merged_seqs

    script:
    """
    # Combine all sequences ------
    cat ${fasta3} >> merged_sequences.fasta

    """
}

/*
  STEP 4: RUN CD-HIT ON FILTERED READS
*/
process cdhit {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "CD-HIT/4.8.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Cluster_CDHIT/${params.assembler}/",mode:"copy"

    input:
    file(merged) from merged_seqs

    output:
    file('clustered.fasta') into clustered

    script:
    """
    # Use cd-hit to cluster and remove redundancy ------
    cd-hit-est \
    -i ${merged} \
    -o clustered.fasta \
    -c ${params.cdhit_identity} \
    -n ${params.cdhit_wordsize} \
    -T ${task.cpus}

    """
}