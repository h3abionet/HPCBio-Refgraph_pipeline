#!/usr/bin/env nextflow

/*
  LIST OF TOOLS USED IN THIS PIPLINE:
  1. repeatmasker  -->  RepeatMasker/4.1.2
  2. quast         -->  quast/5.0.0
  3. cd-hit        -->  CD-HIT/4.8.1
*/

/*Parameters that are specified at the command line or via config file*/
params.samplePath             = false          /*input folder, must specify complete path. Required parameter*/
params.samplePath1            = false          /*input fasta files path for GRCh38 + decoy + alt, must specify complete path. Required parameter*/
params.samplePath2            = false          /*input fasta files path for GRCh38.p0, must specify complete path. Required parameter*/
params.samplePath3            = false          /*input fasta files path for CHM13, must specify complete path. Required parameter*/
params.skipRepeatMasker       = false          /* If set to true in config file, RepeatMasker would not be run. */

/*Parameters to be used inside the pipeline */
params.outputDir              = "./results"    /*output folder, must specify path from current directory. Required parameter*/
params.assembler              = 'megahit'      /*options: megahit|masurca. Default is masurca*/

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

/* Create channcel for input files */
Channel.fromFilePairs("${params.samplePath}", size: 1).into { filtered_for_RepMaster;filtered_for_quast }
filtered_GRCH38_decoys_Ch3 = Channel.fromPath("${params.samplePath1}").map { file -> tuple(file) }
filtered_GRCH38_p0_Ch3 = Channel.fromPath("${params.samplePath2}").map { file -> tuple(file) }
filtered_CHM13_Ch3 = Channel.fromPath("${params.samplePath3}").map { file -> tuple(file) }


if (params.skipRepeatMasker == 'false') {
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
    tuple val(id), file(fasta) from filtered_for_RepMaster 

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

}

/*
  STEP 2-1: RUN QUAST ON FILTERED READS
*/
process quast {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                  "quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             "${resultsPath}/QUAST/${id}/",mode:"copy",overwrite: true

    input:
    tuple val(id), file(fasta2) from filtered_for_quast

    output:
    file '*' into multiqc_ch

    script:
    """
    # Run QUAST ------
    quast.py ${fasta2} \
    --threads ${task.cpus}
    """
}

/*
  STEP 2-2: RUN MultiQC on QUAST
*/
process multiqc {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "MultiQC/1.11-IGB-gcc-8.2.0-Python-3.7.2"
    publishDir             "${resultsPath}/MultiQC/",mode:"copy",overwrite: true

    input:
    file('./QUAST/*/*') from multiqc_ch.collect().ifEmpty([])

    output:
    file "multiqc*" 

    script:
    """
    # Run MultiQC ------
    multiqc .

    """
}

/*
  STEP 3-1: MERGE ALL SEQUENCES FROM ALL SAMPLES 
  (GRCH38_decoys)
*/
process merge_reads_GRCH38_decoys {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Merged_Reads/${params.assembler}/",mode:"copy",overwrite: true

    input:
    file fasta3 from filtered_GRCH38_decoys_Ch3.collect()

    output:
    file('merged_sequences_GRCH38_decoys.fasta') into merged_seqs_GRCH38_decoys

    script:
    """
    # Combine all sequences ------
    cat ${fasta3} >> merged_sequences_GRCH38_decoys.fasta

    """
}

/*
  STEP 3-2: MERGE ALL SEQUENCES FROM ALL SAMPLES
  (GRCH38_p0)
*/
process merge_reads_GRCH38_p0 {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Merged_Reads/${params.assembler}/",mode:"copy",overwrite: true

    input:
    file fasta33 from filtered_GRCH38_p0_Ch3.collect()

    output:
    file('merged_sequences_GRCH38_p0.fasta') into merged_seqs_GRCH38_p0

    script:
    """
    # Combine all sequences ------
    cat ${fasta33} >> merged_sequences_GRCH38_p0.fasta

    """
}

/*
  STEP 3-3: MERGE ALL SEQUENCES FROM ALL SAMPLES
  (CHM13)
*/
process merge_reads_CHM13 {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Merged_Reads/${params.assembler}/",mode:"copy",overwrite: true

    input:
    file fasta333 from filtered_CHM13_Ch3.collect()

    output:
    file('merged_sequences_CHM13.fasta') into merged_seqs_CHM13

    script:
    """
    # Combine all sequences ------
    cat ${fasta333} >> merged_sequences_CHM13.fasta

    """
}

/*
  STEP 4: RUN CD-HIT ON FILTERED READS
*/
process cdhit {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "CD-HIT/4.8.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Cluster_CDHIT/${params.assembler}/",mode:"copy"

    input:
    file merged_GRCH38_decoys from merged_seqs_GRCH38_decoys
    file merged_GRCH38_p0 from merged_seqs_GRCH38_p0
    file merged_CHM13 from merged_seqs_CHM13

    output:
    file('clustered_GRCH38_decoys.fasta') into clustered_GRCH38_decoys
    file('clustered_GRCH38_p0.fasta') into clustered_GRCH38_p0
    file('clustered_CHM13.fasta') into clustered_CHM13

    script:
    """
    # Use cd-hit to cluster and remove redundancy ------

    # GRCH38_decoys -----
    cd-hit-est \
    -i ${merged_GRCH38_decoys} \
    -o clustered_GRCH38_decoys.fasta \
    -c ${params.cdhit_identity} \
    -n ${params.cdhit_wordsize} \
    -T ${task.cpus}

    # GRCH38_p0 -----
    cd-hit-est \
    -i ${merged_GRCH38_p0} \
    -o clustered_GRCH38_p0.fasta \
    -c ${params.cdhit_identity} \
    -n ${params.cdhit_wordsize} \
    -T ${task.cpus}

    # CHM13 -----
    cd-hit-est \
    -i ${merged_CHM13} \
    -o clustered_CHM13.fasta \
    -c ${params.cdhit_identity} \
    -n ${params.cdhit_wordsize} \
    -T ${task.cpus}

    """
}