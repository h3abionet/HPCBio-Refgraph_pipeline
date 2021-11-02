#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome1               = false         /*genome fasta file GRCh38, must specify complete path. Required parameters*/
params.genome2               = false         /*genome fasta file GRCh38.p0, must specify complete path. Required parameters*/
params.genome3               = false         /*genome fasta file, CHM13, must specify complete path. Required parameters*/
params.krakendb              = false         /*kraken database file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = "./results"    /*output folder, must specify complete path. Required parameters*/
params.assembler             = 'masurca'      /*options: megahit|masurca. Default is masurca*/

/* parameters for readprep = seqkit */
params.min_read_length       = '500'          /*minimum length of read to be kept after trimming for downstream analysis. Default is 500*/

/*Stage*/
stage = "annotation"

/*Results path*/
resultsPath = "${params.outputDir}/${stage}"


/*Cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbio'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '6'
assemblerMemory              = '100'
params.clusterAcct           = " -A h3bionet "

/*Prepare input*/
genome_file1                 = file(params.genome1)
genome_file2                  = file(params.genome2)
genome_file3                  = file(params.genome2)
genomeStore1                  = genome_file1.getParent()
genomeStore2                  = genome_file2.getParent()
genomeStore3                  = genome_file3.getParent()
krakendb_file                 = file(params.krakendb)

// Sanity checks
if( !genome_file1.exists() ) exit 1, "Missing reference genome file: ${genome_file1}"
if( !genome_file2.exists() ) exit 1, "Missing reference genome file: ${genome_file2}"
if( !genome_file3.exists() ) exit 1, "Missing reference genome file: ${genome_file3}"
//if( params.assembler != "megahit" || params.assembler != "masurca" ) exit 1, "Unknown assembler: ${params.assembler}"


/* Introduce input files */
fasta_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)


/*
  STEP 0: FILTER THE ASSEMBLY FILES -----
*/
process filter_seqkit {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "seqkit/0.12.1"
    publishDir             "${resultsPath}/seqkit/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(fasta) from fasta_Ch1

    output:
    tuple val(id), file('*.filtered.fasta') into filtered_file
    file "${id}.filter-stats.txt"

    script:
    """
    # Run seqkit to filter below 500 length ------
    seqkit seq --min-len ${params.min_read_length} --remove-gaps ${fasta} > ${id}.filtered.fasta
  
    # Create seqkit stats before and after filtering ------
    seqkit stats ${fasta} > ${id}.filter-stats.txt
    seqkit stats ${id}.filtered.fasta | sed -e '1d' >> ${id}.filter-stats.txt

    """
}

/*
  STEP 1.1: CREATE BLAST DATABASE (GRCh38)  
  This process is executed only once
*/
process blastdbGRCh38 {
    tag                    { "PREP:${genome1}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    storeDir               genomeStore1
    
    input:
    file genome1 from genome_file1

    output:
    file "*" into genome1_db_ch
    
    script:
    """
    makeblastdb -in ${genome1} -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

    """
}

/*
  STEP 1.2: RUN BLAST ON FILTERED READS (GRCh38) -----
*/
process blastnGRCh38 {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/blast/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered) from filtered_file
    file genome1 from genome_file1
    file "*" from genome1_db_ch

    output:
    tuple val(id), file('blast_*.GRCh38.decoy.hla.asn')
    tuple val(id), file('blast_*.GRCh38.decoy.hla.txt')

    script:
    """
    # Run blast ------
    blastn -db ${genome1} \
    -query ${filtered} \
    -out blast_${id}.GRCh38.decoy.hla.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

    # Change format to tabular -----
    blast_formatter -archive blast_${id}.GRCh38.decoy.hla.asn \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
    > blast_${id}.GRCh38.decoy.hla.txt
    
    """
}
