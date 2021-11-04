#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome1               = false         /*genome fasta file GRCh38, must specify complete path. Required parameters*/
params.genome2               = false         /*genome fasta file GRCh38.p0, must specify complete path. Required parameters*/
params.genome3               = false         /*genome fasta file, CHM13, must specify complete path. Required parameters*/
params.krakendb              = false         /*kraken database file, must specify complete path. Required parameters*/
params.samplePath            = false         /*input folder, must specify complete path. Required parameters*/
params.taxdbPath          = false         /*location of blast taxa database. Required parameters*/

/* parameters to be used inside the pipeline */
params.outputDir             = "./results"    /*output folder, must specify path from current directory. Required parameters*/
params.assembler             = 'masurca'      /*options: megahit|masurca. Default is masurca*/

/* parameters for seqkit */
params.min_read_length       = '500'          /*minimum length of read to be kept after trimming with seqkit for downstream analysis. Default is 500*/

/* parameters for blast */
params.max_target_seqs       = '5'            /*number of aligned sequences to keep in blast. Default is 5*/
params.max_hsps              = '10'           /*maximum number of HSPs (alignments) to keep for any single query-subject pair in blast. Default is 10*/
params.evalue                = '1e-5'         /*expect value (E) for saving hits in blast. Default is 1e-5*/
params.perc_identity         = '90'           /*percentage of identical matches in blast. Default is 90*/

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
krakendb_file                = file(params.krakendb)


// Sanity checks
if( !genome_file1.exists() ) exit 1, "Missing reference genome file: ${genome_file1}"
if( !genome_file2.exists() ) exit 1, "Missing reference genome file: ${genome_file2}"
if( !genome_file3.exists() ) exit 1, "Missing reference genome file: ${genome_file3}"
//if( params.assembler != "megahit" || params.assembler != "masurca" ) exit 1, "Unknown assembler: ${params.assembler}"


/* Introduce input files */
fasta_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)
taxdb_Ch1 = Channel.fromFilePairs("${params.taxdbPath}", size: 1)

/*
  STEP 0.1: CREATE BLAST DATABASE (GRCh38)  
  --- This process is executed only once ---
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
    file "*" into blast_db1_ch
    
    script:
    """
    makeblastdb -in ${genome1} -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

    """
}
/*
  STEP 0.2: CREATE BLAST DATABASE (GRCh38.p0)
  --- This process is executed only once ---
*/
process blastdbGRCh38p0 {
    tag                    { "PREP:${genome2}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    storeDir               genomeStore2
    
    input:
    file genome2 from genome_file2

    output:
    file "*" into blast_db2_ch
    
    
    script:
    """
    makeblastdb -in ${genome2} -parse_seqids -title "GRCh38.p0.no.decoy.hla" -dbtype nucl

    """
}
/*
  STEP 0.3: CREATE BLAST DATABASE (CHM13)  
  --- This process is executed only once ---
*/
process blastdbCHM13 {
    tag                    { "PREP:${genome3}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    storeDir               genomeStore3
    
    input:
    file genome3 from genome_file3

    output:
    file "*" into blast_db3_ch
    
    script:
    """
    makeblastdb -in ${genome3} -parse_seqids -title "CHM13.v1.1_GRCh38.p13.chrY" -dbtype nucl

    """
}

/*
  STEP 1: FILTER THE ASSEMBLY FILES 
  --- use seqkit to remove low read lengths ---
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
    tuple val(id), file('*.minLength.*.fasta') into filtered_file1,filtered_file2,filtered_file3,filtered_file4,filtered_file5
    file "${id}.filter-stats.txt"

    script:
    """
    # Run seqkit to filter below 500 length ------
    seqkit seq --min-len ${params.min_read_length} --remove-gaps ${fasta} > ${id}.minLength.${params.min_read_length}.fasta
  
    # Create seqkit stats before and after filtering ------
    seqkit stats ${fasta} > ${id}.filter-stats.txt
    seqkit stats ${id}.minLength.${params.min_read_length}.fasta | sed -e '1d' >> ${id}.filter-stats.txt

    """
}

/*
  STEP 2: RUN KRAKEN ON FILTERED READS 
  --- for removing contaminations ---
*/
process kraken {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "Kraken2/2.0.8-beta-IGB-gcc-4.9.4"
    publishDir             "${resultsPath}/kraken/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered1) from filtered_file1
    file krakendb from krakendb_file

    output:
    tuple val(id),file('*_kraken2_report.txt')
    tuple val(id),file('*_kraken2_classified.fasta')
    tuple val(id),file('*_kraken2_unclassified.fasta')
    tuple val(id),file('*_kraken2_output.txt')
    file "${id}_kraken2_output_filtered.txt"
    file "${id}_to_keep.txt" into keep_list

    script:
    """
    # Run Kraken ------
    kraken2 --use-names --threads 6 --quick   \
    --report ${id}_kraken2_report.txt \
    --classified-out ${id}_kraken2_classified.fasta \
    --unclassified-out ${id}_kraken2_unclassified.fasta \
    --db ${krakendb} \
    ${filtered1} > ${id}_kraken2_output.txt

    # Filter the output to remove contamination ------
    grep -E "Homo sapiens|Eukaryota|cellular organisms|unclassified|root" ${id}_kraken2_output.txt > ${id}_kraken2_output_filtered.txt
    cut -f2 ${id}_kraken2_output_filtered.txt > ${id}_to_keep.txt
    
    """
}

/*
  STEP 3: RUN BLAST ON FILTERED READS (nt)
  --- for removing contaminations ---
*/
process blastncontam {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0","ncbi-blastdb/20201212","seqkit/0.12.1"              
    publishDir             "${resultsPath}/blast-contam/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered2) from filtered_file2
    file keep from keep_list

    output:
    tuple val(id), file('blast_*.nt.asn')
    tuple val(id), file('blast_*.nt.txt')

    script:
    """
    # Filter the fasta by the keep_list ------
    seqkit grep -i -f ${keep} ${filtered2} > ${id}_kn_filtered.fasta 

    # Run blast ------
    blastn -db nt \
    -query ${id}_kn_filtered.fasta \
    -out blast_${id}.nt.asn \
    -outfmt 11 \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.perc_identity} \
    -num_threads ${task.cpus}

    # Change format to tabular -----
    blast_formatter -archive blast_${id}.nt.asn \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
    > blast_${id}.masurca.nt.txt
    
    """
}

/*
  STEP 4: RUN BLAST ON FILTERED READS (GRCh38, GRCh38.p0, CHM13) 

process blastn {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/blast/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered3) from filtered_file3
    file genome1 from genome_file1
    file genome2 from genome_file2
    file genome3 from genome_file3
    file "*" from blast_db1_ch
    file "*" from blast_db2_ch
    file "*" from blast_db3_ch
    file "taxdb.bt?" from taxdb_Ch1

    output:
    tuple val(id), file('blast_*.GRCh38.decoy.hla.asn')
    tuple val(id), file('blast_*.GRCh38.decoy.hla.txt')
    tuple val(id), file('blast_*.GRCh38.p0.no.decoy.hla.asn')
    tuple val(id), file('blast_*.GRCh38.p0.no.decoy.hla.txt')
    tuple val(id), file('blast_*.CHM13.v1.1_GRCh38.p13.chrY.asn')
    tuple val(id), file('blast_*.CHM13.v1.1_GRCh38.p13.chrY.txt')

    script:
    """
    # Run blast (GRCh38) ------
    blastn -db ${genome1} \
    -query ${filtered3} \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore mismatch gapopen qstart qend qlen sstart send slen" \
    -out blast_${id}.GRCh38.decoy.hla.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.perc_identity} \
    -num_threads ${task.cpus}
    """

    script:
    """
    # Run blast (GRCh38.p0) ------
    blastn -db ${genome2} \
    -query ${filtered3} \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore mismatch gapopen qstart qend qlen sstart send slen"
    -out blast_${id}.GRCh38.p0.no.decoy.hla.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.perc_identity} \
    -num_threads ${task.cpus}
    """

   script:
    """
    # Run blast (CHM13) ------
    blastn -db ${genome3} \
    -query ${filtered3} \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore mismatch gapopen qstart qend qlen sstart send slen"
    -out blast_${id}.CHM13.v1.1_GRCh38.p13.chrY.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.perc_identity} \
    -num_threads ${task.cpus} 
    """
   
}

/*
  STEP 5: RUN REPEAT MASKER ON FILTERED READS
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
    tuple val(id), file(filtered4) from filtered_file4

    output:
    file '*'

    script:
    """
    # Run Repeat Masker ------
    RepeatMasker -species human ${filtered4} -nolow

    """
}

/*
  STEP 6: RUN QUAST ON FILTERED READS
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
    tuple val(id), file(filtered5) from filtered_file5

    output:
    file '*'

    script:
    """
    # Run QUAST ------
    quast.py ${filtered5} \
    --threads ${task.cpus}

    """
}
*/

/*
  LIST OF TOOLS USED IN THIS PIPLINE:
  1. blastn        -->  BLAST+/2.10.1
  2. seqkit        -->  seqkit/0.12.1
  2. kraken2       -->  Kraken2/2.0.8
  4. repeatmasker  -->  RepeatMasker/4.1.2
  5. quast         -->  quast/5.0.0
*/

