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
    tuple val(id), file('*.filtered.fasta') into filtered_file,filtered_file1,filtered_file2,filtered_file3,filtered_file4,filtered_file5,filtered_file6
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

/*
  STEP 2.1: CREATE BLAST DATABASE (GRCh38.p0)  
  This process is executed only once
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
    file "*" into genome2_db_ch
    
    script:
    """
    makeblastdb -in ${genome2} -parse_seqids -title "GRCh38.p0.no.decoy.hla" -dbtype nucl

    """
}

/*
  STEP 2.2: RUN BLAST ON FILTERED READS (GRCh38.p0) -----
*/
process blastnGRCh38p0 {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/blast/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered1) from filtered_file1
    file genome2 from genome_file2
    file "*" from genome2_db_ch

    output:
    tuple val(id), file('blast_*.GRCh38.p0.no.decoy.hla.asn')
    tuple val(id), file('blast_*.GRCh38.p0.no.decoy.hla.txt')

    script:
    """
    # Run blast ------
    blastn -db ${genome2} \
    -query ${filtered1} \
    -out blast_${id}.GRCh38.p0.no.decoy.hla.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

    # Change format to tabular -----
    blast_formatter -archive blast_${id}.GRCh38.p0.no.decoy.hla.asn \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
    > blast_${id}.GRCh38.p0.no.decoy.hla.txt
    
    """
}

/*
  STEP 3.1: CREATE BLAST DATABASE (CHM13)  
  This process is executed only once
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
    file "*" into genome3_db_ch
    
    script:
    """
    makeblastdb -in ${genome3} -parse_seqids -title "CHM13.v1.1_GRCh38.p13.chrY" -dbtype nucl

    """
}

/*
  STEP 3.2: RUN BLAST ON FILTERED READS (CHM13) -----
*/
process blastnCHM13 {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/blast/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered2) from filtered_file2
    file genome3 from genome_file3
    file "*" from genome3_db_ch

    output:
    tuple val(id), file('blast_*.CHM13.v1.1_GRCh38.p13.chrY.fna.asn')
    tuple val(id), file('blast_*.CHM13.v1.1_GRCh38.p13.chrY.fna.txt')

    script:
    """
    # Run blast ------
    blastn -db ${genome3} \
    -query ${filtered2} \
    -out blast_${id}.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

    # Change format to tabular -----
    blast_formatter -archive blast_${id}.CHM13.v1.1_GRCh38.p13.chrY.fna.asn \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
    > blast_${id}.CHM13.v1.1_GRCh38.p13.chrY.fna.txt
    
    """
}


/*
  STEP 4: RUN REPEAT MASKER ON FILTERED READS -----
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
    tuple val(id), file(filtered3) from filtered_file3

    output:
    file '*'

    script:
    """
    # Run Repeat Masker ------
    RepeatMasker -species human ${filtered3}

    """
}

/*
  STEP 5: RUN QUAST ON FILTERED READS -----
*/
process quast {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             "${resultsPath}/QUAST/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered4) from filtered_file4

    output:
    file '*'

    script:
    """
    # Run QUAST ------
    quast.py ${filtered4} \
    --threads 2 

    """
}

/*
  STEP 6: RUN BLAST ON FILTERED READS -----
*/
process blastncontam {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0"
    module                 "ncbi-blastdb/20201212"
    publishDir             "${resultsPath}/blast-contam/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered5) from filtered_file5

    output:
    tuple val(id), file('blast_*.masurca.nt.asn')
    tuple val(id), file('blast_*.masurca.nt.txt')

    script:
    """
    # Run blast ------
    blastn -db nt \
    -query ${filtered5} \
    -out blast_${id}.masurca.nt.asn \
    -outfmt 11 \
    -max_target_seqs 5 \
    -max_hsps 10 \
    -evalue 1e-5 \
    -perc_identity 90

    # Change format to tabular -----
    blast_formatter -archive blast_${id}.masurca.nt.asn \
    -outfmt "6 qseqid sseqid stitle pident length evalue qcovs bitscore sblastnames mismatch gapopen qstart qend qlen sstart send slen" \
    > blast_${id}.masurca.nt.txt
    
    """
}

/*
  STEP 7: RUN KRAKEN ON FILTERED READS -----
*/
process kraken {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   assemblerCPU
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 "Kraken2/2.0.8-beta-IGB-gcc-4.9.4"
    publishDir             "${resultsPath}/kraken/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered6) from filtered_file6
    file krakendb from krakendb_file

    output:
    tuple val(id),file('*_kraken2_report.txt')
    tuple val(id),file('*_kraken2_classified.fasta')
    tuple val(id),file('*_kraken2_unclassified.fasta')
    tuple val(id),file('*_kraken2_output.txt')

    script:
    """
    # Run Kraken ------
    kraken2 --use-names --threads 6 --quick   \
    --report ${id}_kraken2_report.txt \
    --classified-out ${id}_kraken2_classified.fasta \
    --unclassified-out ${id}_kraken2_unclassified.fasta \
    --db ${krakendb} \
    ${filtered6} > ${id}_kraken2_output.txt
        
    """
}
