#!/usr/bin/env nextflow

/*
  LIST OF TOOLS USED IN THIS PIPLINE:
  1. blastn        -->  BLAST+/2.10.1
  2. seqkit        -->  seqkit/0.12.1
  3. kraken2       -->  Kraken2/2.0.8
  4. cd-hit        -->  CD-HIT/4.8.1
*/

/*Parameters that are specified at the command line or via config file*/
params.genome1                = false          /*genome fasta file GRCh38, must specify complete path. Required parameter*/
params.genome2                = false          /*genome fasta file GRCh38.p0, must specify complete path. Required parameter*/
params.genome3                = false          /*genome fasta file, CHM13, must specify complete path. Required parameter*/
params.krakendb               = false          /*kraken database file, must specify complete path. Required parameter*/
params.samplePath             = false          /*input folder, must specify complete path. Required parameter*/
params.taxdbPath              = false          /*location of blast taxa database. Required parameter*/
params.skipcdhit              = false          /* If set to true in config file, cdhit would not be run. */

/*Parameters to be used inside the pipeline */
params.outputDir              = "./results/bowtie2"    /*output folder, must specify path from current directory. Required parameter*/
params.assembler              = 'masurca'      /*options: megahit|masurca. Default is masurca*/

/*Parameters for seqkit */
params.min_read_length        = '500'          /*minimum length of read to be kept after trimming with seqkit for downstream analysis. Default is 500*/

/*Parameters for blast */
params.max_target_seqs        = '5'            /*number of aligned sequences to keep in blast. Default is 5*/
params.max_hsps               = '10'           /*maximum number of HSPs (alignments) to keep for any single query-subject pair in blast. Default is 10*/
params.evalue                 = '1e-5'         /*expect value (E) for saving hits in blast. Default is 1e-5*/
params.blastnt_pident         = '60'           /*percentage of identical matches in blast NT. Default is 60*/
params.blastnt_filter_pident  = '60'           /*filtering cut off for percentage of identical matches from blast NT. Default is 60*/
params.blastnt_filter_length  = '100'          /*filtering cut off for alignment length from blast NT. Default is 100*/
params.blastr_pident          = '90'           /*percentage of identical matches in blast human reference. Default is 90*/
params.blastr_filter_pident   = '95'           /*filtering cut off for percentage of identical matches from blast ref genome. Default is 95*/
params.blastr_filter_qcov     = '95'           /*filtering cut off for query coverage from blast ref genome. Default is 95*/

/*Parameters for cdhit */
params.cdhit_identity         = '0.99'          /*proportion of idenitity for clustering. Default is 0.99*/
params.cdhit_wordsize         = '10'            /*word size for cdhit. Default is 5*/

/*Stage*/
stage = "filter"

/*Results path*/
resultsPath = "${params.outputDir}/${stage}"

/*Cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbio'
defaultCPU                   = '24'
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

/* Create channcel for input files */
fasta_Ch1 = Channel.fromFilePairs("${params.samplePath}", size: 1)
taxdb_Ch1 = Channel.fromFilePairs("${params.taxdbPath}", size: 1)


/*
PREPATE DATABASE FILES: STEPS 0.1, 0.2, & 0.3
/*
  STEP 0.1: CREATE BLAST DATABASE (GRCh38)  
*/
process blastdbGRCh38 {
    tag                    { "PREP:${genome1}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.13.0-IGB-gcc-8.2.0"
    storeDir               genomeStore1
     
    input:
    file genome1 from genome_file1

    output:
    file "*" into blast_db1_ch
    
    script:
    """
    start=`date +%s` # capture start time 
    echo "Start of maker annotation"
    makeblastdb -in ${genome1} -parse_seqids -title "GRCh38.decoy.hla" -dbtype nucl

    """
}
/*
  STEP 0.2: CREATE BLAST DATABASE (GRCh38.p0)
*/
process blastdbGRCh38p0 {
    tag                    { "PREP:${genome2}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.13.0-IGB-gcc-8.2.0"
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
*/
process blastdbCHM13 {
    tag                    { "PREP:${genome3}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.13.0-IGB-gcc-8.2.0"
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
  STEP 1: FILTER BASED ON READ LENGTH
/*
  1.1 FILTER THE ASSEMBLY FILES 
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
    publishDir             "${resultsPath}/SeqKit/${params.assembler}/",mode:"copy",overwrite: true

    input:
    tuple val(id), file(fasta) from fasta_Ch1

    output:
    tuple val(id), file("*_minLength_${params.min_read_length}.fasta") into filtered_seqkit_for_kraken
    file "${id}_filter_stats.txt" into all_stats

    script:
    """
    # Run seqkit to filter below 500 length ------
    seqkit seq --min-len ${params.min_read_length} --remove-gaps ${fasta} > ${id}_minLength_${params.min_read_length}.fasta
  
    # Create seqkit stats before and after filtering ------
    seqkit stats ${fasta} > ${id}_filter_stats.txt
    seqkit stats ${id}_minLength_${params.min_read_length}.fasta | sed -e '1d' >> ${id}_filter_stats.txt
    
    """   
}

process merge_seqkit_stats {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/SeqKit/${params.assembler}/",mode:"copy",overwrite: true
  
  input:
  file all_seqkit_stats from all_stats.collect()
  
  output:
  file "All_seqkit_stats0.txt" optional true
  file "All_seqkit_stats.txt"
  
  script:
  """
  sed '1~3d' ${all_seqkit_stats} | cat >> All_seqkit_stats0.txt
   
  # Add header -----
  echo -e "file,format,type,num_seqs,sum_len,min_len,avg_len,max_len" \
  | tr ',' '\t'| cat - All_seqkit_stats0.txt > All_seqkit_stats.txt
  rm All_seqkit_stats0.txt

  """
}


/*
  STEP 2: CONTAMINATION REMOVAL
/*
  2.1: RUN KRAKEN ON FILTERED READS 
*/
process kraken2 {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "Kraken2/2.0.8-beta-IGB-gcc-4.9.4","seqkit/0.12.1" 
    publishDir             "${resultsPath}/Kraken2/${params.assembler}/",mode:"copy",overwrite: true

    input:
    tuple val(id), file(seqkit_filtered) from filtered_seqkit_for_kraken
    file krakendb from krakendb_file

    output:
    file "*_kraken2_report.txt"
    file "*_kraken2_classified.fasta"
    file "*_kraken2_unclassified.fasta"
    file "*_kraken2_output.txt"
    tuple val(id),file('*_kraken2_output_filtered.txt')
    tuple val(id),file('*_to_keep.txt') into keep_list_kn
    tuple val(id), file('*_seqkit_kn_filtered.fasta') into seqkit_kn_filtered

    script:
    """
    # Run Kraken ------
    kraken2 --use-names --threads 6 --quick   \
    --report ${id}_kraken2_report.txt \
    --classified-out ${id}_kraken2_classified.fasta \
    --unclassified-out ${id}_kraken2_unclassified.fasta \
    --db ${krakendb} \
    ${seqkit_filtered} > ${id}_kraken2_output.txt

    # Filter the output to remove contamination ------
    grep -E "Homo sapiens|Eukaryota|cellular organisms|unclassified|root" ${id}_kraken2_output.txt > ${id}_kraken2_output_filtered.txt
    cut -f2 ${id}_kraken2_output_filtered.txt > ${id}_to_keep.txt  

    # Filter the fasta by the keep_list ------
    seqkit grep -i -f ${id}_to_keep.txt ${seqkit_filtered} > ${id}_seqkit_kn_filtered.fasta

    """
}

/*
  2.2: RUN BLAST (NT) ON FILTERED READS 
*/
process blastcontam {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.13.0-IGB-gcc-8.2.0","ncbi-blastdb/20220318","seqkit/0.12.1"              
    publishDir             "${resultsPath}/BLAST_Contam/${params.assembler}/",mode:"copy",overwrite: true


    input:
    tuple val(id), file(kraken_filtered) from seqkit_kn_filtered

    output:
    tuple val(id), file('blast_*_nt0.txt') optional true
    tuple val(id), file('blast_*_nt.txt')
    tuple val(id), file('*_IDs_to_remove_nt.txt')
    tuple val(id), file('*_seqkit_kn_blast_filtered.fasta') into seqkit_blast_kn_filtered_for_cdhit

    script:
    """

    # Run blast -----
    blastn -db nt \
    -query ${kraken_filtered} \
    -out blast_${id}_nt0.txt \
    -outfmt "6 qseqid sseqid stitle pident length \
    evalue qcovs bitscore staxids sscinames sblastnames \
    sskingdoms mismatch gapopen qstart qend qlen sstart send slen" \
    -max_target_seqs 1 \
    -max_hsps 1  \
    -evalue ${params.evalue} \
    -perc_identity ${params.blastnt_pident} \
    -num_threads ${task.cpus}

    # Create headers for the blast file -----
    echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,staxids,\
    sscinames,sblastnames,sskingdoms,mismatch,gapopen,qstart,qend,\
    qlen,sstart,send,slen" | tr ',' '\t'| cat - blast_${id}_nt0.txt > blast_${id}_nt.txt
    rm blast_${id}_nt0.txt

    # Filter the blast NT output to remove contamination ------
    grep -E -v "primates|other sequences" blast_${id}_nt.txt | \
    awk -F"\t" '(\$4 > ${params.blastnt_filter_pident} && \$5 > ${params.blastnt_filter_length})' | cut -f1  > ${id}_IDs_to_remove_nt.txt

    # Filter the fasta from last step (Kraken) by the blast keep list ------
    seqkit grep -i -v -f ${id}_IDs_to_remove_nt.txt ${kraken_filtered} > ${id}_seqkit_kn_blast_filtered.fasta
    
    """
}

/*
  STEP 3: RUN CD-HIT ON FILTERED READS
*/
if (params.skipcdhit == 'false') {
process cdhit {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "CD-HIT/4.8.1-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/CD-HIT/${params.assembler}/",mode:"copy"

    input:
    tuple val(id), file(filtered_blst) from seqkit_blast_kn_filtered_for_cdhit

    output:
    tuple val(id), file('*_seqkit_kn_blast_cdhit.fasta') into seqkit_kn_blast_cdhit_for_blastref,seqkit_kn_blast_cdhit_for_final

    script:
    """
    # Use cd-hit to cluster and remove redundancy ------
    cd-hit-est \
    -i ${filtered_blst} \
    -o ${id}_seqkit_kn_blast_cdhit.fasta \
    -c ${params.cdhit_identity} \
    -n ${params.cdhit_wordsize} \
    -T ${task.cpus}

    """
}

} else {
  
 seqkit_blast_kn_filtered_for_cdhit.into {seqkit_kn_blast_cdhit_for_blastref;seqkit_kn_blast_cdhit_for_final}
}

/*
  STEP 4: RUN BLAST ON FILTERED READS (GRCh38, GRCh38.p0, CHM13) 
*/
process blastref {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "BLAST+/2.10.1-IGB-gcc-8.2.0","seqkit/0.12.1"
    publishDir             "${resultsPath}/BLAST/${params.assembler}/",mode:"copy",overwrite: true

    input:
    tuple val(id), file(seqkit_kn_blast_cdhit_filtered) from seqkit_kn_blast_cdhit_for_blastref
    file genome1 from genome_file1
    file genome2 from genome_file2
    file genome3 from genome_file3
    file "*" from blast_db1_ch
    file "*" from blast_db2_ch
    file "*" from blast_db3_ch

    output:
    tuple val(id), file('*.GRCh38.decoys_blast0.txt') optional true
    tuple val(id), file('*.GRCh38.decoys_blast.txt') 
    tuple val(id), file('*.GRCh38.p0_blast0.txt') optional true
    tuple val(id), file('*.GRCh38.p0_blast.txt') 
    tuple val(id), file('*.CHM13_blast0.txt') optional true
    tuple val(id), file('*.CHM13_blast.txt') 
    tuple val(id), file('*.GRCh38.decoys_IDs_to_filter.txt') into filter_GRCH38_decoy_hla
    tuple val(id), file('*.GRCh38.p0_IDs_to_filter.txt') into filter_GRCH38_p0
    tuple val(id), file('*.CHM13_IDs_to_filter.txt') into filter_CHM13_GRCH38
  
    """
    # Run blast (GRCh38) ------
    blastn -db ${genome1} \
    -query ${seqkit_kn_blast_cdhit_filtered} \
    -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
    -out ${id}.GRCh38.decoys_blast0.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.blastr_pident} \
    -num_threads ${task.cpus}

    # Create headers for the blast file -----
    echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
    qlen,sstart,send,slen" | tr ',' '\t'| cat - ${id}.GRCh38.decoys_blast0.txt > ${id}.GRCh38.decoys_blast.txt
    rm ${id}.GRCh38.decoys_blast0.txt

    # Filter the blast (GRCh38.p0) ------
    awk -F"\t" '(\$4 > ${params.blastr_filter_pident} && \$7 > ${params.blastr_filter_qcov})' \
    ${id}.GRCh38.decoys_blast.txt | cut -f1 | uniq > ${id}.GRCh38.decoys_IDs_to_filter.txt

    # Run blast (GRCh38.p0) ------
    blastn -db ${genome2} \
    -query ${seqkit_kn_blast_cdhit_filtered} \
    -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
    -out ${id}.GRCh38.p0_blast0.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.blastr_pident} \
    -num_threads ${task.cpus}

     # Create headers for the blast file -----
    echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
    qlen,sstart,send,slen" | tr ',' '\t'| cat - ${id}.GRCh38.p0_blast0.txt > ${id}.GRCh38.p0_blast.txt
    rm ${id}.GRCh38.p0_blast0.txt
   
    # Filter the blast (GRCh38.po) ------
      awk -F"\t" '(\$4 > ${params.blastr_filter_pident} && \$7 > ${params.blastr_filter_qcov})' \
      ${id}.GRCh38.p0_blast.txt | cut -f1 | uniq > ${id}.GRCh38.p0_IDs_to_filter.txt

    # Run blast (CHM13 + Chr Y of GRCH38) ------
    blastn -db ${genome3} \
    -query ${seqkit_kn_blast_cdhit_filtered} \
    -outfmt "6 qseqid sseqid stitle pident \
    length evalue qcovs bitscore mismatch \
    gapopen qstart qend qlen sstart send slen" \
    -out ${id}.CHM13_blast0.txt \
    -max_target_seqs ${params.max_target_seqs} \
    -max_hsps ${params.max_hsps}  \
    -evalue ${params.evalue} \
    -perc_identity ${params.blastr_pident} \
    -num_threads ${task.cpus} 

     # Create headers for the blast file -----
    echo -e "qseqid,sseqid,stitle,pident,length,evalue,qcovs,bitscore,mismatch,gapopen,qstart,qend, \
    qlen,sstart,send,slen" | tr ',' '\t'| cat - ${id}.CHM13_blast0.txt > ${id}.CHM13_blast.txt
    rm ${id}.CHM13_blast0.txt

    # Filter the blast (CHM13) ------
    awk -F"\t" '(\$4 > ${params.blastr_filter_pident} && \$7 > ${params.blastr_filter_qcov})' \
    ${id}.CHM13_blast.txt | cut -f1 | uniq > ${id}.CHM13_IDs_to_filter.txt

    """ 
}


/*
  STEP 5: FILTER FASTA FILE BASED ON BLAST STEP
*/
process final_filtering {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "seqkit/0.12.1"
    publishDir             "${resultsPath}/Final-Filtered/${params.assembler}/",mode:"copy",overwrite: true

    input:
    tuple val(id),file(cdhit_filtered),file(filter_GRCH38_decoy),file(filter_GRCH38_p0),file(filter_CHM13) from seqkit_kn_blast_cdhit_for_final.join(filter_GRCH38_decoy_hla).join(filter_GRCH38_p0).join(filter_CHM13_GRCH38)

    output:
    tuple val(id), file('*_GRCH38_decoys_filter.final.fasta')
    tuple val(id), file('*_GRCH38_p0_filter.final.fasta')
    tuple val(id), file('*_CHM13_filter.final.fasta')

    script:
    """
    # Filter the fasta using blast output (GRCh38) ------ [grab what is not in the blast output]
    seqkit grep -i -v -f ${filter_GRCH38_decoy} ${cdhit_filtered} > ${id}_GRCH38_decoys_filter.final.fasta

    # Filter the fasta using blast output (GRCh38.p0) ------
    seqkit grep -i -v -f ${filter_GRCH38_p0} ${cdhit_filtered} > ${id}_GRCH38_p0_filter.final.fasta

    # Filter the fasta using blast output (CHM13) ------
    seqkit grep -i -v -f ${filter_CHM13} ${cdhit_filtered} > ${id}_CHM13_filter.final.fasta
    """
}

