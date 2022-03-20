#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = "./results"    /*output folder, must specify complete path. Required parameters*/
params.assembler             = 'masurca'      /*options: megahit|masurca. Default is megahit*/
params.skipKraken2           = true           /*options: true|false. Default is true which means that kraken2 will be skipped*/

/* parameters for readprep = qctrimming and adapter removal */
params.skipTrim              = false           /*qc-trimming of reads. options: true|false. Default is false*/     
params.min_read_length       = '20'            /*minimum length of read to be kept after trimming for downstream analysis. Default is 20*/
params.min_base_quality      = '10'            /*minimum base quality. Default is 20*/
params.guess_adapter         = true            /*auto-detect adapter from input file. options: true|false. Default is true*/

/* Experimental */
params.tmpdir               = '/scratch'       /*primarily for setting up a defined tmp/scratch space for samtools collate*/

/* Not yet implemented */
params.forward_adapter       = false           /*adapter sequence to be clipped off (forward). */
params.reverse_adapter       = false           /*adapter sequence to be clipped off (reverse). Used for paired reads only*.*/

// params.singleEnd             = false          /*options: true|false. true = the input type is single end reads; false = the input type is paired reads. Default is false*/

/*Stage*/
stage = "assembly"

supportedAsm = ["masurca", "megahit"]

asms = params.assembler.split(',')

for (assembler : asms)
    if (!supportedAsm.contains(assembler)) exit 1, "Unknown assembler, ${assembler}"

resultsPath = "${params.outputDir}/${stage}"

/*cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbio'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '100'
params.clusterAcct           = " -A h3bionet "

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()

// Sanity checks
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"

Channel.fromFilePairs("${params.samplePath}", size: 1)
    .into { CRAM_Ch1; CRAM_Ch2 }

/*
  prepare_genome 
  This process is executed only once
*/

process prepare_genome {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0

    container              "https://depot.galaxyproject.org/singularity/mulled-v2-7ef549f04aa19ef9cd7d243acfee913d928d9b88:f5ff855ea25c94266e524d08d6668ce6c7824604-0"
    tag                    { "PREP:${genome}" }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "SAMtools/1.12-IGB-gcc-8.2.0"
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

/*
  qc_input 
  This process mainly checks the inputs for improperly formed CRAM/BAM input files
*/

process qc_input {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    // 
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    errorStrategy          { task.exitStatus=1 ? 'ignore' : 'terminate' }
    
    input:
    tuple val(id), file(CRAM) from CRAM_Ch1

    output:
    tuple val(id), file("${id}.final.cram") optional true into extract_unmapped_ch
    
    script:
    """
    samtools quickcheck ${CRAM}
    """
}    

/*
   Read extraction.  This step is tricky.
   
   For single end reads, we only need to extract unmapped as there isn't a mapped mate.  
   However this also means that the additional workflows downstream cannot be run as they 
   paired end read data (location placement requires having mapped mates)
   
   samtools view -f 4 > se_unmapped.bam
   For paired end reads, we need to extract three classes of reads, keeping the mates
   
   1. Both unmapped
   2. R1 mapped, R2 unmapped
   3. R2 mapped, R1 unmapped
   
   This can be accomplished with samtools though the logic is a little tricky with the bit
   flags.  Best way would be to extract pairs where any of the reads are unmapped, then 
   pull out subtuples using flags. Otherwise we're running through the 
   
   # both unmapped; -f means unmapped, 
   # bit flag 12  = both reads unmapped (bit flag 4 & 8)
   # bit flag 2304 = not primary alignment (this removes these), not supplemental
   samtools view -hb -f 12 -F 2304 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = R1 unmapped
   # bit flag 2312 = Mate unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                tuple
   samtools view -hb -f 4 -F 2312 alignments.bam  > R1_unmapped.bam
   
   # R2 mapped, R1 not
   # bit flag 8    = R2 (mate) unmapped
   # bit flag 2308 = R1 (read) unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                tuple
   samtools view -hb -f 8 -F 2308 alignments.bam  > R2_unmapped.bam
*/

/*
  extract_unmapped
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step
*/

process extract_improper {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { id }
    executor               myExecutor
    cpus                   6
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Read-Prep/Improper",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(cram) from extract_unmapped_ch
    file genome from genome_file
    file index from genome_index_ch

    output:
    // all unaligned + mates
    tuple val(id), file("${id}.improper.bam") into improper_unmapped_ch, improper_clipped_ch
    
    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority and the logic for extracting discordant
    // reads is an issue w/ SE data
    
    script:
    """
    # grab any non-properly paired reads (includes any unmapped and discordant reads) 
    samtools view -@ ${task.cpus} -hbt ${index} -G 2 -o ${id}.improper.bam ${cram}
    """
}

process extract_unmapped {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { id }
    executor               myExecutor
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Read-Prep/Unmapped",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(bam) from improper_unmapped_ch

    output:
    tuple val(id), file("${id}.all-unmapped.R{1,2}.fastq.gz") optional true into fq_pe_unmapped_ch
    tuple val(id), file("${id}.unmapped.bam") optional true into unmapped_bam_ch  // all unaligned + mates

    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority yet.
    
    script:    
    """    
    # both unmapped
    samtools view -@ ${task.cpus} -hb \\
        -f 12 -F 2304 -o ${id}.both-unmapped.bam ${bam}

    samtools fastq -@ ${task.cpus} ${id}.both-unmapped.bam \\
        -1 ${id}.both-unmapped.R1.fastq.gz -2 ${id}.both-unmapped.R2.fastq.gz

    # R1 only unmapped
    samtools view -@ ${task.cpus} -hb \\
        -f 4 -F 2312 -o ${id}.R1-unmapped.bam ${bam}
    samtools fastq -@ ${task.cpus} ${id}.R1-unmapped.bam \\
        -1 ${id}.R1-unmapped.R1.fastq.gz -2 ${id}.R1-unmapped.R2.fastq.gz
    
    # R2 only unmapped
    samtools view -@ ${task.cpus} -hb \\
        -f 8 -F 2308 -o ${id}.R2-unmapped.bam ${bam}    
    samtools fastq -@ ${task.cpus} ${id}.R2-unmapped.bam \\
        -1 ${id}.R2-unmapped.R1.fastq.gz -2 ${id}.R2-unmapped.R2.fastq.gz
    
    # combine unmapped BAM files for later analysis
    samtools merge -@ ${task.cpus} \
        ${id}.unmapped.bam \
        ${id}.both-unmapped.bam \
        ${id}.R1-unmapped.bam \
        ${id}.R2-unmapped.bam
    
    ## NOTE: in the below code we are keeping *all* paired reads if either or both are unmapped. 
    # combine R1 reads
    cat ${id}.both-unmapped.R1.fastq.gz ${id}.R1-unmapped.R1.fastq.gz ${id}.R2-unmapped.R1.fastq.gz > ${id}.all-unmapped.R1.fastq.gz

    # combine R2 reads
    cat ${id}.both-unmapped.R2.fastq.gz ${id}.R1-unmapped.R2.fastq.gz ${id}.R2-unmapped.R2.fastq.gz > ${id}.all-unmapped.R2.fastq.gz

    # combined SE unmapped FASTQ into one file for assembly
    ## cat ${id}.R1-unmapped.R1.fastq ${id}.R2-unmapped.R2.fastq >> ${id}.orphans.unmapped.fastq
    
    # we could combine the mapped reads here, but we'll use the original BAM instead
    ## cat ${id}.R1-unmapped.R2.fastq ${id}.R2-unmapped.R1.fastq >> ${id}.orphans.mapped.fastq
    """
}

// TODO: note the used of the hard-coded '/scratch' space here.  
// We could default to the current work directory, but the current behavior of 
// samtools collate is to use /tmp and does *not* currently use TMPDIR (yeah, pretty crappy)

process extract_clipped {
    // this may need a mulled repo with python and samtools:
    // https://github.com/BioContainers/multi-package-containers
    // best would be a single container w/ bwa, samtools 1.12+, and python 3.7+
    container              null
    tag                    { id }
    executor               myExecutor
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0","Python/3.7.2-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Read-Prep/Clipped",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(bam) from improper_clipped_ch

    output:
    tuple val(id), file("${id}.clipped.bam") into extract_sort_ch // all clipped pairs; we save this for now
    tuple val(id), file("${id}.all-clipped.R{1,2}.fastq.gz") into fq_pe_clipped_ch

    script:
    """
    # capture only discordant clipped reads; note the -G 2
    # the below is to prevent collisions if the pipeline is interrupted and restarted
    tmpdir=\$( mktemp -d collate.XXXXXXXXX )
    samtools collate --output-fmt bam -@ ${task.cpus - 4} -O ${bam} \$tmpdir/${id} | \
        clipped-filter.py > ${id}.clipped.tmp.bam

    samtools fastq -@ ${task.cpus} ${id}.clipped.tmp.bam \
        -1 ${id}.all-clipped.R1.fastq.gz -2 ${id}.all-clipped.R2.fastq.gz
    
    samtools sort -@ ${task.cpus} -o ${id}.clipped.bam ${id}.clipped.tmp.bam
    samtools index ${id}.clipped.bam*
    """
}

process merge_pairs {
    // singularity run https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0
    container              "https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0"
    tag                    { id }
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "seqkit/0.12.1"
    publishDir             "${resultsPath}/Read-Prep/Merged",mode:"copy", overwrite: true
 
    input:
    tuple val(id), file(unmapped), file(clipped) from fq_pe_unmapped_ch.join(fq_pe_clipped_ch)

    output:
    tuple val(id), file("${id}.all-reads.R{1,2}.fastq.gz") into merge_trim_ch, merge_fastqc_ch
    file("*.rmdup.txt")

    """
    cat ${unmapped[0]} ${clipped[0]} > ${id}.all-reads.R1.tmp.fastq.gz
    cat ${unmapped[1]} ${clipped[1]} > ${id}.all-reads.R2.tmp.fastq.gz

    seqkit rmdup -n -j ${task.cpus} -o ${id}.all-reads.R1.fastq.gz ${id}.all-reads.R1.tmp.fastq.gz 2> ${id}.R1.rmdup.txt
    seqkit rmdup -n -j ${task.cpus} -o ${id}.all-reads.R2.fastq.gz ${id}.all-reads.R2.tmp.fastq.gz 2> ${id}.R2.rmdup.txt

    rm ${id}.all-reads.R1.tmp.fastq.gz ${id}.all-reads.R2.tmp.fastq.gz
    """
} 

process fastqc {
    // singularity run https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1
    container              "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1"
    tag "FASTQC-Pretrim ${id}"
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "12 GB"
    // module                 "FastQC/0.11.8-Java-1.8.0_152"
    publishDir             "${resultsPath}/FASTQC-Pretrim"

    input:
    tuple val(id), file(reads) from merge_fastqc_ch

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

/*
  trimming
*/

process trimming {
    // singularity run https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0
    container              "https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0"
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Read-Prep/Trimmed",mode:"copy", overwrite: true
    // module                 "fastp/0.20.0-IGB-gcc-4.9.4"

    input:
    tuple val(name), file(reads) from merge_trim_ch
    
    output:
    tuple val(name), file('*.PE.R{1,2}.trimmed.fastq.gz'), file('*.unpR{1,2}.trimmed.fastq.gz') optional true into trim_megahit_ch, trim_masurca_ch, trim_fastqc, trim_aln_ch
    tuple val(name), file('*.json') optional true into trim_multiqc_ch
    file '*.html'
        
    script:
    trimOptions      = params.skipTrim ? " " :  " -l ${params.min_read_length} --cut_right --cut_right_window_size 3 --cut_right_mean_quality ${params.min_base_quality}"
    adapterOptionsSE = params.guess_adapter ? ' ' : " --adapter_sequence=${params.forward_adapter} "
    adapterOptionsPE = params.guess_adapter ? ' --detect_adapter_for_pe ' : " --adapter_sequence=${params.forward_adapter}  --adapter_sequence_r2=${params.reverse_adapter} "

    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority yet.

    """
    fastp --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 "${name}.PE.R1.trimmed.fastq.gz"  \
        --out2 "${name}.PE.R2.trimmed.fastq.gz" \
        --unpaired1 "${name}.unpR1.trimmed.fastq.gz"\
        --unpaired2 "${name}.unpR2.trimmed.fastq.gz" \
        ${adapterOptionsPE}  ${trimOptions} \
        --thread ${task.cpus} -w ${task.cpus} \
        --html "${name}"_PE_fastp.html \
        --json "${name}"_PE_fastp.json
    """
}

process fastqc_post {
    // https://biocontainers.pro/tools/fastqc
    // singularity run https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1
    tag "FASTQC-Posttrim ${id}"
    container              "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1"
    executor               myExecutor
    clusterOptions         params.clusterAcct     
    cpus                   4
    queue                  params.myQueue
    memory                 "12 GB"
    // module                 "FastQC/0.11.8-Java-1.8.0_152"
    publishDir             "${resultsPath}/FASTQC-Posttrim",mode:"copy", overwrite: true

    input:
    tuple val(id), file(pereads), file(sereads) from trim_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_trimmed_results

    script:
    """
    fastqc --quiet --threads $task.cpus $pereads $sereads
    """
}

/*
  *Megahit for different input data types:
  *megahit -1 pe_1.fq -2 pe_2.fq -o out  # 1 paired-end library
  *megahit --12 interleaved.fq -o out # one paired & interleaved paired-end library
  *megahit -1 a1.fq,b1.fq,c1.fq -2 a2.fq,b2.fq,c2.fq -r se1.fq,se2.fq -o out # 3 paired-end libraries + 2 SE libraries
*/

process megahit_assemble {
    // https://biocontainers.pro/tools/megahit
    tag                    { name }
    container              "https://depot.galaxyproject.org/singularity/megahit:1.2.9--h2e03b76_1"
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   12
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    // module                 "MEGAHIT/1.2.9-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Raw-Assembly/megahit",mode:"copy",overwrite: true
    
    when:
    asms.contains('megahit')

    input:
    tuple val(name), file(pefastqs), file(sefastqs) from trim_megahit_ch

    output:
    tuple val(name), val("megahit"), file("${name}/final.contigs.fa") into megahit_rename_ch
    file("${name}/*")

    script:
    """
    megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]} \
        -r ${sefastqs[0]},${sefastqs[1]} \
        -o ${name}
    # megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]} -o ${name}.megahit_results
    """
}

/*
  masurca
*/

process masurca_assemble {
    container              "docker://quay.io/h3abionet_org/masurca"
    stageInMode            "copy"
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct
    cpus                   12
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 "MaSuRCA/3.4.2-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Raw-Assembly/masurca", mode:"copy", overwrite: true

    when:
    asms.contains('masurca')

    input:
    tuple val(name), file("${name}/pefastq*.trimmed.fastq.gz"), file("${name}/sefastq*.trimmed.fastq.gz") from trim_masurca_ch

    output:
    tuple val(name), val("masurca"), file("${name}/CA/final.genome.scf.fasta") into masurca_rename_ch
    file("${name}/*")

    script:
    """
    cd ${name}

    gunzip *.trimmed.fastq.gz

    cat << EOF > ${name}.masurca_config_file.txt
    DATA
    PE = pe 300 50 pefastq1.trimmed.fastq pefastq2.trimmed.fastq
    PE = s1 300 50 sefastq1.trimmed.fastq
    PE = s2 300 50 sefastq2.trimmed.fastq
    END

    PARAMETERS
    GRAPH_KMER_SIZE = auto
    USE_LINKING_MATES = 1
    LIMIT_JUMP_COVERAGE = 300
    CA_PARAMETERS = cgwErrorRate=0.15
    KMER_COUNT_THRESHOLD = 1
    NUM_THREADS = ${task.cpus}
    JF_SIZE = 200000000
    CLOSE_GAPS=0
    SOAP_ASSEMBLY=0
    DO_HOMOPOLYMER_TRIM=0
    END

    EOF

    masurca ${name}.masurca_config_file.txt

    ./assemble.sh
    """
}

all_assemblies_rename_ch = megahit_rename_ch.mix(masurca_rename_ch)

process assembly_rename {
    // singularity run https://depot.galaxyproject.org/singularity/perl:5.26.2
    container              "https://depot.galaxyproject.org/singularity/perl:5.26.2"
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   1
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    // TODO: a base perl install is fine, but we have this as a placeholder
    // just in case to remind us for Docker/Singularity
    // module                 "Perl/5.24.1-IGB-gcc-4.9.4"
    publishDir             "${resultsPath}/Final-Assembly/${assembler}",mode:"copy", overwrite: true

    input:
    tuple val(name), val(assembler), file(assembly) from all_assemblies_rename_ch

    output:
    tuple val(name), val(assembler), file("${name}.${assembler}.final.fasta") into all_assemblies_metrics_ch,all_assemblies_aln_ch

    script:
    // Strip off anything after the first '.'
    shortname = name.replaceAll(~/\.\S+$/, "")
    """
    # this could be replaced with something else if needed for a simpler container
    perl -p -e 's/^>(\\N+)/>${name}:\$1/' ${assembly} > ${name}.${assembler}.final.fasta
    """
}

// all_assemblies_metrics_ch = megahit_metrics_ch.mix(masurca_metrics_ch)

process assembly_metrics {
    // singularity run https://depot.galaxyproject.org/singularity/quast:5.0.1--py36pl526ha92aebf_0

    // there is apparently an issue with the locale settings for quast in biocontainers:
    // https://github.com/BioContainers/containers/issues/206
    // this doesn't seem to be resolved yet.

    // testing other containers in the meantime

    tag {id}
    container              "staphb/quast:5.0.2"
    executor               myExecutor
    clusterOptions         params.clusterAcct     
    cpus                   4
    queue                  params.myQueue
    memory                 "12 GB"
    // module                 "quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             "${resultsPath}/QUAST/", mode:"copy", overwrite: true
    // errorStrategy          { task.exitStatus=4 ? 'ignore' : 'terminate' }

    input:
    tuple val(id), val(assembler), file(asm) from all_assemblies_metrics_ch

    output:
    file "${id}.${assembler}/*" 
    file "${id}.${assembler}" into metrics_multiqc_ch

    script:
    """
    quast.py \\
        --output-dir ${id}.${assembler} \\
        --threads ${task.cpus} \\
        ${asm}
    """
}

process aln_reads {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    tag {id}
    container              "https://depot.galaxyproject.org/singularity/mulled-v2-7ef549f04aa19ef9cd7d243acfee913d928d9b88:f5ff855ea25c94266e524d08d6668ce6c7824604-0"
    executor               myExecutor
    cpus                   8
    queue                  params.myQueue
    clusterOptions         params.clusterAcct     
    memory                 "24 GB"
    // module                 "BWA/0.7.17-IGB-gcc-8.2.0","SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/BWA-MEM/${assembler}",mode:"copy",overwrite: true
    // errorStrategy          { task.attempt == 5 ? 'retry' : 'ignore' }

    input:
    tuple val(id), val(assembler), file(assembly), file(pereads), file(sereads) from all_assemblies_aln_ch.combine(trim_aln_ch, by:0)

    output:
    file "${id}.${assembler}.sorted.pe.bam*"
    file "${id}.${assembler}.sorted.se.bam*"
    file "${id}.${assembler}.sorted.merged.bam*"

    script:
    """
    bwa index ${assembly}

    # PE
    bwa mem -t ${task.cpus - 4} ${assembly} ${pereads[0]} ${pereads[1]} | \
        samtools sort -@ 4 -O bam -T ${id} -o ${id}.${assembler}.sorted.pe.bam
    samtools index ${id}.${assembler}.sorted.pe.bam

    # SE
    cat ${sereads[0]} ${sereads[1]} > ${id}.se.fastq.gz
    bwa mem -t ${task.cpus - 4} ${assembly} ${id}.se.fastq.gz | \
        samtools sort -@ 4 -O bam -T ${id} -o ${id}.${assembler}.sorted.se.bam
    samtools index ${id}.${assembler}.sorted.se.bam

    # merge both files
    samtools merge -@ ${task.cpus} ${id}.${assembler}.sorted.merged.bam ${id}.${assembler}.sorted.pe.bam ${id}.${assembler}.sorted.se.bam
    samtools index ${id}.${assembler}.sorted.merged.bam 
    """
}

process MultiQC {
    // singularity run https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0
    container              "https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0"
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    // module                 "MultiQC/1.11-IGB-gcc-8.2.0-Python-3.7.2"
    publishDir             "${resultsPath}/MultiQC",mode:"copy",overwrite: true
 
    input:
    file('./QUAST/*') from metrics_multiqc_ch.collect().ifEmpty([])
    file('./FASTQC-Pretrim/*') from fastqc_results.collect().ifEmpty([])
    file('./FASTQC-Posttrim/*') from fastqc_trimmed_results.collect().ifEmpty([])
    file('./FASTP/*') from trim_multiqc_ch.collect().ifEmpty([])

    output:
    file "multiqc*"

    """
    multiqc  .
    """
} 
