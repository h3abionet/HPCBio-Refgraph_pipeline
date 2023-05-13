#!/usr/bin/env nextflow

nextflow.enable.dsl=1
/*parameters that are specified at the command line or via config file*/
params.reference             = false          /* reference genome fasta file, must specify complete path. Required parameter if CRAM is used */
params.samplePath            = false          /* input folder, must specify complete path. Required parameters */
params.outputDir             = "results"      /* output folder, must specify complete path. Required parameters */
params.format                = "cram"         /* input format (default: 'cram', options = 'cram', 'bam') */
params.genome                = false

/* NYI, original alignment filtering criteria (for proper pair distance) */
// params.meanInsSize           = 300
// params.stdevInsSize          = 100

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

/*Stage*/
stage = "read-extraction"

resultsPath = "${params.outputDir}/${stage}"

/* cluster parameters */
myExecutor                   = 'slurm'
params.myQueue               = 'hpcbioamd,hpcbio'
defaultCPU                   = '1'
defaultMemory                = '20'

/* Prepare input if given */
reference_file               = params.reference ? file(params.reference) : Channel.empty()
reference_store              = params.reference ? reference_file.getParent() : ''

/* Alternative genome (only used for additional alignment) */
genome_file                  = params.genome ? file(params.genome) : Channel.empty()
genome_store                 = params.genome ? genome_file.getParent() : ''

Channel.fromFilePairs("${params.samplePath}", size: 1)
    .into { Aln_QC; Aln_Stats; Aln_Improper; Aln_Clipped }

// TODO: add general sanity checks on required data and formats

// Check format to make sure it's 'bam' or 'cram'
// Check that a reference is provided if the format is 'cram'

if (params.reference) {
    // needed for CRAM, see sanity checks
    if( !reference_file.exists() ) exit 1, "Missing reference genome file: ${reference_file}"

    process prepare_reference {
        container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
        tag                    { "PREP_REFERENCE" }
        cpus                   defaultCPU
        memory                 "$defaultMemory GB"
        storeDir               reference_store
        
        input:
        file reference from reference_file

        output:
        file "*.fai" into reference_index_ch
        
        script:
        """
        samtools faidx ${reference}
        """
    }
} else {
    // BAM only
    reference_index_ch = Channel.empty()
}

if (params.genome) {
    if( !genome_file.exists() ) exit 1, "Missing alt genome file: ${genome_file}"

    process SAMTOOLS_FAIDX {
        container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
        tag                    { "PREP_ALT_GENOME:SAMTOOLS" }
        cpus                   defaultCPU
        memory                 "$defaultMemory GB"
        storeDir               genome_store
        
        input:
        file genome from genome_file

        output:
        file "*.fai" into genome_fai_ch
        
        script:
        """
        samtools faidx ${genome}
        """
    }

    process BWA_INDEX {
        container              "https://depot.galaxyproject.org/singularity/bwa:0.7.17--h7132678_9"
        tag                    { "PREP_ALT_GENOME:BWA" }
        cpus                   2
        memory                 "24GB"
        storeDir               genome_store
        
        input:
        file genome from genome_file

        output:
        file "*.{ann,amb,pac,bwt,sa}" into bwa_index_ch_pe
        
        script:
        """
        bwa index ${genome}
        """
    }
} else {
    // BAM only
    Channel.empty().
        into {  bwa_index_ch_pe;genome_fai_ch }
}

/*
  qc_input 
  This process mainly checks the inputs for improperly formed CRAM/BAM input files
*/

process CHECK_INPUT {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    // 
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { "CHECK_INPUT:${id}" }
    cpus                   defaultCPU
    memory                 "$defaultMemory GB"
    errorStrategy          { task.exitStatus=1 ? 'ignore' : 'terminate' }
    
    input:
    tuple val(id), file(aln) from Aln_QC

    output:
    // tuple file("${id}.stats"), file("${id}.idxstat"), file("${id}.flagstat")
    tuple val(id), val("${task.exitStatus}") into Aln_QC_Check_Stats, Aln_QC_Check_Extract, Aln_QC_Check_Clipped

    script:
    ref = params.reference ? "--reference ${reference_file}" : ""

    """
    samtools quickcheck ${ref} ${aln}
    """
}

process ALIGNMENT_STATS {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { "ALN_STATS:${id}" }
    cpus                   defaultCPU
    memory                 "$defaultMemory GB"
    errorStrategy          { task.exitStatus=1 ? 'ignore' : 'terminate' }
    
    input:
    tuple val(id), file(aln), val(status) from Aln_Stats.join(Aln_QC_Check_Stats)

    output:
    file("${id}.original.stats")
    file("${id}.original,idxstat")
    file("${id}.original.flagstat")

    script:
    ref = params.genome ? "--reference ${genome_file}" : ""

    """
    samtools stats -@ ${task.cpus} ${aln} > ${id}.original.stats
    # TODO: check for the presence of an index file
    # samtools idxstat -@ ${task.cpus} ${aln} > ${id}.original.idxstat
    samtools flagstat -@ ${task.cpus} ${aln} > ${id}.original.flagstat
    """
}

/*
   Read extraction.  This step is tricky.
   
   For single end reads, we only need to extract unmapped as there isn't a mapped mate.  
   However this also means that the additional workflows downstream cannot be run as they 
   paired end read data (location placement requires having mapped mates)
   
   samtools view -f 4 > se_unmapped.bam

   For paired end reads, we have two options now.  

   Originally we would need to extract three classes of reads, keeping the mates for each:
   
   1. Both unmapped
   2. R1 mapped, R2 unmapped
   3. R2 mapped, R1 unmapped
   
   This can be accomplished with samtools though the logic is a little tricky with the bit
   flags.  Best way would be to extract pairs where any of the reads are unmapped, then 
   pull out subtuples using flags.
   
   # You can get all improper pairs, which seem to include any unmapped reads, by using -G 2 

   # both unmapped; -f means unmapped, 
   # bit flag 12  = both reads unmapped (bit flag 4 & 8)
   # bit flag 2304 = not primary alignment (this removes these), not supplemental
   samtools view -hb -f 12 -F 2304 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = first in pair (R1), read unmapped (this could be bit flag 68)
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

   This is largely because the original versions of samtools didn't have an option 
   for retaining or filtering reads that matched *any* one FLAG (for example, if either the 
   read is unmapped or it's mate is unmapped).  More recent versios of samtools now have the `--rf` 
   option, so you can do this in one step to capture any unmapped reads:

   samtools view -hb --rf 12 -F 2304 alignments.bam  > R2_unmapped.bam
*/

// Aln_Ch1.join(Aln_QC_Check).dump()

process EXTRACT_IMPROPER {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { "EXTRACT_BAM:${id}" }
    cpus                   12
    publishDir             "${resultsPath}/Read-Prep/Improper",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(aln), val(status) from Aln_Improper.join(Aln_QC_Check_Extract)
    // path(genome) from reference_file
    // path(index) from reference_index_ch

    output:
    // all unaligned + mates
    tuple val(id), file("${id}.improper.bam") into improper_unmapped_ch, improper_clipped_ch
    
    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority and the logic for extracting discordant
    // reads is an issue w/ SE data
    
    script:
    ref = params.format == 'cram' ? "-t ${genome}" : ''
    """
    # grab any non-properly paired reads (includes any unmapped and discordant reads) 
    samtools view -@ ${task.cpus} \
        -hb \
        ${ref} \
        -G 2 \
        -o ${id}.improper.bam \
        ${aln}

    # samtools view -@ ${task.cpus} \
    #    -uhb \
    #    ${ref} \
    #    -f 4 -F2312 \
    #    -o ${id}.R1-unmapped.tmp.bam \
    #    ${aln}

    # samtools view -@ ${task.cpus} \
    #    -uhb \
    #    ${ref} \
    #    -f 8 -F2308 \
    #    -o ${id}.R2-unmapped.tmp.bam \
    #    ${aln}

    # samtools view -@ ${task.cpus} \
    #    -uhb \
    #    ${ref} \
    #    -f 12 -F2304 \
    #    -o ${id}.both-unmapped.tmp.bam \
    #    ${aln}

    # samtools merge -@ ${task.cpus-6} \
    #    -u - \
    #    ${id}.R1-unmapped.tmp.bam \
    #    ${id}.R2-unmapped.tmp.bam \
    #    ${id}.both-unmapped.tmp.bam |\
    #    samtools sort -@ 6 \
    #    -o ${id}.unmapped.bam 
    """
}

process EXTRACT_UNMAPPED {
    // singularity run https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0
    container              "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    tag                    { "EXTRACT_UNMAPPED:${id}" }
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Read-Prep/Unmapped",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(bam) from improper_unmapped_ch

    output:
    tuple val(id), file("${id}.any-unmapped.R{1,2}.fastq.gz") optional true into fq_pe_unmapped_ch
    file("${id}.any-unmapped.singletons.fastq.gz") optional true
    tuple val(id), file("${id}.unmapped.bam") optional true into unmapped_bam_ch  // all unaligned + mates
    
    script:    
    """    
    # any unmapped
    samtools view -@ ${task.cpus} -hb \\
        --rf 12 -F 2304 -o ${id}.any-unmapped.bam ${bam}

    # note the singletons flag, this captures 
    samtools fastq -@ ${task.cpus} \\
        -1 ${id}.any-unmapped.R1.fastq.gz \\
        -2 ${id}.any-unmapped.R2.fastq.gz \\
        -s ${id}.any-unmapped.singletons.fastq.gz \\
        ${id}.any-unmapped.bam

    ## R1 only unmapped
    #samtools view -@ ${task.cpus} -hb \\
    #    -f 4 -F 2312 -o ${id}.R1-unmapped.bam ${bam}

    #samtools fastq -@ ${task.cpus} \\
    #    -1 ${id}.R1-unmapped.R1.fastq.gz \\
    #    -2 ${id}.R1-unmapped.R2.fastq.gz \\
    #    -s ${id}.R1-unmapped.singletons.fastq.gz \\
    #    ${id}.R1-unmapped.bam
    
    ## R2 only unmapped
    #samtools view -@ ${task.cpus} -hb \\
    #    -f 8 -F 2308 -o ${id}.R2-unmapped.bam ${bam}

    #samtools fastq -@ ${task.cpus} \\
    #    -1 ${id}.R2-unmapped.R1.fastq.gz \\
    #    -2 ${id}.R2-unmapped.R2.fastq.gz \\
    #    -s ${id}.R2-unmapped.singletons.fastq.gz \\
    #    ${id}.R2-unmapped.bam 
    
    # combine unmapped BAM files for later analysis
    #samtools merge -@ ${task.cpus} \\
    #    ${id}.both-unmapped.bam \\
    #    ${id}.R1-unmapped.bam \\
    #    ${id}.R2-unmapped.bam
    
    ## NOTE: in the below code we are keeping *all* paired reads if either or both are unmapped. 
    # combine R1 reads
    #cat ${id}.both-unmapped.R1.fastq.gz \\
    #    ${id}.R1-unmapped.R1.fastq.gz \\
    #    ${id}.R2-unmapped.R1.fastq.gz \\
    #    > ${id}.all-unmapped.R1.fastq.gz

    # combine R2 reads
    #cat ${id}.both-unmapped.R2.fastq.gz \\
    #    ${id}.R1-unmapped.R2.fastq.gz \\
    #    ${id}.R2-unmapped.R2.fastq.gz \\
    #    > ${id}.all-unmapped.R2.fastq.gz

    # cat singletons
    #cat ${id}.both-unmapped.singletons.fastq.gz \\
    #    ${id}.R1-unmapped.singletons.fastq.gz \\
    #    ${id}.R2-unmapped.singletons.fastq.gz \\
    #    > ${id}.all-unmapped.singletons.fastq.gz
    """
}

process EXTRACT_CLIPPED {
    // this may need a mulled repo with python and samtools:
    // https://github.com/BioContainers/multi-package-containers
    // best would be a single container w/ bwa, samtools 1.12+, and python 3.7+
    container              null
    tag                    { "EXTRACT_CLIPPED:${id}" }
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.17-IGB-gcc-8.2.0","Python/3.7.2-IGB-gcc-8.2.0"
    publishDir             "${resultsPath}/Read-Prep/Clipped",mode:"copy", overwrite: true
    
    input:
    tuple val(id), file(bam), val(status) from Aln_Clipped.join(Aln_QC_Check_Clipped)

    output:
    tuple val(id), file("${id}.clipped.bam") into extract_sort_ch // all clipped pairs; we save this for now
    tuple val(id), file("${id}.all-clipped.R{1,2}.fastq.gz") into fq_pe_clipped_ch

    script:
    """
    # the below is to prevent collisions if the pipeline is interrupted and restarted
    tmpdir=\$( mktemp -d collate.XXXXXXXXX )
    samtools collate --output-fmt bam \\
        -@ ${task.cpus - 4} \\
        -O ${bam} \\
        \$tmpdir/${id} | \\
        clipped-filter.py > ${id}.clipped.tmp.bam

    samtools fastq -@ ${task.cpus} \\
        -1 ${id}.all-clipped.R1.fastq.gz \\
        -2 ${id}.all-clipped.R2.fastq.gz \\
        ${id}.clipped.tmp.bam
    
    samtools sort -@ ${task.cpus} \\
        -o ${id}.clipped.bam \\
        ${id}.clipped.tmp.bam
    samtools index ${id}.clipped.bam
    """
}

process MERGE_PAIRS {
    // singularity run https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0
    container              "https://depot.galaxyproject.org/singularity/seqkit:2.1.0--h9ee0642_0"
    tag                    { "MERGE_READS:${id}" }
    cpus                   2
    memory                 "$defaultMemory GB"
    publishDir             "${resultsPath}/Read-Prep/Merged",mode:"copy", overwrite: true
 
    input:
    tuple val(id), file(unmapped), file(clipped) from fq_pe_unmapped_ch.join(fq_pe_clipped_ch)

    output:
    tuple val(id), file("${id}.all-reads.R{1,2}.fastq.gz") into merge_trim_ch, merge_fastqc_ch
    file("*.rmdup.txt")

    """
    cat ${unmapped[0]} ${clipped[0]} > ${id}.all-reads.R1.tmp.fastq.gz
    cat ${unmapped[1]} ${clipped[1]} > ${id}.all-reads.R2.tmp.fastq.gz

    seqkit rmdup -n \
        -j ${task.cpus} ${id}.all-reads.R1.tmp.fastq.gz 2> ${id}.R1.rmdup.txt \
        | seqkit sort -n -o ${id}.all-reads.R1.fastq.gz 
    
    seqkit rmdup -n \
        -j ${task.cpus} ${id}.all-reads.R2.tmp.fastq.gz 2> ${id}.R2.rmdup.txt \
        | seqkit sort -n -o ${id}.all-reads.R2.fastq.gz

    rm ${id}.all-reads.R1.tmp.fastq.gz ${id}.all-reads.R2.tmp.fastq.gz
    """
} 

process FASTQC {
    // singularity run https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1
    container              "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1"
    tag                    {"FASTQC-PRETRIM:${id}"}
    cpus                   2
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

if (!params.skipTrim) {

    process FASTP_TRIM {
        // singularity run https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0
        container              "https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0"
        tag                    "FASTP_TRIM:${id}"
        cpus                   2
        memory                 "$defaultMemory GB"
        publishDir             "${resultsPath}/Read-Prep/Trimmed",mode:"copy", overwrite: true
        // module                 "fastp/0.20.0-IGB-gcc-4.9.4"

        input:
        tuple val(id), file(reads) from merge_trim_ch
        
        output:
        tuple val(id), file('*.PE.R{1,2}.trimmed.fastq.gz'), file('*.unpR{1,2}.trimmed.fastq.gz') optional true into trim_fastqc, trim_aln_ch
        tuple val(id), file('*.json') optional true into trim_multiqc_ch
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
            --unpaired1 "${name}.unpR1.trimmed.fastq.gz" \
            --unpaired2 "${name}.unpR2.trimmed.fastq.gz" \
            ${adapterOptionsPE}  ${trimOptions} \
            --thread ${task.cpus} -w ${task.cpus} \
            --html "${name}"_PE_fastp.html \
            --json "${name}"_PE_fastp.json
        """
    }
} else {
    merge_trim_ch.into{ trim_fastqc;trim_aln_ch }
    trim_multiqc_ch = Channel.empty()
}

process FASTQC_POST {
    // https://biocontainers.pro/tools/fastqc
    // singularity run https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1
    tag                    {"FASTQC-POSTTRIM:${id}"}
    container              "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1"
    cpus                   4
    memory                 "12 GB"
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

process MULTIQC {
    // singularity run https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0
    container              "https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0"
    cpus                   2
    memory                 "$defaultMemory GB"
    // module                 "MultiQC/1.11-IGB-gcc-8.2.0-Python-3.7.2"
    publishDir             "${resultsPath}/MultiQC",mode:"copy",overwrite: true
 
    input:
    file('./FASTQC-Pretrim/*') from fastqc_results.collect().ifEmpty([])
    file('./FASTQC-Posttrim/*') from fastqc_trimmed_results.collect().ifEmpty([])
    file('./FASTP/*') from trim_multiqc_ch.collect().ifEmpty([])

    output:
    file "multiqc*"

    """
    multiqc .
    """
} 

if (params.genome) {

    process BWA_MEM {
        tag                    { "BWA_MEM:${id}" }
        container              "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0"
        cpus                   12
        memory                 "24 GB"
        publishDir             "${resultsPath}/Aln", mode: 'link'
        stageOutMode           'copy'
        
        input:
        tuple val(id), file(reads) from trim_aln_ch
        file(idx) from bwa_index_ch_pe

        output:
        file("${id}.sorted.pe.bam*")
        file("*.log")
        
        script:
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa mem -t ${task.cpus - 2} \\
            \$INDEX \\
            $reads 2> ${id}.bwa.pe.log \\
            | samtools sort -@ 2 -T . -o ${id}.sorted.pe.bam -
        
        samtools index ${id}.sorted.pe.bam
        """
    }
}