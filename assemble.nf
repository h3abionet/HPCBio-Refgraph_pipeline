#!/usr/bin/env nextflow

/*parameters that are specified at the command line or via config file*/
params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = "./results"    /*output folder, must specify complete path. Required parameters*/
params.singleEnd             = false          /*options: true|false. true = the input type is single end reads; false = the input type is paired reads. Default is false*/
params.assembler             = 'masurca'      /*options: megahit|masurca. Default is megahit*/
params.skipKraken2           = true           /*options: true|false. Default is true which means that kraken2 will be skipped*/

/* parameters for readprep = qctrimming and adapter removal */
params.skipTrim              = false           /*qc-trimming of reads. options: true|false. Default is false*/     
params.min_read_length       = '20'            /*minimum length of read to be kept after trimming for downstream analysis. Default is 20*/
params.min_base_quality      = '20'            /*minimum base quality. Default is 20*/
params.guess_adapter         = true            /*auto-detect adapter from input file. options: true|false. Default is true*/
params.forward_adapter       = false           /*adapter sequence to be clipped off (forward). */
params.reverse_adapter       = false           /*adapter sequence to be clipped off (reverse). Used for paired reads only*.*/

/*Stage*/
stage = "assembly"

/*output folder paths*/
readPrepPath                 = "${params.outputDir}/read_prep"
trimPath                     = "${params.outputDir}/trimmed"
megahitPath                  = "${params.outputDir}/megahit"
masurcaPath                  = "${params.outputDir}/masurca"
multiqcPath                  = "${params.outputDir}/multiqc"
metricsPath                  = "${params.outputDir}/assembly_metrics"

// Moving this to annotation
// kraken2Path                  = "${params.outputDir}/kraken2"

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

/*
  qc_input 
  This process mainly checks the inputs for improperly formed CRAM/BAM input files
*/

process qc_input {
    tag                    { id }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    errorStrategy          { task.exitStatus=1 ? 'ignore' : 'terminate' }
    
    input:
    set val(id), file(CRAM) from CRAM_Ch1

    output:
    set val(id), file('*_ok.cram') optional true into extract_unmapped_ch,extract_clipped_ch
    
    script:
    """
    samtools quickcheck ${CRAM}
    if [ \$? -eq 0 ]
    then
        cp ${CRAM} ${id}_ok.cram
    fi
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
   pull out subsets using flags. Otherwise we're running through the 
   
   # both unmapped; -f means unmapped, 
   # bit flag 12  = both reads unmapped (bit flag 4 & 8)
   # bit flag 2304 = not primary alignment (this removes these), not supplemental
   samtools view -hb -f 12 -F 2304 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = R1 unmapped
   # bit flag 2312 = Mate unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 4 -F 2312 alignments.bam  > R1_unmapped.bam
   
   # R2 mapped, R1 not
   # bit flag 8    = R2 (mate) unmapped
   # bit flag 2308 = R1 (read) unmapped and not primary alignment (removes these), not supplemental
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 8 -F 2308 alignments.bam  > R2_unmapped.bam
*/

/*
  extract_unmapped
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step
*/

process extract_improper {
    tag                    { id }
    executor               myExecutor
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             readPrepPath
    
    input:
    set val(id), file(cram) from extract_unmapped_ch
    file genome from genome_file
    file index from genome_index_ch

    output:
    // all unaligned + mates
    set val(id), file("${id}.improper.bam") into improper_unmapped_ch, improper_clipped_ch
    
    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority and the logic for extracting discordant
    // reads is an issue w/ SE data
    
    script:
    if(params.singleEnd) {
    
    """
    echo "SE support NYI"
    exit 1

    # Uncomment below when we test SE
    ## TODO: UNTESTED!!!!
    ## we only need to extract reads that are unmapped, no worries about pairing
    # samtools view -@ ${task.cpus} -hbt ${index} -f 4 -o ${id}.unmapped.bam ${cram} 
    """
    
    } else {
    
    """
    # grab any non-properly paired reads (includes any unmapped and discordant reads) 
    samtools view -@ ${task.cpus} -hbt ${index} -G 2 -o ${id}.improper.bam ${cram}
    """
    }
}

process extract_unmapped {
    tag                    { id }
    executor               myExecutor
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             readPrepPath
    
    input:
    set val(id), file(bam) from improper_unmapped_ch

    output:
    set val(id), file("${id}.all-unmapped.R{1,2}.fastq.gz") optional true into fq_pe_unmapped_ch
    set val(id), file("${id}.unmapped.bam") optional true into unmapped_bam_ch  // all unaligned + mates

    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority yet.
    
    script:
    if(params.singleEnd) {
    
    """
    echo "SE support NYI"
    exit 1

    # Uncomment below when we test SE
    ## TODO: UNTESTED!!!!
    ## convert to FASTQ
    # samtools fastq -@ ${task.cpus} ${bam} > ${id}.orphans.unmapped.fastq
    """
    
    } else {
    
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
}

// TODO: note the used of the hard-coded '/scratch' space here.  
// We could default to the current work directory, but the current behavior of 
// samtools collate is to use /tmp and does *not* currently use TMPDIR (yeah, pretty crappy)

process extract_clipped {
    tag                    { id }
    executor               myExecutor
    cpus                   12
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "SAMtools/1.12-IGB-gcc-8.2.0","Python/3.7.2-IGB-gcc-8.2.0"
    publishDir             readPrepPath, mode: "copy"    
    stageOutMode           'copy'
    
    input:
    set val(id), file(bam) from improper_clipped_ch

    output:
    set val(id), file("${id}.clipped.bam") into extract_sort_ch // all clipped pairs; we save this for now
    set val(id), file("${id}.all-clipped.R{1,2}.fastq.gz") into fq_pe_clipped_ch

    script:
    if(params.singleEnd) {
    """
    echo "SE support NYI"
    exit 1
    """
    } else {
    """
    # capture only discordant clipped reads; note the -G 2
    # the below is to prevent collisions if the pipeline is interrupted and restarted
    tmpdir=\$( mktemp -d /scratch/collate.XXXXXXXXX )
    samtools collate --output-fmt bam -@ ${task.cpus - 4} -O ${bam} \$tmpdir/${id} | \
        clipped-filter.py > ${id}.clipped.tmp.bam

    samtools fastq -@ ${task.cpus} ${id}.clipped.tmp.bam \
        -1 ${id}.all-clipped.R1.fastq.gz -2 ${id}.all-clipped.R2.fastq.gz
    
    samtools sort -@ ${task.cpus} -o ${id}.clipped.bam ${id}.clipped.tmp.bam
    """
    }
}

process merge_pairs {
    tag                    { id }
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "seqkit/0.12.1"
    publishDir             readPrepPath
 
    input:
    set val(id), file(unmapped), file(clipped) from fq_pe_unmapped_ch.join(fq_pe_clipped_ch)

    output:
    set val(id), file("${id}.all-reads.R{1,2}.fastq.gz") into merge_trim_ch, merge_fastqc_ch
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
    tag "FASTQC-Pretrim ${id}"
    executor               myExecutor
    cpus                   2
    queue                  params.myQueue
    memory                 "12 GB"
    module                 "FastQC/0.11.8-Java-1.8.0_152"
    publishDir             "${params.outputDir}/FASTQC-Pretrim"

    input:
    set val(id), file(reads) from merge_fastqc_ch

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
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    publishDir             trimPath
    module                 "fastp/0.20.0-IGB-gcc-4.9.4"

    input:
    set val(name), file(reads) from merge_trim_ch
    
    output:
    set val(name), file('*.PE.R{1,2}.trimmed.fastq.gz'), file('*.unpR{1,2}.trimmed.fastq.gz') optional true into trim_megahit_ch, trim_masurca_ch, trim_fastqc, trim_aln_ch
    set val(name), file('*.json') optional true into trim_multiqc_ch
    file '*.html'
        
    script:
    trimOptions      = params.skipTrim ? ' ' :  ' -l 20 --cut_right --cut_right_window_size 3 --cut_right_mean_quality 10'
    adapterOptionsSE = params.guess_adapter ? ' ' : " --adapter_sequence=${params.forward_adapter} "
    adapterOptionsPE = params.guess_adapter ? ' --detect_adapter_for_pe ' : " --adapter_sequence=${params.forward_adapter}  --adapter_sequence_r2=${params.reverse_adapter} "

    // TODO: leaving this switch in but note that none of the samples in the workflow are SE data;
    // we can prep for this but it's not high priority yet.

    if(params.singleEnd){
    """
    echo "SE support NYI"
    exit 1
    # fastp --in1 ${reads[0]} \
    #   --out1 "${name}.SE.R1.trimmed.fastq.gz" \
    #   ${adapterOptionsSE} ${trimOptions} \
    #   --thread ${task.cpus} \
    #   -w ${task.cpus} \
    # --html "${name}"_SE_fastp.html \
    # --json "${name}"_SE_fastp.json
    """
    
    } else {
    
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
}

process fastqc_post {
    tag "FASTQC-Posttrim ${id}"
    executor               myExecutor
    clusterOptions         params.clusterAcct     
    cpus                   4
    queue                  params.myQueue
    memory                 "12 GB"
    module                 "FastQC/0.11.8-Java-1.8.0_152"
    publishDir             "${params.outputDir}/FASTQC-Posttrim"

    input:
    set val(id), file(pereads), file(sereads) from trim_fastqc

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
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   24
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 "MEGAHIT/1.2.9-IGB-gcc-8.2.0" 
    publishDir             megahitPath
    
    input:
    set val(name), file(pefastqs), file(sefastqs) from trim_megahit_ch

    output:
    set val(name), val("megahit"), file("${name}/final.contigs.fa") into megahit_rename_ch
    file("${name}/*")

    script:
    if(params.singleEnd){
    """
    echo "SE support NYI"
    exit 1
    # megahit -1 ${fastqs[0]} -o ${name}.megahit_results
    # 
    # perl $params.assemblathon ${name}.megahit_results/final.contigs.fa > ${name}.megahit_results/final.contigs.fa.stats
    """
    } else {
    """
    megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]} \
        -r ${sefastqs[0]},${sefastqs[1]} \
        -o ${name}
    # megahit -1 ${pefastqs[0]} -2 ${pefastqs[1]}  -o ${name}.megahit_results
    """
    }
}

/*
  masurca
*/

process masurca_assemble {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct
    cpus                   16
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    module                 "MaSuRCA/3.4.2-IGB-gcc-8.2.0"
    publishDir             masurcaPath

    input:
    set val(name), file(pefastqs), file(sefastqs) from trim_masurca_ch

    output:
    set val(name), val("masurca"), file("${name}/CA/final.genome.scf.fasta") into masurca_rename_ch
    file("${name}/*")

    script:
    """
    mkdir ${name}
    cd ${name}

    cat << EOF > ${name}.masurca_config_file.txt
    DATA
    PE = pe 150 50 ../${pefastqs[0]} ../${pefastqs[1]}
    PE = s1 150 50 ../${sefastqs[0]}
    PE = s2 150 50 ../${sefastqs[1]}
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
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   1
    queue                  params.myQueue
    memory                 "$assemblerMemory GB"
    // TODO: a base perl install is fine, but we have this as a placeholder
    // just in case to remind us for Docker/Singularity
    // module                 "Perl/5.24.1-IGB-gcc-4.9.4"
    publishDir             "${params.outputDir}/FinalAssemblies/${assembler}"

    input:
    set val(name), val(assembler), file(assembly) from all_assemblies_rename_ch

    output:
    set val(name), val(assembler), file("${name}.${assembler}.final.fasta") into all_assemblies_metrics_ch,all_assemblies_aln_ch

    script:
    // Strip off anything after the first '.'
    shortname = name.replaceAll(~/\.\S+$/, "")
    """
    perl -p -e 's/^>(\\N+)/>${name}:\$1/' ${assembly} > ${name}.${assembler}.final.fasta
    """
}

// all_assemblies_metrics_ch = megahit_metrics_ch.mix(masurca_metrics_ch)

process assembly_metrics {
    tag {id}
    executor               myExecutor
    clusterOptions         params.clusterAcct     
    cpus                   4
    queue                  params.myQueue
    memory                 "12 GB"
    module                 "quast/5.0.0-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             "${params.outputDir}/QUAST/${assembler}"

    input:
    set val(id), val(assembler), file(asm) from all_assemblies_metrics_ch

    output:
    file "${id}/*" 
    file "${id}.${assembler}.report.txt" into metrics_multiqc_ch

    script:
    """
    quast.py \\
        --output-dir ${id} \\
        --threads ${task.cpus} \\
        ${asm}
    ln -s ${id}/report.txt ${id}.${assembler}.report.txt
    """
}

// not sure this will work
process aln_reads {
    tag {id}
    executor               myExecutor
    cpus                   8
    queue                  params.myQueue
    clusterOptions         params.clusterAcct     
    memory                 "12 GB"
    module                 "BWA/0.7.17-IGB-gcc-8.2.0","SAMtools/1.12-IGB-gcc-8.2.0"
    publishDir             "${params.outputDir}/BWA-MEM/${assembler}"

    input:
    set val(id), val(assembler), file(assembly), file(pereads), file(sereads) from all_assemblies_aln_ch.combine(trim_aln_ch, by:0)

    output:
    file "${id}.${assembler}.sorted.pe.bam*"
    file "${id}.${assembler}.sorted.se.bam*"
    file "${id}.${assembler}.sorted.merged.bam*"

    script:
    """
    bwa index ${assembly}

    # PE
    bwa mem -t ${task.cpus - 4} ${assembly} ${pereads[0]} ${pereads[1]} | samtools sort -@ 4 -o ${id}.${assembler}.sorted.pe.bam
    samtools index ${id}.${assembler}.sorted.pe.bam

    # SE
    cat ${sereads[0]} ${sereads[1]} > ${id}.se.fastq.gz
    bwa mem -t ${task.cpus - 4} ${assembly} ${id}.se.fastq.gz | samtools sort -@ 4 -o ${id}.${assembler}.sorted.se.bam
    samtools index ${id}.${assembler}.sorted.se.bam

    # merge both files
    samtools merge -@ ${task.cpus} ${id}.${assembler}.sorted.merged.bam ${id}.${assembler}.sorted.pe.bam ${id}.${assembler}.sorted.se.bam
    samtools index ${id}.${assembler}.sorted.merged.bam    
    """
}


process MultiQC {
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   2
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 "MultiQC/1.7-IGB-gcc-4.9.4-Python-3.6.1"
    publishDir             multiqcPath
 
    input:
    file('./QUAST/*') from metrics_multiqc_ch.collect().ifEmpty([])
    file('./FASTQC-Pretrim/*') from fastqc_results.collect().ifEmpty([])
    file('./FASTQC-Posttrim/*') from fastqc_trimmed_results.collect().ifEmpty([])
    file('./FASTP/*') from trim_multiqc_ch.collect().ifEmpty([])

    output:
    file "multiqc*"

    """
    multiqc . 
    """
} 
