#!/usr/bin/env nextflow

params.genome                = false          /*genome fasta file, must specify complete path. Required parameters*/
params.samplePath            = false          /*input folder, must specify complete path. Required parameters*/
params.outputDir             = false          /*output folder, must specify complete path. Required parameters*/

readPrepPath                 = "${params.outputDir}/read_prep"
trimPath                     = "${params.outputDir}/trimmed"

myExecutor                   = 'slurm'
params.myQueue                      = 'normal'
defaultCPU                   = '1'
defaultMemory                = '20'
assemblerCPU                 = '12'
assemblerMemory              = '500'

params.samtoolsMod           = 'SAMtools/1.10-IGB-gcc-8.2.0'

/*Prepare input*/
genome_file                  = file(params.genome)
genomeStore                  = genome_file.getParent()
if( !genome_file.exists() ) exit 1, "Missing reference genome file: ${genome_file}"
CRAM_Ch1 = Channel.fromPath("${params.samplePath}")

/*

  prepare_genome 
  This process is executed only once
  NOTE: This process should remain the same as in the original workflow, will be merged in

*/


process prepare_genome{
    tag                    { genome }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    storeDir               genomeStore
    validExitStatus        0
    
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
   # bit flag 256 = not primary alignment (this removes these)
   samtools view -hb -f 12 -F 256 alignments.bam > both_unmapped.bam
    
   # R1 mapped, R2 not
   # bit flag 4   = R1 unmapped
   # bit flag 264 = Mate unmapped and not primary alignment (removes these). 
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 4 -F264 alignments.bam  > R1_unmapped.bam
   
   # R2 mapped, R1 not
   # bit flag 8   = R2 (mate) unmapped
   # bit flag 260 = R1 (read) unmapped and not primary alignment (removes these). 
   #                Note this is to make sure we're not keeping reads *also* in the first  
   #                set
   samtools view -hb -f 8 -F 260 alignments.bam  > R2_unmapped.bam

*/

/*

  extract_unmap 
  The input files are in CRAM format which have been previously aligned to genome prepared in previous step

*/

process extract_both_unmapped {
    tag                    { name }
    executor               myExecutor
    clusterOptions         params.clusterAcct 
    cpus                   defaultCPU
    queue                  params.myQueue
    memory                 "$defaultMemory GB"
    module                 params.samtoolsMod
    publishDir             'extracted', mode: "copy"	    
    validExitStatus        0,1
    errorStrategy          'finish'
    scratch                '/scratch'
    stageOutMode           'copy'
    
    input:
    set val(name), file(CRAM) from CRAM_Ch1	
    file genome from genome_file
    file index from genome_index_ch

    output:
    set val(name), file('*.both-unmapped.fastq') optional true into fq_paired_ch
    set val(name), file('*.SE-unmapped.fastq') optional true into fq_se_ch
    file '*'
    
    script:
    if(params.singleEnd){
    """
    # we only need to extract reads that are unmapped, no worries about pairing
    samtools view -@ ${defaultCPU} -hbt ${index} -f 4 ${name} > ${name.baseName}.unmap.bam
    samtools fastq -@ ${defaultCPU} -f 4 ${name.baseName}.unmap.bam > ${name.baseName}_SE_R1.fastq
    """
    } else {
    """
    # two stages; grab any non-properly paired reads (includes unmapped)
    samtools view -@ ${task.cpus} -hbt ${index} -G 2 -o ${name.baseName}.improper.bam ${name}
    
    # now capture the three classes; this is much faster than parsing the full BAM each time
    
    # both unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 12 -F 256 -o ${name.baseName}.both-unmapped.bam ${name.baseName}.improper.bam

    samtools fastq-@ ${task.cpus} ${name.baseName}.both-unmapped.bam \\
        -1 ${name.baseName}_PE_R1-both-unmapped.fastq -2 ${name.baseName}_PE_R2-both-unmapped.fastq
        
    # R1 unmapped
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 4 -F 264 -o ${name.baseName}.R1-unmapped.bam ${name.baseName}.improper.bam

    samtools fastq -@ ${task.cpus} ${name.baseName}.R1-unmapped.bam \\
        -1 ${name.baseName}_PE_R1-unmapped.fastq -2 ${name.baseName}_PE_R2-mapped.fastq

    # R2 unmapped    
    samtools view -@ ${task.cpus} -hbt ${index} \\
        -f 8 -F 260 -o ${name.baseName}.R2-unmapped.bam ${name.baseName}.improper.bam
        
    samtools fastq -@ ${task.cpus} ${name.baseName}.R1-unmapped.bam \\
        -1 ${name.baseName}_PE_R1-mapped.fastq -2 ${name.baseName}_PE_R2-unmapped.fastq
    
    # combined SE unmapped FASTQ into one file for assembly
    
    cat ${name.baseName}_PE_R1-unmapped.fastq ${name.baseName}_PE_R2-unmapped.fastq >> ${name.baseName}_SE-unmapped.fastq
    """    
    }

}

