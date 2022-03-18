/*
The logic in this step of the workflow follows the Chinese Han publication
logic for placing contigs up to clustering.  We'll include the relevant
python code in 'bin' but possibly start by using the nextflow script location
and maybe check in as a submodule for now.
*/

// Set up input channels from 'assembly' folder; the relevant files are in the 
// subfolders


// Step 1 involved assemblies and basic filtering
// Step 2: Get placement

process AlignToAssembly {

  // TODO: not sure why this switches to bowtie2 here from bwa mem in the initial alignment,
  // maybe strictness?
  module "bowtie2"

  input:
  // TODO: Maybe the input needs to be a CSV, or two channels with matching IDs? 
  // Here it's two channels joined by the ID
  tuple(id), file(ctg), file(fastq) from filteredContigs.join(unaln)

  output:
  tuple(id), file(sam) into alnToAssembly

  """
  bowtie2-build ${ctg} ${id} 

  # note this aligns the reads as unpaired, similar to the assembly; this would be modified for our workflow
  bowtie2 -x ${id} -U R1_alignedmate.fq, R2_alignedmate.fq  -S readtocontig.sam 
  """
}

process PlaceRegions {

  module "samtools", "bedtools"

  input:
  // TODO: alnToRef comes from the original alignment:
  // samtools view -f 8 -F 4 alignment.bam > alignedmate_GRCh38.sam
  // We could get this from the extracted discordant data 
  // from the assembly step
  tuple(id), file(aln), file from alnToAssembly.join(alnToRef)

  output

  """
  #
  samtools view -h -F 2304 readtocontig.sam  \
    | samtools sort -n -O bam \
    | bedtools bamtobed -i stdin \
    | awk '{OFS="\t"} {print \$4,\$1,\$6,\$2,\$3}' \
    | sed -e "s/\/[1-2]//g" \
    | sort > readtocontig.txt 

  # TODO: Split into another process

  samtools view -H alignedmate_GRCh38.sam \
    | cat - <(awk 'NR==FNR{ a[\$1]; next }\$1 in a{ print \$0 ; delete a[\$1]; next }' readtocontig.txt \
    <( samtools view alignedmate_GRCh38.sam )) \
    | samtools sort -n -O bam \
    | bedtools bamtobed -i stdin \
    | awk '{OFS="\t"}{print \$4,\$1,\$6,\$2,\$3}' \
    | sed -e "s/\/[1-2]//g" \
    | sort > pass_mates.txt 

  join -j 1 readtocontig.txt pass_mates.txt > mates_region.txt 
  """
}

process ExtractEnds {

  module "samtools"
  // TODO: notice here they use the no_alt reference (we use alts+decoys)
  """
  samtools faidx GRCh38_no_alt.fa unambiguou_placement > GRCh38_Region.fa 
  samtools faidx contig_ID.fa LEP/REP_region > LEP/REP_contig.fa
  """
}

process AlignContigsToRegion {

  module "mummer"

  """
  nucmer  --maxmatch -l 15 -b 1 -c 15 -p contig_ID GRCh38_Regions.fa REP_contig.fa/LEP_contig.fa  
  delta-filter -q -r -o 0 -g contig_ID.delta > filtered_info.delta 
  """
}

process ClassifyContigs {

  """
  python contig_type.py  --ref_name_id GRCH38.fa.fai --alignment_info PATH_filtered_info.delta  --LEP_contigs LEP_folder --REP_contigs REP_folder --BEP_contigs BEP_folder --BEP_contigs_all all_BEP_folder
  """
}

// Step 3: Cluster placed contigs

process ClusterPlacedContigs {

  module "bedtools"

  """
  # For BEP contigs, 
  awk '{OFS="\t"} {split(FILENAME,b,"."); if($7=="reverse") print $2,$3-1,$5,$1"_"b[1],"-";  else print $2,$3-1,$5,$1"_"b[1],"+"}' BEP_folder/* |bedtools sort -i > BEP_contigs.bed 

  # For LEP/REP contigs, 
  awk '{OFS="\t"} {split(FILENAME,b,"."); if($4=="reverse") print $2,$7-1,$8,$1"_"b[1],"-";  else print $2,$7-1,$8,$1"_"b[1],"+"}' LEP/REP_folder/* |bedtools sort -i > LEP/REP_contigs.bed

  bedtools merge -d 20 -c 4 -o distinct -i  placed_contigs.sorted.bed > merge_contigs.bed 
  """
}

process PickRepresentativeCluster {

  """
  python rep_obtain.py --seq_path LEP/REP/BEP_seq_path --path_merge_bed merge_contigs.bed --path_rep save_rep_folder --path_contig save_cluster_folder
  """
}

process RemoveNAContigs {

  module "mummer"
  """
  nucmer -p align_info  rep.fa cluster.fa<br>
  """
}

process AlignContigsToClusterContigs {

  module "blast"
  """
  makeblastdb -in remaining_cluster.fa -dbtype nucl -out remainingcontigs_Id 
  blastn -db remainingcontigs_Id -query othertype_contig.fa -outfmt "6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch gapopen gaps evalue bitscore" -max_target_seqs 1  -max_hsps 1  -out  othertype_contig.tsv 

  # TODO: Maybe split into a separate process
  # Two types of contigs
  awk '{OFS="\t"}{if(\$3>99 && (\$6-\$13)/\$4>=0.99 && (\$6-\$13) /\$5>=0.8 ) print \$2,\$1}' othertype_contig.tsv \
    > Ensure_contigs.txt
  awk '{OFS="\t"}{if(\$3>99 && (\$6-\$13)/\$5<0.8 && (\$6-\$13)/\$4>=0.99 ) print \$2,\$1}' othertype_contig.tsv \
    > candidate_contigs.txt

  # TODO: Maybe split into a separate process
  python move_contigs.py \
    --ensure_contigs Ensure_contigs.txt \
    --pass_contigs pass_contigs.txt \
    --cluster_folder remain_cluster_folder \
    --contig_path othertype_contig_path
  """
}

process AlignCloseContigs {

  module "mummer"

  """
  nucmer -f  -p align_info left_placed.fa  right_placed.fa 
  delta-filter -q  -r -g -m -1 align_info > filterdalign_info.delta 
  show-coords -H -T -l -c -o filterdalign_info.delta > filterdalign_info.coords  
  """
}

process UpdateAlnResult {



  """
  python reorg_align_info.py  --LEP_bed LEP_contigs.bed --REP_bed REP_contig.bed --Identity_path identity.coords --Contained_path contained.coords  --Overlap_path overlap.coords --Part_align_path part.coords --save_folder updated_alignment_folder
  """
}

process FilterFPMergingClusters {

  module "mummer"

  """
  nucmer -p Lrep_Rcluster  REP_cluster.fa LEP_rep.fa   
  nucmer -p Rrep_Lcluster LEP_cluster.fa  REP_rep.fa 
  delta-filter  -r -q -g LEP_rep_REP_cluster.delta > LEP_rep_REP_cluster_filter.delta 
  delta-filter  -r -q -g REP_rep_LEP_cluster.delta > REP_rep_LEP_cluster_filter.delta 
  show-coords -H -T -l -c -o LEP_rep_REP_cluster_filter.delta > LEP_rep_REP_cluster_filter.coords 
  show-coords -H -T -l -c -o REP_rep_LEP_cluster_filter.delta > REP_rep_LEP_cluster_filter.coords  
  """
}

process MergeLEP_REP {
  module "popins"
  """
  popins merge -c LEP_REP.fa <br>
  """
}

process UpdatePlacedReps {

  """
  python update_ref.py --LEP_folder LEP_folder/ --REP_folder  REP_folder/  --contigs_fai  contigs_fai_path  --LEP_cluster_folder LEP_cluster_folder --REP_cluster_folder  --LEP_rep_folder LEP_rep/ --REP_rep_folder REP_rep/ --align_folder updated_alignment_folder/  --LEP_cluster_update_folder LEP_cluster_update/ --LEP_rep_update_folder LEP_rep_update --REP_cluster_update_folder REP_cluster_update/ --REP_rep_update_folder REP_rep_update/ --BEP_rep_folder BEP_rep/ --BEP_cluster_folder BEP_cluster/ --merged_contig_folder merged/
  """
}

