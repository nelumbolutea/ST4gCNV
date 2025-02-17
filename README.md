# ST4gCNV

# The ST4gCNV pipeline leverages Next-Generation Sequencing (short reads) data (>=10x) to assess gene copy number varations between (homoploid) populations and closely related species.

# Step 1. Megahit for rapid assembly of NGS

##conda install -c bioconda megahit

megahit -1 input_1.fastq.gz -2 input_2.fastq.gz -o input -t 10  ###or use your customized parameters




# Step 2. Simulation of 10x reads based on genome assemblies
 
## Generate the adjacent 'windows' of megahit assembly 

samtools faidx input.final.megahit.contigs.fa

cut -f1,2 input.final.megahit.contigs.fa.fai > input.final.megahit.contigs.fa.sizes

bedtools makewindows -g  input.final.megahit.contigs.fa.sizes -w 200 -s 200 > input.final.megahit.contigs.windows200.bed

## Process the .bed file to merge the last window in a contig if it's shorter than 200bp
awk 'BEGIN {OFS="\t"} {contig=$1; start=$2; end=$3; if (contig != prev_contig) {if (NR > 1) {if (prev_end - prev_start < 200) {print prev_contig, prev_start - 200, prev_end} else {print prev_contig, prev_start, prev_end}} prev_contig=contig; prev_start=start; prev_end=end} else {if (end - start < 200) {prev_end=end} else {if (prev_start != "") {print prev_contig, prev_start, prev_end} prev_start=start; prev_end=end}}} END {print prev_contig, prev_start, prev_end}' input.final.megahit.contigs.windows200.bed > input.final.megahit.contigs.windows200_merged.bed

bedtools getfasta -fi input.final.megahit.contigs.fa -bed iput.final.megahit.contigs.windows200_merged.bed -fo input.final.megahit.contigs.windows200_merged.fasta

## create a new FASTA file where each sequence from the original is repeated 10 times##
awk '/^>/ {header=$0; next} {seq=$0; for(i=0;i<10;i++) {print header "_" i; print seq}}' input.final.megahit.contigs.windows200_merged.fasta> input.final.megahit.contigs.windows200_merged10x.fasta




# Step 3. Mapping simulated reads

bowtie2-build Reference_genome.fa Reference_genome

bowtie2 --local -x Reference_genome -f input.final.megahit.contigs.windows200_merged10x.fasta -S input.final.megahit.sam -p 20

samtools view -b input.final.megahit.sam -@ 10 -o input.final.megahit.bam

samtools sort -@ 10 -o input.final.megahit.sort.bam input.final.megahit.bam 

samtools index input.final.megahit.sort.bam



# Step 4. Gene CNV estimation based on simulated read coverages
## Simulated read depth on gene CDS
##Note that You can also utilize genomic coordinates (GFF or BED format) for any regions of interest, including pre-miRNAs, conserved non-coding regions, or even sliding windows across chromosomes for scanning CNVs. However, it is advisable to mask low-complexity and highly repetitive regions beforehand###

grep "CDS" Reference_genome.gff3 > Reference_genome.cds.gff

sed 's/\S\+Parent=//g' Reference_genome.cds.gff

##coverage per base###

bedtools coverage -a Reference_genome.cds.gff  -b input.final.megahit.sort.bam -d > input.final.megahit.sort.perbase.depth

##Aggregate and Calculate Average Coverage on CDS##

awk '{ sum[$9] += $11; count[$9]++ } END { for (id in sum) print id "\t" sum[id] / count[id] }' input.final.megahit.sort.perbase.depth > input.final.megahit.sort.perCDS.depth


# Step 5. BH-corrected t-tests for detecting gene CNV differentiation between populations in R platform

##suppose you have n samples, merge input_1.final.megahit.sort.perCDS.depth ... input_n.final.megahit.sort.perCDS.depth into a single tab-delimit table, for example "input_01_n.perCDS.depth.txt"

##Perform downstream t-tests with Benjamini-Hochberg (BH) correction between populations by running 'Bulk_t_test_gCNV.R' file###

Rscript Bulk_t_test_gCNV.R








