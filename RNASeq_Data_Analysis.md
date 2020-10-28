# RNA seq data analysis 
This doc contains a step-by-step account of the analysis of RNA seq data using the funannotate pipeline https://funannotate.readthedocs.io/en/latest/ 

Here we have provided examples of code for each step. The code chunks will not run automatically "as is" it will need to be edited by the user.

Here is a list of the software used (no including all dependencies)
funnanotate
fastqc v0.11.4
cufflinks
star v2.5.3a
stringtie v1.3.3b
repeatmasker v4.0.6 

## Check RNA sequence read quality with ```fastqc```
Example code

```
while read p; do
fastqc path_to_seq_data/RNAseq_${name}.fastq.gz -o fastqc_out/
done<RNAseq_sample.list

```
Sequences were not trimmed based on recommendations from this paper
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0956-2

## Generate files for ```star``` and ```funannotate```

Convert reference genome gff to gtf using ```cufflinks```
```
gffread /path_to_file/Epichloe_clarkii.gff3 -T -F -o Epichloe_clarkii.gtf

```
## Clean the genome
```
funannotate clean -i Ety1756_Epichloe_typhina_1756_33930528_v4.fna -o Ety1756_clean.fna 
```

## Generate genome files for ```star``` 

```
fasta="genome_sequence.fna"
gtf="file.gtf"
GDIR="new_directory"

STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $GDIR  \
--genomeFastaFiles $fasta \
--sjdbGTFfile      $gtf \
--sjdbOverhang 124

```

## Map RNAseq reads to the reference genome using ```star```
```
GDIR="path_to_genome_files"

STAR --runMode alignReads --runThreadN 20 --genomeDir $GDIR \
--readFilesIn /path_to_seq_data/RNAseq_${name}_R1.fastq.gz /path_to_seq_data/RNAseq_${name}_R2.fastq.gz \
--readFilesCommand zcat	\
--outFileNamePrefix ${name}_ \
--outSAMtype BAM SortedByCoordinate	\
--outSAMstrandField intronMotif 

```

## Make transcript asssembly (reference guided) using ```stringtie```
https://ccb.jhu.edu/software/stringtie/

```
while read name; do
stringtie /path_to_mapped_bam/${name}_Aligned.sortedByCoord.out.bam -o ${name}.gtf -G file.gtf 
done<RNAseq_sample.list

```
## Merge all sample gtf files to create a single non redundant file 
First create a list of all the sample.gtf files to merge and then merege gfts using ```stringtie``` 
```
stringtie --merge Ec_transcript.list -G /path_to_file/Epichloe_clarkii.gtf -o Ec_merged.gtf

```

## Compare the merged gtf with the ref genome gtf using ```gffcompare```
https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
```
./gffcompare -r /paht_to_file/Epichloe_clarkii.gtf /path_to_file/Ec_merged.gtf -o Compare_Ec_Ref
```

## Use funannotate to predict genes
Use RNAseq data for training and then use this to predict genes across the genome using funannotate

## Create a Repeat masked geneome file using ```repeatmasker```
First we shorten the chromosome names using ```funannotate sort```, then get Fungi TE library and run RepeatMasker
```
funannotate sort -i Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna  -b Chr -o Ecl1605_newScaffname.fna
queryRepeatDatabase.pl  -species Fungi > FungiDatabase.lib
RepeatMasker Ecl1605_newScaffname.fna -species fungi -dir /out_dir/

```

## Use ```funannotate train``` to perform a genome-guided Trinity RNA-seq assembly followed by PASA assembly
```
funannotate train -i Ecl1605_newScaffname.fna.masked -o ${TMPDIR}  \
    --left /path_to_RNAseq_data/RNAseq_sample1_R1.fastq.gz /path_to_RNAseq_data/RNAseq_sample2_R1.fastq.gz /path_to_RNAseq_data/RNAseq_sample3_R1.fastq.gz  \
    --right /path_to_RNAseq_data/RNAseq_sample1_R2.fastq.gz /path_to_RNAseq_data/RNAseq_sample2_R2.fastq.gz /path_to_RNAseq_data/RNAseq_sample3_R2.fastq.gz \
    --stranded RF --jaccard_clip --species "Epichloe clarkii" --no_trimmomatic  --no_normalize_reads --cpus 12
    
scp -r ${TMPDIR} ./Ec_train

```

## Get gene predictons with ```funannotate predict```
This creates gft and protiens.fa 
```
bsub -W 24:00 -n 4 funannotate predict -i Ecl1605_newScaffname.fna.masked \
	-o /cluster/scratch/stapleyj/Epichloe/funannotate/Ec_train -s "Epichloe_clarkii" --name Ecl_ \
	--transcript_alignments /paht_to_merged_transcripts/Ec_merged.gtf --cpus 12

```

## Run ```Iterproscan``` to annotate the gene models
using protiens.fa files created in previous step

```
interproscan.sh -dp -iprlookup --goterms --pathway \
	-b Ec_itps \
	-i /path_to_predict_results/Epichloe_clarkii.proteins.fa \
	-u /out_dir/
	

```
