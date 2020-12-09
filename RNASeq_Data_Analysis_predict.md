# RNA seq data analysis 
This document contains a step-by-step account of the analysis of RNA seq data using the funannotate pipeline https://funannotate.readthedocs.io/en/latest/ and annotation of putative proteins using Inerproscan, signalP and effectorP.

Here we have provided examples of code for each step. The code chunks will not run automatically "as is" it will need to be edited by the user.

Here is a list of the software used (no including all dependencies)
funnanotate v1.7.0, fastqc v0.11.4, cufflinks, star v2.5.3a, stringtie v1.3.3b, repeatmasker v4.0.6 , signalP v4.1, effectorP v2.0

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
stringtie --merge Ec_transcript.list -o Ec_merged.gtf

```

## Use funannotate to predict genes
Use RNAseq data for training and then use this to predict genes across the genome using funannotate

## Create a Repeat masked geneome file using ```repeatmasker```
First we shorten the chromosome names using ```funannotate sort```. Change "Ecl_1605_22_1" to "Chr1"
```
funannotate sort -i Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna  -b Chr -o Ecl1605_newScaffname.fna
funannotate mask -i Ecl1605_newScaffname.fna -o Ecl1605_newScaffname_FUNmasked_noLib.fna 
```

## Get gene predictons with ```funannotate predict```
First convert the chromosome names back to original names on the funannotate masked file (FUNmasked.fna), so the information in the bams is matching the genome.
Then run predict with as follows, this creates gft and protiens.fa 
A merged bam file eas created using the mapped bams and samtools merge (e.g. samtools merge -b RNAseq_bam.list merged.out.bam).
The gff file was produced by another group using (Winter et al) using the funnannotate pipeline.

```
bsub -W 24:00 -n 4 funannotate predict -i /cluster/work/gdc/shared/p427/GenAccFiles/FunMask/Ecl1605_FUNmasked_noLib.fna \
-o /cluster/scratch/stapleyj/Epichloe/funannotate/Ec \
-s "Epichloe_clarkii" --name Ecl_ \
--rna_bam /cluster/work/gdc/shared/p427/RNAseq/merged_bams/H_merged.out.bam \
--augustus_gff  /cluster/work/gdc/shared/p427/Genomes/Ec/Epichloe_clarkii_Hl.gff3 \
--stringtie /cluster/work/gdc/shared/p427/RNAseq/transcripts/Ec_merged_noG.gtf \
--busco_db sordariomycetes \
--cpus 12

```

## Run ```Iterproscan``` to annotate the gene models
using protiens.fa files created in previous step

```
interproscan.sh -dp -iprlookup --goterms --pathway \
	-b Ec_itps \
	-i /path_to_predict_results/Epichloe_clarkii.proteins.fa \
	-u /out_dir/
	

```
## Identify putative effector proteins
We used signalP and effectoP to identify possible effector protiens that may play an important role in plant-pathogen interaction. SignalP searches for signal peptides. The fasta output from signalP was used to run the online version of ```effectorP 2.0```<http://effectorp.csiro.au/>  to find fungal effectors.

```
/path_to_signalp/4.1/signalp -m Ec_signalP.fa -n Ec_signalP.gff /path_to_predict_results/Epichloe_clarkii.proteins.fa > Ec_signalP.out

```

