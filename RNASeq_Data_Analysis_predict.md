# RNA seq data analysis 
This document contains a step-by-step account of how we analyzed the RNA seq data, predicted putative genes using the funannotate pipeline https://funannotate.readthedocs.io/en/latest/ and annotated putative proteins using Inerproscan, signalP and effectorP.

Here we have provided examples of code for each step. The code chunks will not run automatically "as is" it will need to be edited by the user.

Here is a list of the software used (not including all dependencies)
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
First we create a list of all the sample.gtf files to merge and then merege gfts using ```stringtie``` 
```
stringtie --merge Ec_transcript.list -o Ec_merged.gtf

```
# Predict and annotate putative genes
We used funnanotate pipeline to predict putative protiens and then annotated these protiens using Interproscan.

## Create a Repeat masked geneome file using ```funannotate mask```
First we shorten the chromosome names using ```funannotate sort```. Change "Ecl_1605_22_1" to "Chr1"
```
funannotate sort -i Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna  -b Chr -o Ecl1605_newScaffname.fna
funannotate mask -i Ecl1605_newScaffname.fna -o Ecl1605_newScaffname_FUNmasked_noLib.fna 
```

## Get gene predictons with ```funannotate predict```
First we converted the chromosome names back to original names on the funannotate masked file (FUNmasked.fna), so the information in the bams is matching the genome. Then we ran funnanotate predict using: a merged bam file (created using the mapped bams and samtools merge, e.g. samtools merge -b RNAseq_bam.list merged.out.bam), a gff file that was produced by another group (Winter et al) using the funnannotate pipeline and the merged transcripts file produced using stringtie (see above). This step outputs a gft and protiens.fa. 

```
bsub -W 24:00 -n 4 funannotate predict -i /path_to_file/Ecl1605_FUNmasked_noLib.fna \
-o /path_to_out_dir/ \
-s "Epichloe_clarkii" --name Ecl_ \
--rna_bam /path_to_bam/H_merged.out.bam \
--augustus_gff  /path_to_gff/Epichloe_clarkii_Hl.gff3 \
--stringtie /path_to_transcripts/Ec_merged.gtf \
--busco_db sordariomycetes \
--cpus 12

```

## Run ```Iterproscan``` to annotate the gene models
Using the protiens.fa files created in the predict step

```
interproscan.sh -dp -iprlookup --goterms --pathway \
	-b Ec_itps \
	-i /path_to_predict_results/Epichloe_clarkii.proteins.fa \
	-u /out_dir/
	

```
# Identify putative effector proteins
We used signalP and effectoP to identify possible effector protiens that can play an important role in plant-pathogen interactions. SignalP searches for signal peptides. The fasta output from signalP was used to run the online version of ```effectorP 2.0``` (<http://effectorp.csiro.au/>)  to find fungal effectors.

```
/path_to_signalp/4.1/signalp -m Ec_signalP.fa -n Ec_signalP.gff /path_to_predict_results/Epichloe_clarkii.proteins.fa > Ec_signalP.out

```

