# Comparative analysis of genome structure and synteny 

This file provides example code lines of how we analyzed Epichloe genome sequences

## Identify genome compartments with ```OcculterCut```

We identified RIP-affected AT-rich regions using OcculterCut v1.1. The software requires a genome sequence fasta file and, by providing an annotation file (gff), can also output how many genes are located in which compartment respectively.

```
OcculterCut -f path_to_genome_fasta/Ety1756_Epichloe_typhina_1756_33930528_v4.fna -a path_to_annotation_gff/Epichloe_typhina.gff3
```


GC-content plots with cut-offs to classify genomic regions belonging to distinct compartments are generated using the software output with gnuplot

Example for E. typhina:

```
gnuplot
set xlabel "GC (%)"
set ylabel "Proportion of genome"
set sample 1000
set xrange[0:100]
set yrange[0:]
set boxwidth 1
set style fill solid
set key off
set style line 1 lt 1 lc rgb "#0000FF" lw 3
set style line 2 lt 2 lc rgb "#32CD32" lw 3
Cauchy(x,xo,wi) = (1./pi) * wi / ((x - xo)**2 + wi**2)
plot 'compositionGC.txt' w boxes, 0.314439*Cauchy(x, 24.393,1.69404) + 0.685561*Cauchy(x, 52.6862, 1.15478) ls 2
set yrange[0:GPVAL_Y_MAX]
set arrow from 37.5569,0 to 37.5569,GPVAL_Y_MAX nohead front ls 1
set terminal postscript eps enhanced color "Helvetica" 20
set output 'Et_plot.eps'
replot
```

## Whole-genome alignments to analyze synteny

We performed pairwise alignments of three genomes, the two new genomes and the E. festucae Fl1 reference genome, with minimap2 v2.12 using the preset option for cross-species full genome alignment.

```
# between E. clarkii and E. typhina
minimap2 -cx asm20 --cs path_to_genome_fasta/Ety1756_Epichloe_typhina_1756_33930528_v4.fna path_to_genome_fasta/Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna > EcEt_ali_asm20.paf 

# with Ef as the target
minimap2 -cx asm20 --cs   path_to_genome_fasta/EfFl1_v0.2.fna  path_to_genome_fasta/Ety1756_Epichloe_typhina_1756_33930528_v4.fna > EtEf_ali_asm20.paf 
minimap2 -cx asm20 --cs   path_to_genome_fasta/Ef/EfFl1_v0.2.fna  path_to_genome_fasta/Ecl1605_22_Epichloe_clarkii_1605_22_45692596_v2.fna> EcEf_ali_asm20.paf 
```
