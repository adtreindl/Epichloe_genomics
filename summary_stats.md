# Calculating summary statistics for Epichloe population data

We are using vcftools with an added --haploid switch to allow computation of pi and Fst and Tajima's D for haploid data.
Details are here: https://github.com/vcftools/vcftools/pull/69

# What files we are working with

We have a filtered vcf file containing all individuals (populations) for each of two focal species with biallelic, haploid SNPs. 
We also have files with the sample names for every population within each species to split the vcf file by populations. They all end in {popname}_{speciesname}.pop

# calculating nucleotide diversity Pi

We are using the --window-pi function from vcftools (haploid mode) with a window sike of 10kb
```
for p in {popname}_{speciesname}.pop
do
 vcftools --vcf file.vcf --haploid --keep $p --window-pi 10000 --out $p
done
```
