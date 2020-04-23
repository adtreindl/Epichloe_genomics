# Calculating summary statistics for Epichloe population data

We are using vcftools with an added --haploid switch to allow computation of pi and Fst and Tajima's D for haploid data.
Details are here: https://github.com/vcftools/vcftools/pull/69

# What files we are working with

We have a filtered vcf file containing all individuals (populations) for each of two focal species with biallelic, haploid SNPs. 
We also have files with the sample names for every population within each species to split the vcf file by populations. They all end in {popname}_{speciesname}.pop

# Calculating nucleotide diversity Pi

We are subsetting the main vcf file for each popualtion and using the --window-pi function from vcftools (haploid mode) with a window size of 10kb.
Calcualting Pi per site (--site-pi) doesn't seem to work because it doesn't take into account the monomorphic sites (values are too high).
```
for p in {popname}_{speciesname}.pop
do
 vcftools --vcf file.vcf --haploid --keep $p --window-pi 10000 --out $p
done
```

# Calculating Weir and Cockerhams Fst 

We are using the --weir-fst-pop function from vcftools (haploid mode) to calculate Fst for every site in every population pair
```
for pop1 in $(find {popname}_{speciesname}.pop)
do
    for pop2 in $(find {popname}_{speciesname}.pop)
    do
        vcftools --vcf file.vcf --haploid --weir-fst-pop ${pop1}  --weir-fst-pop ${pop2}  --out ${pop1}_${pop2}
    done
done
```
This calculates Fst for every pairwise comaprison so there will be a file named popX_popY and one names popY_popX. Obviously estimates should be exactly the same and we only need to keep one file at the end.
This also outputs a .log file for every comparison from which we can get the mean and mean weighted Weir and Cockerham Fst estimate. We can also calculate it ourselves though.

# Calculating Tajima's D

We did this for two different window sizes: 10kb and 40kb

```
for p in {popname}_{speciesname}.pop
do
 vcftools --vcf file.vcf --haploid --keep $p --TajimaD 10000 --out 10kb_$p
done

for p in {popname}_{speciesname}.pop
do
 vcftools --vcf file.vcf --haploid --keep $p --TajimaD 40000 --out 10kb_$p
done
```
