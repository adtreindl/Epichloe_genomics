# Calculating summary statistics for Epichloe population data

We are using vcftools (v. 0.1.16) with an added --haploid switch to allow computation of pi and Fst and Tajima's D for haploid data.
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
# Calculating LD along genome

We split our vcf file into the different populations

```
for p in {popname}_{speciesname}.pop
do
 vcftools --vcf file.vcf --keep $p --recode --recode-INFO-all --out $p
done
```

We create map and ped files using plink (v. 1.90) and filter for minor allele frequency (maf 0.05)

```
for p in {popname}_{speciesname}.pop
do
 plink --vcf $p.recode.vcf  --allow-extra-chr --const-fid --keep-allele-order  --recode --out $p
 plink --file $p --allow-extra-chr --maf 0.05 --recode --out $p_maf05
done
```

We need to edit the map file to make the second column have seq numbers and chromosome names 1-7

```
cp $p_maf05.map $p_maf05.mapOLD
awk '$2=NR' OFS="\t" $p_maf05.mapOLD > $p_maf05.map

```
We calculate LD in 10kb windows
```
plink --file $p_maf05 --allow-extra-chr --r2 --ld-window-kb 10000 --ld-window-r2 0  --out $p_maf05
```
