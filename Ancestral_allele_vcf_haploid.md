# Setting the ancestral allele in a vcf file with haploid data
This document describes how to annotate a vcf file to indicate which allele is the ancestral allele. 

The appraoch closely follows instructions posted here for diploid data
http://wasabiapp.org/vbox/data/ngg2016/21/Day2Session1.Genomicalignmentancestralalleles.html

You need a basic working knowledge of the command line and experience with vcftools and bcftools would be beneficial.

You need a biallellic SNP vcf file (remove indels and non-biallelic SNPs)
e.g biallelic SNPs only

```
vcftools --vcf file.vcf  --max-alleles 2 --recode --recode-INFO-all --out file_BI
```
e.g SNPs only
```
vcftools --vcf file_BI.recode.vcf  ---remove-indels --recode --recode-INFO-all --out file_BI_SNPS
```

Then create a vcf file for the outgroup individual only
```
vcftools --vcf file_BI_SNPS.recode.vcf --indv outgroupIDname --recode --recode-INFO-all --out outgroupIDname
```
Create a table of SNP positions and alleles using bcftools 
```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' outgroupIDname.recode.vcf  > file.tab
```
Create a table file with ancestral allele information.
Here we select the reference allele (REF) or the variant allele (ALT) and put it in the 5th column. If the site is not covered by our outgroup reads, i.e. missing '.', then we record it as missing.

for Haploid    
```
awk '{OFS="\t";if($5=="0"){print $1,$2,$3,$4,$3} \
	if($5=="1"){print $1,$2,$3,$4,$4} \
	if($5=="."){print $1,$2,$3,$4,$5}}' file.tab > file_aa.tab
```

Compress and index the table file and the original vcf file
```
bgzip file_aa.tab
tabix -s1 -b2 -e2 file_aa.tab.gz
bgzip file_BI_SNPS.recode.vcf
```
Create an INFO file line for the new vcf file
```
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt
```

Use bcftools to annotate the vcf file with the ancestral allele information 
```
bcftools annotate -a file_aa.tab.gz \
 -c CHROM,POS,REF,ALT,INFO/AA -h hdr.txt -Oz \
 -o newfile_aa.vcf.gz file_BI_SNPS.recode.vcf.gz
```

Check that it has worked. There should be an info field AA
```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' newfile_aa.vcf.gz | less
```
Count how many sites have ancestral allele information
```
bcftools view -e 'INFO/AA=="."' newfile_aa.vcf.gz -H | wc -l
```