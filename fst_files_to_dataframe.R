# make a combined table with the Fst estimates files from each population pair
# we are using R version 3.4.0

filelist <- read.table("{list_of_filenames}")

# we use one of the files to create a dataframe with the chromsome information and position of the loci	

loci <- read.table("{fst_filename}", header=TRUE)
loci <- loci[,1:2] # CHROM, POS

# now we add the data from the population files
df.fst <- loci 
for (i in filelist$V1){
	fstdat <- read.table(paste(i, sep=""), header=TRUE)
	df.fst <- cbind(df.fst,fstdat$WEIR_AND_COCKERHAM_FST)
	}

# all columns have the same name now so make sure to name the columns before writing out the table
#pop.nm <- gsub(".pop.weir.fst","",as.character(filelist$V1))
#names(df.fst)[3:30] <- pop.nm

write.table(df.fst, file="{filename}.txt", quote=FALSE, row.names=FALSE) 
