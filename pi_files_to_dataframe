# make a combined table with the diversity estimates files from each population (Pi in 10kb windows)
# we are using R version 3.4.0

chr.names <- read.csv("chr_length", header=FALSE) # a file with the chromosome names and their length in bp
chr.names <- chr.names[-8,] # there are 7 chromosomes, we remove the mitochondrion (8)

filelist <- read.table("{list_of_filenames}")

# we create a dataframe with three columns: CHROM, BIN_START, BIN_END	

df.pi <- NULL
for(i in chr.names[,1]){
	ln <- chr.names[chr.names$V1==i,2]
	BIN_START <-  seq(1,ln, by=10000) # we want 10kb windows 
	BIN_END <- BIN_START+9999
	df <- data.frame(cbind(BIN_START,BIN_END))
	df$CHROM <- rep(i, dim(df)[1])
	df <- df[,c(3,1,2)]
	df.pi <- rbind(df.pi, df)
	}

# now we add the data from the population files
for (i in filelist$V1){
	pidat <- read.table(paste(i, sep=""), header=TRUE)
	pidat <- pidat[,-4]
	ff <- merge(df.pi,pidat, by = c("CHROM","BIN_START","BIN_END"), all.x=TRUE)
	df.pi <- cbind(df.pi,ff$PI)
	}

# all columns have the same name now so make sure to name the columns before writing out the table
write.table(df.pi, file="{filename}.txt", quote=FALSE, row.names=FALSE) 
