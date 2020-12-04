---
title: "script_bioinformatics"
author: "Katja Kozjek"
---

# Clear environment
rm(list = ls())

# Prepare the environment
## load libraries
library(dada2); packageVersion("dada2") 
library(phyloseq); packageVersion("phyloseq") 

# Start with DADA2 pipeline, create ASVs from fastq files
## set path

list.files(path_fastq) #listing all files in the path
fns <- sort(list.files(path_fastq, pattern=".fq", full.names = TRUE)) #selection of all fastq files
sample.names <- sapply(strsplit(basename(fns), "_"), `[`, 2)
print(sample.names)

## visualize sequence quality 

plotQualityProfile(fns[1:6])
plotQualityProfile(fns[7:12])
plotQualityProfile(fns[13:18])
plotQualityProfile(fns[19:24])
plotQualityProfile(fns[25:30])
plotQualityProfile(fns[31:36])
plotQualityProfile(fns[37:42])
plotQualityProfile(fns[43:48])

## length distribution of the original sequences 

#inspect length distribution of sequences to determine how to filter the sequences
length_fns <- lapply(path_fastq, function(fn) nchar(getSequences(fn)))
length <- do.call(c, length_fns)
hist(length, 100)
mean(length); median(length) 

## filter and trim original sequences 

#place filtered files in filtered/ subdirectory 
filt_path <- file.path(path_fastq, "filtered") 
filt_fns <- file.path(path_fastq, "filtered", paste0(sample.names, "_filt.fq.gz"))
names(filt_fns) <- sample.names

#filter from 1400 to 1650
filt_out <- filterAndTrim(fns, filt_fns, minQ=2, minLen=1450, maxLen=1650, maxN=0, rm.phix=FALSE, maxEE=3, multithread = FALSE) #we retain only reads with fewer than three exp. errors 
print(filt_out)

#csv file containing number of input and output/filtered sequences
write.csv(filt_out, "results/filtered_seq.csv") 

## inspect length distribution of filtered sequences

length_fns_filt <- lapply(filt_path, function(fn) nchar(getSequences(fn)))
length_filt <- do.call(c, length_fns_filt)
hist(length_filt, 100)
mean(length_filt); median(length_filt)

## dereplication of identical reads
#dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence
#collapsing together all reads that encode the same sequence

derep <- derepFastq(filt_fns, verbose=TRUE)
names(derep) <- sample.names #name the derep-class objects by the sample names

## error rates for each possible transition

set.seed(100)
err <- learnErrors(derep, errorEstimationFunction = PacBioErrfun, BAND_SIZE = 32, verbose = TRUE)  

plotErrors(err)
plotErrors(err, nominalQ=TRUE) #red line shows the error rates, the estimated error rates (black line) are a good fit to the observed rates (points), error rates should drop with the increased quality 

## denoise 
#infer the sequence variants in each sample
#separate true biological sequences from errors

dada_dd <- dada(derep, err = err, BAND_SIZE = 32, multithread=TRUE, verbose = TRUE, pool = "pseudo") #number of true sequences
dada_dd[[1]] #check only the first sample

# Construct sequence variant (ASV) table, a higher-resolution version of the OTU table produced by traditional methods

seqtab <- makeSequenceTable(dada_dd); dim(seqtab) 
table(nchar(getSequences(seqtab))) #distibution of sequence lengths 

#save R object containing all original sequences
saveRDS(seqtab, "data/rds/original_seq.rds")

# extract original sequences

asv_seqs_original <- colnames(seqtab)
asv_headers_original <- vector(dim(seqtab)[2], mode = "character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers_original[i] <- paste(">ASV", i, sep = "_")
}
asv_fasta_original <- c(rbind(asv_headers_original, asv_seqs_original))
write(asv_fasta_original, "data/fasta/original_seq.fasta")

## remove chimeras 

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
#chimeras removed 
#seqtab.nochim file is the final file, samples and sequences that represents ASVs

saveRDS(seqtab.nochim, "data/rds/nochim_seq.rds")

#what proportion of the sequence variants are chimeras
sum(seqtab.nochim)/sum(seqtab) 
1-sum(seqtab.nochim)/sum(seqtab)

## Summary, track retained sequences

getN <- function(x) sum(getUniques(x))
track <- cbind(filt_out, sapply(dada_dd, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "nonchim")
rownames(track) <- sample.names
print(track)
class(track) #it is a matrix

#save summary 
write.csv(track, "results/summary_seq.csv")

#plot the distribution of sequence length 
plot(sort(unname(rowSums(seqtab.nochim))))

#species richness
rarecurve(seqtab.nochim)

# Export fasta file 

#give our sequence headers more manageable names (ASV_1, ASV_2...) and write out a FASTA of our ASV seqs
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode = "character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))

#export fasta file with processed ASVs
write(asv_fasta, "/data/fasta/processed_seq.fasta")

# Export count table

asv_tab <- (seqtab.nochim)
colnames(asv_tab) <- sub(">", "", asv_headers)

write.table(asv_tab, "data/csv/ASVs_counts.csv", sep = ",", quote = FALSE, col.names = NA)

#matrix: rows corresponding to samples, columns with ASV, the value of each entry is the number of times that ASV was observed in the sample

# Taxonomy with UNITE 2020

set.seed(NULL) #return to the original state by unsetting the seed
set.seed(62)

#take UNITE database

taxonomy_unite <- assignTaxonomy(seqtab.nochim, "reference databases/sh_general_release_dynamic_04.02.2020_dev.fasta", multithread=TRUE, tryRC = TRUE)

tax_unite2020 <- taxonomy_unite # Removing sequence rownames for display only
rownames(tax_unite2020) <- NULL
rownames(tax_unite2020) <- getSequences(asv_tab)
head(tax_unite2020)

write.csv(tax_unite2020, "data/csv/tax_unite2020.csv")

#add sequences
write.csv(taxonomy_unite, "data/csv/tax_unite2020_seq.csv")
