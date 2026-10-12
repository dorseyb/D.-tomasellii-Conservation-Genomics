############################################################################################################################################################
# R script to filter data after process_radtags in Stacks and run analyses found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
#
# 2/8/24
# Modified for use with D. tomasellii samples
#############################################################################################################################################################



# Libraries ---------------------------------------------------------------
{ # this block loads all packages and functions
  
setwd("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/Rworking/Revision")

library(adegenet)
library(pcadapt)
library(qvalue)

medianFunc <- function(data, indices) {
  data <- as.matrix(data)
  dt <- data[indices,]
  med <- apply(dt,2,median, na.rm = T)
  return(med)
}

## Import bayescan R code
source("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/bayescan/plot_R.r")

#  From https://rdrr.io/github/romunov/zvau/src/R/writeGenPop.R
writeGenPop <- function(gi, file.name, comment) {
  
  if (is.list(gi)) {
    # do all genind objects have the same number of loci?
    if (length(unique(sapply(gi, nLoc))) != 1) stop("Number of loci per individual genind object in a list is not equal for all.")
    gi.char <- gi[[1]]
    loc.names <- locNames(gi[[1]])
  } else {
    gi.char <- gi
    loc.names <- locNames(gi)
  }
  
  # Calculate the length of two alleles.
  lng <- as.character(na.omit(genind2df(gi.char)[, locNames(gi.char)[1]]))
  lng <- unique(nchar(lng))
  
  stopifnot(length(lng) == 1)
  
  cat(paste(comment, "\n"), file = file.name)
  cat(paste(paste(loc.names, collapse = ", "), "\n"), file = file.name, append = TRUE)
  
  if (is.list(gi)) {
    pop.names <- seq_len(length(gi))
  } else {
    pop.names <- popNames(gi)
  }
  
  for (i in pop.names) {
    cat("pop\n", file = file.name, append = TRUE)
    if (is.list(gi)) {
      intm <- gi[[i]]
      loc.names <- locNames(gi[[i]])
    } else {
      intm <- gi[pop(gi) == i, drop = FALSE]
    }
    ind.names <- indNames(intm)
    intm <- genind2df(intm, sep = "")
    intm[is.na(intm)] <- paste(rep("0", lng), collapse = "")
    out <- cbind(names = paste(ind.names, ",", sep = ""), intm[, loc.names])
    write.table(out, file = file.name, row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
  }
  
  return(NULL)
}

source("/Users/bdorsey/Library/CloudStorage/OneDrive-TheHuntingtonLibrary,ArtMuseum,andBotanicalGardens/Documents/Rfunctions/genind2structure.R")
}


# Import data from Stacks ---------------------------------------------------

# D tom only
stackdir <- '/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullSet/M1n1/populations/mac3r5p6BL'

# With OG
stackdir <- '/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullsetOG/m7M1n1/populations/mac3r5p7BL'

hapsgenepop <- paste(stackdir, "populations.haps.genepop.gen", sep='/')
hapstsv <- paste(stackdir, "populations.haplotypes.tsv", sep='/')
snpsgenepop <- paste(stackdir, "populations.snps.genepop.gen", sep='/')

# haplotypes in genepop format
haps <- read.genepop(file=hapsgenepop)

# Check pop names vs levels and substitute  ---------------------------------------------------
levels(haps@pop)
# Dtom only
levels(haps@pop) <- c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad")

# with og
levels(haps@pop) <- c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad", "Dholm")

pops <- levels(haps@pop)

# snps in genepop format
snps <- read.genepop(file=snpsgenepop)
popNames(snps)

# Dtom only
popNames(snps) <- c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad")

# with og
levels(snps@pop) <- c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad", "Dholm")





mafAllhaps <- minorAllele(haps)
mafAllsnps <- minorAllele(snps)
length(which(mafAllsnps<0.01))
length(which(mafAllhaps<0.05))

which(mafAllhaps == min(mafAllhaps))

lowFreqSnps <- which(mafAllsnps < 0.05)


allPerLoc <- nAll(haps)
max(allPerLoc)
min(allPerLoc)
alHist <- hist(allPerLoc, breaks=(c(1:24)))
allGT2 <- which(allPerLoc > 2)
length(allGT2)
length(which(allPerLoc == 1))

allFreqs <- apply(tab(haps, freq=T),2,mean, na.rm=T)

allCounts <- apply(tab(haps, freq=F),2,sum, na.rm=T)
length(which(allCounts < 3))

snpAllCounts <- apply(tab(snps, freq=F),2,sum, na.rm=T)
length(which(snpAllCounts < 3))



tab(haps)[1:10,1:10]

# Filtering missing data  ---------------------------------------------------
# get proportion of missing data 
# across loci
haps.typedLoc <- propTyped(haps, by="loc")

# Filter by loci, keep with > 0.7 coverage
haps70.loci <- which(haps.typedLoc >= 0.7) # indices of loci with >=0.7
haps70 <- haps[,loc=haps70.loci] # new genind object

keep70loci <- locNames(haps70)

mafAllhaps70 <- minorAllele(haps70)
length(which(mafAllhaps70<0.01))
length(which(mafAllhaps70<0.05))

which(mafAllhaps70 == min(mafAllhaps70))

# pcadapt for non-neutral loci----------------------------------------------------------

# use plink to convert vcf to bed
# In terminal:
# create new file for plink:
#   sed 's/^\([0-9]\)/s\1/g' ./populations.snps.vcf > snps.forPlink.vcf
# 
# use in plink
# plink2 --vcf snps.forPlink.vcf \
# --allow-extra-chr \
# --max-alleles 2 \
# --set-all-var-ids @:# \
#   --make-bed \
# --out mac3r5p6BL.all.snps

# read in bed file for pcadapt
bedFile <- "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullSet/M1n1/populations/mac3r5p6BL/mac3r5p6BL.all.snps.bed"

#og bed file
bedFile <- "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullsetOG/m7M1n1/populations/mac3r5p7BL/mac3r5p6BL.all.snps.bed"

snps70.pcadapt <- read.pcadapt(bedFile, type = "bed")

# run with many pca axes to decide how many to use
ktest1 <- pcadapt(snps70.pcadapt, K=20, min.maf = 0)

# screeplot
plot(ktest1, option = "screeplot")

# score plot
plot(ktest1, option = "scores")
plot(ktest1, option = "scores", i=3, j=4)
plot(ktest1, option = "scores", i=5, j=6)
plot(ktest1, option = "scores", i=7, j=8)
plot(ktest1, option = "scores", i=9, j=10)

# Choose K=4

# set min maf
maf <- 0


# run with K=4
adapt.K4.res <- pcadapt(snps70.pcadapt, K=4, min.maf = maf)

# Check res
plot(adapt.K4.res, option = "manhattan")
plot(adapt.K4.res, option = "qqplot")
hist(adapt.K4.res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# choose outliers
# using qvalues
qval <- qvalue(adapt.K4.res$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha) # indices of outliers
length(outliers)

# select outlier snps
select.snps <- snps[,outliers]

# Get outlier names - as snps
select.snps.names <- locNames(select.snps)

# as loci
# Remove allele numbers
select.snps.LociNames <- gsub("_[0-9]+", "", x=select.snps.names)
head(select.snps.LociNames)

# Get unique loci
select.snps.LociNames <- sort(unique(select.snps.LociNames))
length(select.snps.LociNames)

# get names of neutral loci
NeutLoci <- setdiff(locNames(haps), select.snps.LociNames)
length(NeutLoci)

# get neutral snp names
neutral.snps.names <- locNames(snps[,-outliers])
head(neutral.snps.names)

# Final loci to keep  ---------------------------------------------------
NeutralKeepers <- intersect(keep70loci,NeutLoci)
haps70neut <- haps[loc=NeutralKeepers]

SelectKeepers <- intersect(keep70loci, select.snps.LociNames)
haps70Select <- haps[loc=SelectKeepers]
haps70maf05Select <- haps[loc=SelectKeepers]

# Convert and write files ---------------------------------------------------

# raw for finRAD
haps.raw <- read.delim(hapstsv, check.names = F)
# haps.raw2 <- haps.raw[which(haps.raw[,1] %in% locNames(haps),arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column

haps.raw.neut <- haps.raw[which(haps.raw[,1] %in% NeutralKeepers,arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column
# replace '-' with '' for missing data
haps.raw.neut[haps.raw.neut == '-'] <- ''

haps.raw.select <- haps.raw[which(haps.raw[,1] %in% SelectKeepers,arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column
haps.raw.select[haps.raw.select == '-'] <- ''

write.table(haps.raw.neut[,2:length(haps.raw.neut)], file = "Dtom.haps70.neut.frs.txt", quote = F, sep = "\t", na = "", row.names = F)

write.table(haps.raw.select[,2:length(haps.raw.select)], file = "Dtom.haps70.select.frs.txt", quote = F, sep = "\t", na = "", row.names = F)

# Genepop files

writeGenPop(haps70neut, "Dtom.haps70.neut.gen", "Dtom.filtered.neutral")

writeGenPop(haps70Select, "Dtom.haps70.select.gen", "Dtom.filtered.select")

writeGenPop(haps70maf05Select, "Dtom.haps70.maf05.select.gen", "Dtom.filtered.select")

haps70neut <- read.genepop("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/Rworking/Revision/data/Dtom.haps70.neut.gen")

# locus lists

write.table(NeutralKeepers, file = "neut70filterWL.txt", quote = F, row.names = F, col.names = F)

write.table(SelectKeepers, file = "select70filterWL.txt", quote = F, row.names = F, col.names = F)

write.table(select.snps.names, file = "select.snps.txt", quote = F, row.names = F, col.names = F)

write.table(neutral.snps.names, file = "neutral.snps.txt", quote = F, row.names = F, col.names = F)


# use WL to filter original snps using populations

# write structure file
snpfilt <- read.genepop("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullSet/M1n1/populations/mac3r5p6BL/revision/onesnp/populations.snps.genepop.gen")
genind2structure(snpfilt, "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/structure/revision/Dtom.neut.snps.structure")

# Go back and filter OG data

# OG set
# raw for finRAD
haps.raw <- read.delim(hapstsv, check.names = F)
# haps.raw2 <- haps.raw[which(haps.raw[,1] %in% locNames(haps),arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column

haps.raw.neut <- haps.raw[which(haps.raw[,1] %in% NeutralKeepers,arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column
# replace '-' with '' for missing data
haps.raw.neut[haps.raw.neut == '-'] <- ''

haps.raw.select <- haps.raw[which(haps.raw[,1] %in% SelectKeepers,arr.ind=T),c(1,3:length(haps.raw))] # remove filtered haps and count column
haps.raw.select[haps.raw.select == '-'] <- ''

write.table(haps.raw.neut[,2:length(haps.raw.neut)], file = "DtomOG.haps70.neut.frs.txt", quote = F, sep = "\t", na = "", row.names = F)

write.table(haps.raw.select[,2:length(haps.raw.select)], file = "DtomOG.haps70.select.frs.txt", quote = F, sep = "\t", na = "", row.names = F)

# Genepop files

writeGenPop(haps70neut, "DtomOG.haps70.neut.gen", "DtomOG.filtered.neutral")

writeGenPop(haps70Select, "DtomOG.haps70.select.gen", "DtomOG.filtered.select")

haps70neutOG <- read.genepop("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/Rworking/Revision/data/DtomOG.haps70.neut.gen")

# locus lists

write.table(NeutralKeepers, file = "OGneut70filterWL.txt", quote = F, row.names = F, col.names = F)

write.table(SelectKeepers, file = "OGselect70filterWL.txt", quote = F, row.names = F, col.names = F)

write.table(select.snps.names, file = "OGselect.snps.txt", quote = F, row.names = F, col.names = F)

write.table(neutral.snps.names, file = "OGneutral.snps.txt", quote = F, row.names = F, col.names = F)


