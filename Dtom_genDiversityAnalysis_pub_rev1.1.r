############################################################################################################################################################
# R script to filter data after process_radtags in Stacks and run analyses found in Dorsey et al 2023 Conservation genomics of Dioon holmgrenii (Zamiaceae) reveals
# a history of range expansion, fragmentation, and isolation
# of populations. Conservation Genetics
#
# 2/8/24
# Modified for use with D. tomasellii samples
#############################################################################################################################################################



# Libraries ---------------------------------------------------------------
{ # this block loads libraries and functions
  
setwd("/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/Rworking/Revision")
library(adegenet)
library(pegas)
library(ggplot2)
library(ggrepel)
library(hierfstat)
library(poppr)
library(mmod)

medianFunc <- function(data, indices) {
  data <- as.matrix(data)
  dt <- data[indices,]
  med <- apply(dt,2,median, na.rm = T)
  return(med)
}

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

}
# Filtered Data -----------------------------------------------------------

dir <- getwd()

# read from filtered files
# Neutral Loci
h <- paste(dir,"data/Dtom.haps70.neut.gen", sep = "/")

# random neutral snps
s <- "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/stacks/fullSet/M1n1/populations/mac3r5p6BL/revision/onesnp/populations.snps.genepop.gen"

haps <- read.genepop(h)

levels(haps@pop)
# "DtoNmr-9" "DtoMva-9" "DtoReu-9" "DtoSdq-9" "DtoPmo-9" "DtoHxa-9"

pops <- c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad")
popsOG <- c(pops, "Dholm")

# Dtom only
levels(haps@pop) <- pops

haps

popNames(haps)

# convert to genpop object
haps.pop <- genind2genpop(haps)

# convert to heirfstat
haps.hier <- genind2hierfstat(haps)

snps <- read.genepop(s)
levels(snps@pop)

separatePops <- seppop(snps)
snpIn <- repool(separatePops[1:6])
levels(snpIn@pop)


## Non-Neutral Loci

hs <- paste(dir,"data/Dtom.haps70.select.gen", sep = "/")

hs <- paste(dir,"data/Dtom.haps70.maf05.select.gen", sep = "/")

haps <- read.genepop(hs)

# DAPC --------------------------------------------------------------------
## Neutral Loci -----------------------------------------------------------

dat <- tab(haps, freq=T, NA.method="mean")

# find clusters
grp <- find.clusters(dat, n.start=100, scale=F, n.pca=200, n.grp=4)
# 200 pc's kept
# 4 clusters chosen

#xvalidation to determine proper # pc to retain
crossvalTest <- xvalDapc(dat, grp$grp, n.pca.max = 200, center = T, scale = F, n.pca = NULL, n.rep = 30, xval.plot = T)

crossvalTest[6]
# $`Number of PCs Achieving Lowest MSE`
# [1] "140"

# dapc
haps.da <- dapc(dat, n.pca=140, n.da=3, grp=grp$grp, var.loadings=T)

## Non-neutral Loci -----------------------------------------------------------

dat <- tab(haps, freq=T, NA.method="mean")

# find clusters
grp <- find.clusters(dat, n.start=100, scale=F)
# 200 pc's kept
# 4 clusters chosen

#xvalidation to determine proper # pc to retain
crossvalTest <- xvalDapc(dat, grp$grp, n.pca.max = 200, center = T, scale = F, n.pca = NULL, n.rep = 30, xval.plot = T)

crossvalTest[6]
# $`Number of PCs Achieving Lowest MSE`
# [1] "40"

# maf05 select loci
haps.da <- dapc(dat, n.pca=40, n.da=3, grp=grp$grp, var.loadings=T)

summary(haps.da)



# Plot DAPC ---------------------------------------------------------------
dev.off()
co <- data.frame(haps.da$ind.coord, cluster=haps.da$grp, pop=haps$pop)
# co$pop <- haps$pop

cents <- as.data.frame(ddply(co, ~pop, summarize, mean1=mean(LD1), mean2=mean(LD2), mean3=mean(LD3)))

grp.coords <- as.data.frame(haps.da$grp.coord)



# ld12 <- ggplot(co) + geom_point(aes(LD1, LD2, color=pop)) + geom_text_repel(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop) + theme(legend.position = "none") + geom_text_repel(data=grp.coords,  x=grp.coords$LD1, y=grp.coords$LD2, label=rownames(grp.coords))

ld12 <- ggplot(co) + geom_point(aes(LD1, LD2, color=pop)) + stat_ellipse(aes(x=LD1,y=LD2, group = cluster)) + geom_text_repel(data=cents, x=cents$mean1, y=cents$mean2, label = cents$pop) + theme(legend.position = "none") + geom_text_repel(data=grp.coords,  x=grp.coords$LD1, y=grp.coords$LD2, label=rownames(grp.coords))

ld23 <- ggplot(co) + geom_point(aes(LD2, LD3, color=pop)) + stat_ellipse(aes(x=LD2,y=LD3, group = cluster)) + geom_text_repel(data=cents, x=cents$mean2, y=cents$mean3, label = cents$pop) + theme(legend.position = "none") + geom_text_repel(data=grp.coords,  x=grp.coords$LD2, y=grp.coords$LD3, label=rownames(grp.coords))

pdf("Dtom.select.dapc2.pdf",8,12)
# pdf("Dtom.dapc.pdf",8,12)

arrplot <- grid.arrange(ld12,ld23)

dev.off()

# pdf("Dtom.compoplat.pdf",12,8)
pdf("Dtom.compoplat.pdf",12,8)
compoplot(haps.da)
dev.off()

# pdf("Dtom.assignplot.pdf",8,8)
pdf("Dtom.assignplot.pdf",8,8)
assignplot(haps.da)
dev.off()


# c("Navoa", "Miravalles", "Rancho Eureka", "Sierra de Quila", "Pedro Moreno", "Soledad")

# Apply cluster results to Genind object
# Check cluster levels in plot if redoing this as they are arbitrary
# cluster.levels <- c("Sierra De Quila", "Jalisco West", "Nayarit", "Soledad") # haps
# cluster.levels <- c("Sierra De Quila","Nayarit", "Jalisco West", "Soledad") # snps

clusters <- grp$grp
clusters <- snps@pop
levels(clusters)[snps@pop["DtoHxa-9"]] <- "Soledad"
levels(clusters)[snps@pop["DtoNmr-9"]] <- "Western Jalisco"
levels(clusters)[snps@pop["DtoMva-9"]] <- "Nayarit"
levels(clusters)[snps@pop["DtoSdq-9"]] <- "Sierra de Quila"

levels(clusters) <- c("Santuario", "Nayarit", "Nayarit", "Sierra de Quila", "Santuario", "Soledad")

levels(clusters)[grp$grp["DtoHxa-1"]] <- "Soledad"
levels(clusters)[grp$grp["DtoNmr-10"]] <- "Western Jalisco"
levels(clusters)[grp$grp["DtoMva-10"]] <- "Nayarit"
levels(clusters)[grp$grp["DtoSdq-1"]] <- "Sierra de Quila"

strata.df <- data.frame(local=haps@pop, cluster=clusters)
strata.df <- data.frame(local=snps@pop, cluster=clusters)

# redefine populations by dapc results
strata(haps) <- strata.df
setPop(haps) <- ~cluster

strata(snps) <- strata.df
setPop(snps) <- ~cluster

popNames(haps)

# write files
# haps neut clust
writeGenPop(haps, "Dtom_haps_neut_clusters.gen", "Dtom neut dapc clusters")
# snps neut clust
writeGenPop(snpIn, "Dtom_snps_neut_clusters.gen", "Dtom nuet snps dapc clusters")
# haps select clust
writeGenPop(hapsSelect, "Dtom_haps_maf05.select_clusters.gen", "Dtom select haps dapc clusters")

# load from file (if needed)
h <- paste(dir,"/data/Dtom_haps_neut_clusters.gen", sep="")

haps <- read.genepop(h)
popNames(haps)

# pops <- cluster.levels # defined above
popNames(haps) <- c("WJalisco", "Nayarit", "Sierra de Quila", "Soledad")
strata(haps) <- haps@pop


# convert formats
haps.pop <- genind2genpop(haps)
haps.hier <- genind2hierfstat(haps)


# filter loci with >4 alleles for BA3
BAloci <- which(haps@loc.n.all < 5)

hapsBA <- haps[loc=BAloci]

writeGenPop(hapsBA, "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/ba3/revision/Dtom_clust_BA3.gen", "Dtom loci with max 4 alleles")

hba <- "/Volumes/CatDisk/Dioon_RADseq_analysis/Dtom/ba3/Dtom_clust_BA3.gen"

bagen <- read.genepop(hba)

length(which(bagen@loc.n.all > 4))

# Within pop diversity-------------

## Private Alleles ---------------------------------------------------------
# poppr

haps <- snps

priv.alls.dose <- poppr::private_alleles(haps)
# priv.alls.dose[1:4, 1:5]
#total count of private allele copies per population
private.count <- rowSums(priv.alls.dose)

#count of alleles that are private to each population
priv.alls <- poppr::private_alleles(haps,count.alleles = F)
# priv.alls[,1:10]
priv.alls.sums <- rowSums(priv.alls)

# number of distinct alleles present in each pop
distinctallelesInpops <- apply(haps.pop@tab,1,function(l){length(l[l != 0])})

# percentage of alleles present in each pop that are private
perc.priv <- round((priv.alls.sums/distinctallelesInpops*100),2)

# number of individuals per population
Indcounts <- summary(haps@pop)

# possible number of genotypes per population
all.poss <- Indcounts*nLoc(haps)

# number of missing genotypes per individual
n.na.gt <- apply(is.na(haps.hier),1,sum)

# number of missing genotypes per pop
pop.na.gt <- as.numeric(t(rowsum(n.na.gt,haps.hier$pop)))

# number of present genotypes
actual.genos <- all.poss-pop.na.gt

# actual number of alleles (copies) per pop
total.alls <- actual.genos*2

# percentage of copies that are private alleles
perc.private.copies <- as.numeric(round((private.count/total.alls)*100,2))
# as proportion (freq)
freq.private.copies <- as.numeric(round((private.count/total.alls),4))
names(perc.private.copies) <- names(perc.priv)

# put into table
private.allele.df <- data.frame(names(private.count),Indcounts,actual.genos,total.alls,distinctallelesInpops,private.count,perc.private.copies,priv.alls.sums,perc.priv)

pa.df.sort <- arrange(private.allele.df,perc.priv)

# write.csv(pa.df.sort, file = "Dtom.privateAlleles.csv")
# write.csv(pa.df.sort, file = "DtomHWfilt.privateAlleles.csv")
# write.csv(pa.df.sort, file = "Dtom_HWfilt_clusters_privateAlleles.csv")
write.csv(pa.df.sort, file = "Dtom_neut_clusters_privateAlleles.csv")
# write.csv(pa.df.sort, file = "Dtom_HWfilt_clusters_privateAlleles.csv")
write.csv(pa.df.sort, file = "Dtom_select_clusters_privateAlleles.csv")

# non-neutral loci
write.csv(pa.df.sort, file = "Dtom_selectMaf05_pop_privateAlleles.csv")
write.csv(pa.df.sort, file = "Dtom_selectMaf05_clust_privateAlleles.csv")

## Basic Stats -------------------------------------------------------------

# hier.fstat for diversity

basics <- basic.stats(haps, diploid = T, digits = 3)

basics$overall

write.csv(basics$overall, file = "Dtom_neut_clust_basics_overall.csv")

write.csv(basics$overall, file = "Dtom_selectMaf05_clust_basics_overall.csv")

arHier <- allelic.richness(haps)
apply(arHier$Ar,2, mean, na.rm=T)

basics.pop.meds <- matrix( c(
  apply(basics$Fis,2,median, na.rm=T),
  apply(basics$Ho,2,median, na.rm=T),
  apply(basics$Hs,2,median, na.rm=T), 
  apply(arHier$Ar,2, median, na.rm=T) ),
  nrow = 4, ncol = length(popNames(haps)), byrow = T, dimnames = list(c("Fis", "Ho", "Hs", "Ar"), popNames(haps))
)

write.csv(basics.pop.meds, file = "Dtom_neut_clust_basic_meds.csv")

write.csv(basics.pop.meds, file = "Dtom_selectMaf05_clust_basic_meds.csv")

# bootstraps
FisBootHier <- boot(basics$Fis, medianFunc, R=10000, parallel = "multicore")
HoBootHier <- boot(basics$Ho, medianFunc, R=10000, parallel = "multicore")
HsBootHier <- boot(basics$Hs, medianFunc, R=10000, parallel = "multicore")
ArBootHier <- boot(arHier$Ar, medianFunc, R=10000, parallel = "multicore")

FhierCIlist <- list()
for (i in 1:length(levels(haps@pop))) {
  FhierCIlist[[i]] <- boot.ci(FisBootHier, index = i)
}
allFisCI <- do.call(rbind, (lapply(1:length(FhierCIlist), function (x) FhierCIlist[[x]]$normal)))

# write.csv(allFisCI, file="Dtom_neut_FisCI.csv")
# write.csv(allFisCI, file="Dtom_HWfilt_FisCI.csv")
write.csv(allFisCI, file="Dtom_neut_clustFisCI.csv")
# write.csv(allFisCI, file="Dtom_HWfilt_clustFisCI.csv")
write.csv(allFisCI, file="Dtom_selectMaf05_clustFisCI.csv")


HohierCIlist <- list()
for (i in 1:length(levels(haps@pop))) {
  HohierCIlist[[i]] <- boot.ci(HoBootHier, index = i)
}
allHoCI <- do.call(rbind,(lapply(1:length(HohierCIlist), function (x) HohierCIlist[[x]]$normal)))

write.csv(allHoCI, file="Dtom_neut_clustHoCI.csv")
write.csv(allHoCI, file="Dtom_selectMaf05_clustHoCI.csv")

HshierCIlist <- list()
for (i in 1:length(levels(haps@pop))) {
  HshierCIlist[[i]] <- boot.ci(HsBootHier, index = i)
}
allHsCI <- do.call(rbind,(lapply(1:length(HshierCIlist), function (x) HshierCIlist[[x]]$normal)))

write.csv(allHsCI, file="Dtom_neut_clustHsCI.csv")
write.csv(allHsCI, file="Dtom_selectMaf05_clustHsCI.csv")


ArhierCIlist <- list()
for (i in 1:length(levels(haps@pop))) {
  ArhierCIlist[[i]] <- boot.ci(ArBootHier, index = i)
}
allArCI <- do.call(rbind,(lapply(1:length(ArhierCIlist), function (x) ArhierCIlist[[x]]$normal)))

write.csv(allArCI, file="Dtom_neut_clustArCI.csv")
write.csv(allArCI, file="Dtom_selectMaf05_clustArCI.csv")


# Between Pop Diversity-------------

hierwc <- wc(haps)
setPop(haps) <- ~cluster
popNames(haps)
## mmod for Gst and D---------------------

mmodDiff <- diff_stats(haps)
mmodDiff$global
pwD <- pairwise_D(haps)
pwGstP <- pairwise_Gst_Hedrick(haps)
pwGstN <- pairwise_Gst_Nei(haps)

## poppr for amova----------------

amova.result <- poppr.amova(haps, ~cluster)
amova.sig <- randtest(amova.result, nrepet = 999)
amova.sig
plot(amova.sig)
