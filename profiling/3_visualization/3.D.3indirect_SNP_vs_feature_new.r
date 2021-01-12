library(data.table)
library(dplyr)

f <- list.files(path="/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/", pattern="indirect_", full.names = T)
fmic <- list.files(path="/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/", pattern="indirect_T[[0-9]", full.names = T)
id = which(f %in% fmic)
f = f[-c(id)]
id = grep("P368", f)
f = f[-c(id)]



out <- data.frame()
for(i in 1:length(f)){
  tem <- read.csv(f[i], header=TRUE)
  tem$file <- f[i]
  out <- rbind(out, tem)
}

out$file  <- gsub(".*\\/", "", out$file)
# method
out$method <- gsub("_.*", "", out$file)
# tissue type
out$tissue <- gsub(".*_ts_|.csv", "", out$file)
table(out$tissue)
# trait
out$trait <- gsub(".*indirect_|_ts.*", "", out$file)
table(out$trait)


snps_p7msnps <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/hmp321_282_agpv4_maf005_miss03_pruned.bim", header = FALSE , data.table=FALSE)
out <- merge(out, snps_p7msnps[,c(1,2,4)], by.x = "snps_for_medi", by.y = "V2", sort = F, all.x = T)
colnames(out)[c(8,9)] = c("chr", "pos")

write.table(out, "largedata/isnps_vs_feature/isnps_633005rows.csv", sep=",", row.names=FALSE, quote=FALSE)

isnps <- fread("largedata/isnps_vs_feature/isnps_633005rows.csv", header=TRUE, data.table=FALSE)

gff <- fread("largedata/isnps_vs_feature/Zea_mays.B73_RefGen_v4.46.chr.txt", header=FALSE, data.table=FALSE)
names(gff) <- c("seq", "source", "feature", "start", "end", "score", "strand", "phase", "att")
table(gff$feature)


### Get genes and upstream and downstream 5kb regions
g <- subset(gff, feature %in% "gene")
g$geneid <- gsub(".*gene:|;biotype.*", "", g$att)

### Get 5' utr and upstream 5kb regions
utr5 = subset(gff, feature %in% "five_prime_UTR")
utr5$geneid = gsub(".*transcript:|_T.*", "", utr5$att)
utr5$seq = as.integer(utr5$seq)

### + strand
gp <- subset(g, strand %in% "+") 


### get the 5k upstream of the + strand gene model
gp_up <- gp
gp_up$end <- gp_up$start - 1
gp_up$start <- gp_up$end - 5000 

### get the 5k downstream of the + strand gene model
gp_down <- gp
gp_down$start <- gp_down$end + 1
gp_down$end <- gp_down$start + 5000 


### - strand
gm <- subset(g, strand %in% "-") 


### get the 5k upstream of the - strand gene model
gm_up <- gm
gm_up$start <- gm_up$end + 1
gm_up$end <- gm_up$start + 5000 

### get the 5k downstream of the - strand gene model
gm_down <- gm
gm_down$end <- gm_down$start - 1
gm_down$start <- gm_down$end - 5000 



### + strand
utr5p <- subset(utr5, strand %in% "+") 


### get the 5k upstream of the + strand utr5
utr5p_up <- as.data.frame(utr5p %>% group_by(geneid) %>% summarise(seq = mean(seq), start = min(start), end = max(end)) )
utr5p_up$end <- utr5p_up$start - 1
utr5p_up$start <- utr5p_up$end - 5000 
utr5p_up$strand = "+"


### - strand
utr5m <- subset(utr5, strand %in% "-") 

### get the 5k upstream of the - strand utr5
utr5m_up <- as.data.frame(utr5m %>% group_by(geneid) %>% summarise(seq = mean(seq), start = min(start), end = max(end)) )
utr5m_up$start <- utr5m_up$end + 1
utr5m_up$end <- utr5m_up$start + 5000 
utr5m_up$strand = "-"



###output

gup <- rbind(gp_up, gm_up)
fwrite(gup, "largedata/isnps_vs_feature/gene_up5k.txt", sep="\t", row.names = FALSE, quote=FALSE)

gdown <- rbind(gp_down, gm_down)
fwrite(gdown, "largedata/isnps_vs_feature/gene_down5k.txt", sep="\t", row.names = FALSE, quote=FALSE)

fwrite(g, "largedata/isnps_vs_feature/gene.txt", sep="\t", row.names = FALSE, quote=FALSE)


utr5up = rbind(utr5p_up, utr5m_up)
fwrite(utr5up, "largedata/isnps_vs_feature/utr5up5k.txt", sep="\t", row.names = FALSE, quote=FALSE)


library("GenomicRanges")
library("plyr")

gup$start = as.integer(gup$start)
gup$end = as.integer(gup$end)

colnames(isnps)[8] = "seq"
colnames(isnps)[9] = "Pos"
isnps$seq = as.character(isnps$seq)

idx <- which(is.na(isnps$Pos))
isnps = isnps[-idx,]
isnps_uni = data.frame(snps_for_medi = isnps$snps_for_medi, seq = isnps$seq, Pos = isnps$Pos)
isnps_uni = unique(isnps_uni)

grc <- with(isnps_uni, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(gup, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))

### find overlaps between the two
tb <- findOverlaps(query=grf, subject=grc)
tb <- as.matrix(tb)

out1 <- as.data.frame(grf[tb[,1]])
out2 <- as.data.frame(grc[tb[,2]])
### for each genomic feature, find the sites with non-missing data
colnames(out2)[1:2] = c("seq", "pos")
out <- cbind(out1, out2[, 1:2]) 
out$seqnames = as.integer(as.character(out$seqnames))
out$seq = as.integer(as.character(out$seq))

isnps_uni$seq = as.integer(as.character(isnps_uni$seq))
isnps_uni$snps_for_medi = as.character(isnps_uni$snps_for_medi)

isnps$seq = as.integer(as.character(isnps$seq))

out$snps_for_medi = "not_exist"
for (i in 1:nrow(out)) {
  t = intersect(which(isnps_uni$seq == out$seq[i]), which(isnps_uni$Pos == out$pos[i])) 
  out$snps_for_medi[i] = isnps_uni$snps_for_medi[t]
}


fwrite(out, "largedata/isnps_vs_feature/isnps_in_gene_up5k.txt", sep="\t", row.names = FALSE, quote=FALSE)

fwrite(isnps, "largedata/isnps_vs_feature/isnps_47143rows.csv", sep=",", row.names=FALSE, quote=FALSE)



### null probability

colnames(snps_p7msnps)[c(1,2,4)] = c("seq", "snp", "Pos")
snps_p7msnps$seq = as.character(snps_p7msnps$seq)
grc1 <- with(snps_p7msnps, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(gup, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))

### find overlaps between the two
tb1 <- findOverlaps(query=grf, subject=grc1)
tb1 <- as.matrix(tb1)

out_p7m1 <- as.data.frame(grf[tb1[,1]])
out_p7m2 <- as.data.frame(grc1[tb1[,2]])
### for each genomic feature, find the sites with non-missing data
colnames(out_p7m2)[1:2] = c("seq", "pos")
out_p7m <- cbind(out_p7m1, out_p7m2[, 1:2]) 
out_p7m$seqnames = as.integer(as.character(out_p7m$seqnames))
out_p7m$seq = as.integer(as.character(out_p7m$seq))

snps_p7msnps$seq = as.integer(as.character(snps_p7msnps$seq))



print("p7msnps_gup5k")
print(nrow(unique(out_p7m[,7:8])))
print("p7msnps_gup5k_ratio")
print(nrow(unique(out_p7m[,7:8]))/nrow(snps_p7msnps))
