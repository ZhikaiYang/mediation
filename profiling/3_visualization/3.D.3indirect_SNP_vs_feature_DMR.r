setwd("/common/jyanglab/zhikaiyang/projects/mediation")
library(data.table)


isnps <- fread("largedata/isnps_vs_feature/isnps_55257rows.csv", header=TRUE, data.table=FALSE)

dmr <- fread("largedata/isnps_vs_feature/DMR.txt", header=T, data.table=FALSE)
colnames(dmr)[2:4] = c("chr", "start", "end")
dmr$chr = as.character(dmr$chr)

library("GenomicRanges")
library("plyr")

colnames(isnps)[8] = "seq"
colnames(isnps)[9] = "Pos"
isnps$seq = as.character(isnps$seq)

isnps_uni = data.frame(snps_for_medi = isnps$snps_for_medi, seq = isnps$seq, Pos = isnps$Pos)
isnps_uni = unique(isnps_uni)

isnps_uni$seq = as.character(isnps_uni$seq)
isnps_uni$snps_for_medi = as.character(isnps_uni$snps_for_medi)


grc <- with(isnps_uni, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(dmr, GRanges(seqnames=chr, IRanges(start=start, end=end)))

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


fwrite(out, "largedata/isnps_vs_feature/isnps_in_dmr.txt", sep="\t", row.names = FALSE, quote=FALSE)



### null probability
snps_p7msnps <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/gen/allchr_bisnp_n282_snpid_maf01_geno2_pruned.bim", header = FALSE , data.table=FALSE)

colnames(snps_p7msnps)[c(1,2,4)] = c("seq", "snp", "Pos")
snps_p7msnps$seq = as.character(snps_p7msnps$seq)
grc1 <- with(snps_p7msnps, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(dmr, GRanges(seqnames=chr, IRanges(start=start, end=end)))

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


print("p7msnps_dmr")
print(nrow(unique(out_p7m2[,1:2])))
print("p7msnps_dmr_ratio")
print(nrow(unique(out_p7m2[,1:2]))/nrow(snps_p7msnps))
