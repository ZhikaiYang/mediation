setwd("/common/jyanglab/zhikaiyang/projects/mediation")
library(data.table)


isnps <- fread("largedata/isnps_vs_feature/isnps_55257rows.csv", header=TRUE, data.table=FALSE)

gff <- fread("largedata/isnps_vs_feature/Zea_mays.B73_RefGen_v4.46.chr.txt", header=FALSE, data.table=FALSE)
names(gff) <- c("seq", "source", "feature", "start", "end", "score", "strand", "phase", "att")
table(gff$feature)


### Get features


f_gene <- subset(gff, feature %in% "gene")
f_gene$geneid <- gsub(".*gene:|;biotype.*", "", f_gene$att)

#####change by features########################################################################################
int327shoot = fread("largedata/isnps_vs_feature/interactions_from_H3K27ac-ChIA-PET_in_shoot.txt", header=TRUE, data.table=FALSE)
#####change by features########################################################################################

int327shoot$start1 = as.integer(int327shoot$start1)
int327shoot$start2 = as.integer(int327shoot$start2)
int327shoot$end1 = as.integer(int327shoot$end1)
int327shoot$end2 = as.integer(int327shoot$end2)

fwrite(f_gene, "largedata/isnps_vs_feature/f_gene.txt", sep="\t", row.names = FALSE, quote=FALSE)


library("GenomicRanges")
library("plyr")
f_gene$start = as.integer(f_gene$start)
f_gene$end = as.integer(f_gene$end)

colnames(isnps)[8] = "seq"
colnames(isnps)[9] = "Pos"
isnps$seq = as.character(isnps$seq)

isnps_uni = data.frame(snps_for_medi = isnps$snps_for_medi, seq = isnps$seq, Pos = isnps$Pos)
isnps_uni = unique(isnps_uni)


#################################################################################
idx = which(is.na(int327shoot$start1))
int327shoot = int327shoot[-idx,]
idx = which(is.na(int327shoot$end1))

grc1 <- with(int327shoot, GRanges(seqnames=chr1, IRanges(start=start1, end=end1)))

grf <- with(f_gene, GRanges(seqnames=seq, IRanges(start=start, end=end), geneid=geneid))

### find overlaps between the two
tb <- findOverlaps(query=grf, subject=grc1)
tb <- as.matrix(tb)
tb1 = tb

grc2 <- with(int327shoot, GRanges(seqnames=chr2, IRanges(start=start2, end=end2)))
tb <- findOverlaps(query=grf, subject=grc2)
tb <- as.matrix(tb)
tb2 = tb


length(intersect(tb1[,2], tb2[,2]))
#13632
idx1 = which(!(tb1[,2] %in% intersect(tb1[,2], tb2[,2])))
enhancer_id2 = unique(tb1[,2][idx1])

idx2 = which(!(tb2[,2] %in% intersect(tb1[,2], tb2[,2])))
enhancer_id1 = unique(tb2[,2][idx2])

length(intersect(enhancer_id1, enhancer_id2))
# 0
tem1 =int327shoot[enhancer_id1,1:3]
tem2 = int327shoot[enhancer_id2,4:6]
colnames(tem1) = c("chr", "start", "end")
colnames(tem2) = c("chr", "start", "end")

enhancer = rbind(tem1, tem2)



grc <- with(isnps_uni, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(enhancer, GRanges(seqnames=chr, IRanges(start=start, end=end)))

### find overlaps between the two
tb <- findOverlaps(query=grf, subject=grc)
tb <- as.matrix(tb)



#out1 <- as.data.frame(grf[tb[,1]])
out2 <- as.data.frame(grc[tb[,2]])
### for each genomic feature, find the sites with non-missing data
colnames(out2)[1:2] = c("seq", "pos")
#out <- cbind(out1, out2[, 1:2]) 
out2$seqnames = as.integer(as.character(out2$seqnames))
#out$seq = as.integer(as.character(out$seq))

isnps_uni$seq = as.integer(as.character(isnps_uni$seq))
isnps_uni$snps_for_medi = as.character(isnps_uni$snps_for_medi)

isnps$seq = as.integer(as.character(isnps$seq))

out2$snps_for_medi = "not_exist"
for (i in 1:nrow(out2)) {
  t = intersect(which(isnps_uni$seq == out2$seq[i]), which(isnps_uni$Pos == out2$pos[i])) 
  out2$snps_for_medi[i] = isnps_uni$snps_for_medi[t]
}


fwrite(out2, "largedata/isnps_vs_feature/isnps_in_int327shoot.txt", sep="\t", row.names = FALSE, quote=FALSE)



### null probability
snps_p7msnps <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/gen/allchr_bisnp_n282_snpid_maf01_geno2_pruned.bim", header = FALSE , data.table=FALSE)

colnames(snps_p7msnps)[c(1,2,4)] = c("seq", "snp", "Pos")
snps_p7msnps$seq = as.character(snps_p7msnps$seq)
grcnull <- with(snps_p7msnps, GRanges(seqnames=seq, IRanges(start=Pos, end=Pos)))

grf <- with(enhancer, GRanges(seqnames=chr, IRanges(start=start, end=end)))

### find overlaps between the two
tbnull <- findOverlaps(query=grf, subject=grcnull)
tbnull <- as.matrix(tbnull)

#out_p7m1 <- as.data.frame(grf[tbnull[,1]])
out_p7m2 <- as.data.frame(grcnull[tbnull[,2]])
### for each genomic feature, find the sites with non-missing data
colnames(out_p7m2)[1:2] = c("seq", "pos")
#out_p7m <- cbind(out_p7m1, out_p7m2[, 1:2]) 
#out_p7m$seqnames = as.integer(as.character(out_p7m$seqnames))
out_p7m2$seq = as.integer(as.character(out_p7m2$seq))

#snps_p7msnps$seq = as.integer(as.character(snps_p7msnps$seq))


print("p7msnps_int327shoot")
print(nrow(unique(out_p7m2[,1:2])))
print("p7msnps_int327shoot_ratio")
print(nrow(unique(out_p7m2[,1:2]))/nrow(snps_p7msnps))




