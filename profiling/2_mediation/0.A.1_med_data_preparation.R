setwd('/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020')
library(data.table)



rna <- fread("/common/jyanglab/shared/dbcenter/Kremling_Nature3RNASeq282_March2018/Expression_matrix/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm_rounded_origNames_and_Altnames.txt", data.table=FALSE)
gs <- subset(rna, Tissue %in% c("GShoot"))
dim(gs) # 295 37136
# remove the duplicated ones
gs0 <- gs[!duplicated(gs$HMP32Name), ]
dim(gs0) #278 37136

#determine the first column of the gene expression
idx1 <- which(names(gs0) == "AC148152.3_FG001")

X <- gs0[, c(6,10:ncol(gs0))]
X$HMP32Name <- gsub("282set_", "", X$HMP32Name)


variances<-apply(X[,2:ncol(X)], 2, var)
X <- X[, -(which(variances<=0.05)+1)] # 278 25038
means<-apply(X[,2:ncol(X)], 2, mean)
X <- X[, -(which(means<=1)+1)] #278 19357


p <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/geno_282/allchr_bisnp_n282_snpid_maf01_geno2_pruned_NA_0_matrix.txt", data.table=FALSE)
p <- p[!is.na(p$GDDDaystoSilk),] #274 111963


pc <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/geno_282/allchr_bisnp_n282_snpid_maf01_geno2_pruned.eigenvec", data.table = FALSE)
dim(pc) #277

X2 <- merge(pc[, 2:4], X[, 1:3],  by.x="V2", by.y="HMP32Name")
df <- merge(X2, p[, 1:3], by.x="V2", by.y="Genotype")
dim(df)#250
df <- df[order(df$V2),]
gid <- df$V2


p <- p[order(p$Genotype),]            
y <- p[p$Genotype %in% gid, 1:2]
sum(y$Genotype != gid) # must be 0
y <- as.matrix(y[,-1])

row.names(y) <- gid

### Z
Z <- as.matrix(p[p$Genotype %in% gid, -1:-2])
row.names(Z) <- gid

### X
X <- X[order(X$HMP32Name),]
X <- X[X$HMP32Name %in% gid, ]
sum(X$HMP32Name != gid) # must be 0
X <- as.matrix(X[, -1])
row.names(X) <- gid

### X0
pc <- pc[order(pc$V2), ]
pc <- pc[pc$V2 %in% gid, -1]
sum(pc$V2 != gid) # must be 0
X0 <- pc[, 2:4]
row.names(X0) <- gid

###remove genes with expression correlation with phenotype lower than 0.05
cors<-apply(X, 2, function(X) cor(X, y))
X  <- X[, -(which((cors<=0.05) & (cors>=-0.05) ))] #255 12570

library(dplyr)
X <- as.data.frame(X)
X <- mutate_all(X, scale)
row.names(X) <- gid
X <- as.matrix(X)

X0 <- as.matrix(X0)



fwrite(as.list(colnames(X)), "largedata/geno/colnamesX_gs.csv", sep=",", row.names=FALSE, quote=FALSE)
fwrite(as.list(colnames(Z)), "largedata/geno/colnamesZ_gs.csv", sep=",", row.names=FALSE, quote=FALSE)




