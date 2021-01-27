setwd("/common/jyanglab/zhikaiyang/projects/mediation") 
library(circlize)
library(data.table)
library(Ropt)
library(dplyr)
data<-read.table("./largedata/Chromosome_v4.txt",head=T,stringsAsFactors=FALSE,sep='\t') 

p_chr <- vector()
p_start <- vector()
p_end <- vector()
for (i in 1:nrow(data)) {
  p_chr <- c(p_chr, rep(data[i,1],50))
  t <- sample(c(data[i,2]:(data[i,3] - 4000000)), 50)
  p_start <- c(p_start, t ) 
  p_end <- c(p_end, (t + 4000000))
}
p_lines <- data.frame(chr= p_chr, start = p_start, end = p_end, p_5 = rep(5, 50*nrow(data)), p_10 = rep(10, 50*nrow(data)) )


V3_V4 <- fread("./largedata/V3_V4.csv", header = FALSE , data.table=FALSE)
V3_V4 <- V3_V4[,c(1,2)]
Zea_mays.AGPv4_anno <- fread("./largedata/Zea_mays.AGPv4_anno.txt", header = TRUE, data.table = FALSE)
Zea_mays.AGPv4_anno <- data.frame(chr = Zea_mays.AGPv4_anno$chromosome_name, gene4 = Zea_mays.AGPv4_anno$ensembl_gene_id_v4, start4 = Zea_mays.AGPv4_anno$start_position, end4 = Zea_mays.AGPv4_anno$end_position, name = Zea_mays.AGPv4_anno$description)
Zea_mays.AGPv34_anno <- merge(V3_V4, Zea_mays.AGPv4_anno, by.x = "V2", by.y = "gene4")
idx <- which ( Zea_mays.AGPv34_anno$chr %in% c(1:10))
Zea_mays.AGPv34_anno <- Zea_mays.AGPv34_anno[idx,]
Zea_mays.AGPv34_anno <- arrange(Zea_mays.AGPv34_anno, chr, start4)
Zea_mays.AGPv34_anno <- data.frame(chr = Zea_mays.AGPv34_anno$chr, gene3 = Zea_mays.AGPv34_anno$V1, gene4 = Zea_mays.AGPv34_anno$V2, start4 = Zea_mays.AGPv34_anno$start4, end4 = Zea_mays.AGPv34_anno$end4, name = Zea_mays.AGPv34_anno$name)

snps_p7msnps <- fread("./largedata/geno/hmp321_282_agpv4_maf005_miss03_pruned.bim", header = FALSE , data.table=FALSE)

####################################################################################################################################
#GWAS_TRAIT_FILE_PROCESSING


pval <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/gwas/output/TasselPrimaryBranches.assoc.txt", data.table=FALSE)

pval$p <- -log10(pval$p_wald)

res <- subset(pval, p > 4)
fwrite(res, "./largedata/circos/TasselPrimaryBranches.assoc_pval4.txt", sep=",", row.names=FALSE, quote=FALSE)
inputdf <- fread("./largedata/circos/TasselPrimaryBranches.assoc_pval4.txt", data.table=FALSE)                                                  
gwas <- data.frame(Chr = qq("Chr{inputdf$chr}"), Start = inputdf$ps, End = inputdf$ps, p = inputdf$p)




####################################################################################################################################

####################################################################################################################################


#GWAS_VISUAL


tiff("./largedata/circos/fixed_TasselPrimaryBranches_p7msnps_L3Base.tiff",res=600,units = "mm",height = 120,width = 120)
par(mar=c(0,0,1,0))
circos.genomicInitialize(data,plotType="NULL")
#chr_col=colours()[c(12,41,46,52,60,79,125,414,429,190)]
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, border = NA,col = colours()[407])
  circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = "black",font = 2,
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.1, bg.border = NA)

title(main = "TasselPrimaryBranches_L3Base" , cex.main = 1)

bg.col <- rep(colours()[c(407,140)], 5)

ch=gwas[,1]
ch=gsub("Chr","",ch)
col=ifelse(as.numeric(ch)%%2==1,"darkblue", "darkred")
circos.genomicTrackPlotRegion(gwas,panel.fun = function(region, value, ...){
  #circos.genomicLines(region, value, type = "h")
  circos.genomicPoints(region, value, pch = 16, cex = 0.25, ...)}
  ,bg.col =bg.col, bg.border = "white",track.height = 0.40
) 

res_fixed_TasselPrimaryBranches_p7msnps_L3Base <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/res.fixed.pcut_sum_TasselPrimaryBranches_ts_L3Base.csv", header = T , data.table=FALSE)

if (res_fixed_TasselPrimaryBranches_p7msnps_L3Base$n.direct >= 1) {
  
  #Direct SNPs
  dsnps_fixed_TasselPrimaryBranches_p7msnps_L3Base <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/res.fixed_dsnps_TasselPrimaryBranches_ts_L3Base.csv", header = T , data.table=FALSE)
  
  hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps <- filter(snps_p7msnps, V2 %in% dsnps_fixed_TasselPrimaryBranches_p7msnps_L3Base$snp)
  hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps <- select(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps, V1, V4, V2)
  colnames(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps)<-c("chr", "pos","name")
  
  
  fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps <- data.frame(Chr = qq("Chr{hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps$chr}"), Start = hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps$pos , End = hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps$pos, value = rep(0.4, nrow(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps)))
  fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps$value <- c(rep(10.5, nrow(fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps)))
  fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps <- rbind(fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps[1,], fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps)
  fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps$value[1] = 1
  circos.genomicTrack(fixed_TasselPrimaryBranches_p7msnps_L3Base_dsnps, 
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, type = "h", col="red")
                      }, track.index =2)
  
  
  circos.genomicTrack(p_lines, stack = TRUE, 
                      panel.fun = function(region, value, ...) {
                        i = getI(...)
                        circos.genomicLines(region, value, col = "#FFFF00", ...)
                      }, track.index =2)
  
  
  
}



####################################################################################################################################
#MED_VISUAL

if (res_fixed_TasselPrimaryBranches_p7msnps_L3Base$n.med >= 1) {
  
  fixed_TasselPrimaryBranches_p7msnps_L3Base <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/res.fixed_mediator_TasselPrimaryBranches_ts_L3Base.csv", header = T , data.table=FALSE)
  
  
  hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34 <- filter(Zea_mays.AGPv34_anno, gene3 %in% fixed_TasselPrimaryBranches_p7msnps_L3Base$id) 
  
  if (length(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$gene4) >= 1) {
    
    color = c("#CC0033", "#FF9900", "#CCCC00", "#00FF00", "#3399FF", "#00FFFF", "#FF00FF", "#990066", "#999999", "#000000")
    hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$col = color[1:nrow(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34)]   
    
    fwrite(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34, "./largedata/circos/hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34.csv", sep=",", row.names=FALSE, quote=FALSE)
    
    fixed_TasselPrimaryBranches_p7msnps_L3Base_v4 <- data.frame(chr=hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$chr, start=hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$start4, end = hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$end4, value = rnorm(nrow(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34), 0, 0.5) )
    fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$chr =qq("Chr{fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$chr}")
    
    
    
    fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$col = hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34$col
    
    circos.track(ylim = c(0, 0.05), track.height = 0.05, bg.border = "black")
    
    for(i in 1:nrow(fixed_TasselPrimaryBranches_p7msnps_L3Base_v4))
    {
      circos.rect((fixed_TasselPrimaryBranches_p7msnps_L3Base_v4[i,2]-1000000),0,(fixed_TasselPrimaryBranches_p7msnps_L3Base_v4[i,3]+1000000),0.05,sector.index=fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$chr[i], col= fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$col[i], border=fixed_TasselPrimaryBranches_p7msnps_L3Base_v4$col[i],track.index =3)
    }
    
    
    
    
    isnps_fixed_TasselPrimaryBranches_p7msnps_L3Base <- fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/med_output/res.fixed_indirect_TasselPrimaryBranches_ts_L3Base.csv")
    isnps_fixed_TasselPrimaryBranches_p7msnps_L3Base = isnps_fixed_TasselPrimaryBranches_p7msnps_L3Base[,-3]
    genev4_isnps <- merge(hdf_fixed_TasselPrimaryBranches_p7msnps_L3Base_v34, isnps_fixed_TasselPrimaryBranches_p7msnps_L3Base, by.x = "gene3", by.y = "medi")
    
    genev4_isnps_pos <- merge(genev4_isnps, snps_p7msnps, by.x = "snps_for_medi", by.y = "V2")
    
    genev4_isnps_pos <- arrange(genev4_isnps_pos, chr, start4 )
    
    
    
    genev4_isnps_pos$chr <- qq("Chr{genev4_isnps_pos$chr}")
    genev4_isnps_pos$V1 <- qq("Chr{genev4_isnps_pos$V1}")
    
    for (i in 1 : nrow(genev4_isnps_pos)) {
      circos.link(genev4_isnps_pos$chr[i], c(genev4_isnps_pos$start4[i], genev4_isnps_pos$end4[i]), genev4_isnps_pos$V1[i], genev4_isnps_pos$V4[i], col = genev4_isnps_pos$col[i],  border = genev4_isnps_pos$col[i])
      
    }
    
    
  }
  
  
  
  
  
}


dev.off()











