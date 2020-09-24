slurm_med_wrapper <- function(ti=2, cutoff_pm=0.05, ncores=12,
                              phenofile = "data/geno_trait.txt",
                              genofile = "largedata/geno/allchr_bisnp_n282_snpid_maf01_geno2_pruned_NA_0_matrix.txt",
                              rnafile = "largedata/RNA-seq/filtered_GRoot.csv",
                              pcfile="largedata/allchr_bisnp_n282_snpid_maf01_geno2_pruned.eigenvec",
                              outdir="largedata/med_output"){
  # ti: the ith trait. [int, =2]
  # cutoff_pm: remove genes with expression correlation with phenotype lower than this cutoff. [num, = 0.05]
  
  
  ### phentoypic data
  pheno <- read.delim(phenofile, header=T)
  ### working on the ith phenotype
  ph <- pheno[, c(1, ti)]
  ph <- ph[!is.na(ph[,2]),]
  
  rna <- fread(rnafile, data.table=FALSE)
  
  geno <- fread(genofile, data.table=FALSE)
  names(geno)[1] <- "genotype"
  # 275 696549
  
  pc <- fread(pcfile, data.table=FALSE)
  names(pc)[2] <- "genotype"
  
  ### find the common set of genotypes by merging them togather
  df <- Reduce(function(x, y) merge(x, y, by="genotype"), list(ph, rna[, 1:2], geno[, 1:2], pc[,1:2]))
  
  ### stat for output
  sout <- data.frame(trait=names(pheno)[ti], rna=rnafile, ngeno=nrow(df))
  
  df1 <- merge(df[, 1:2], rna, by="genotype")
  # must be 0
  if(sum(df$genotype != df1$genotype) != 0){
    stop("Error! RNA-seq genotype id not matching df")
  } 
  
  ### remove genes with expression correlation with phenotype lower than 0.05
  y <- df1[,2]
  X <- df1[, -1:-2]
  cors <- apply(X, 2, function(X) cor(X, y))
  X  <- X[, -(which((cors <= cutoff_pm) & (cors >= -cutoff_pm) ))] #255 12570
  
  X <- as.data.frame(X)
  X <- mutate_all(X, scale)
  #row.names(X) <- gid
  X <- as.matrix(X)
  sout$nM <- ncol(X) #num of M (mediators) after normalization!
  
  ### first 3 PC were used for the analysis
  X0 <- merge(df[,1:2], pc[, 2:5], by="genotype")
  # must be 0
  if(sum(df$genotype != X0$genotype) != 0){
    stop("Error! PC genotype id not matching df")
  } 
  X0 <- as.matrix(X0[, -1:-2])
  sout$nPC <- ncol(X0)
  
  Z <- merge(df[, 1:2], geno[, -2], by="genotype")
  # must be 0
  if(sum(df$genotype != Z$genotype) != 0){
    stop("Error! SNP genotype id not matching df")
  } 
  Z <- as.matrix(Z[, -1:-2])
  sout$nZ <- ncol(Z)
  
  ###' run the estimation step for MedFix
  vs.fixed = adlasso2typeWrapper(y, X0, X, Z, ncores=ncores) 
  # extract the model that minimizes BIC
  mod.fixed = vs.fixed[['model']]  
  # extract the model that assign equal penalty on the two data types.
  mod.eq = vs.fixed[['mod.eq']] 
  # for mod.fixed, run the mediator models for the mediators with non-zero estimated
  # coefficient in the outcome model
  e2m.fixed= e2mFixed(mod.fixed, ncores=ncores) 
  # for mod.eq, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
  e2m.eq= e2mFixed(mod.eq, ncores=ncores) 
  
  ###' perform the hypotheses testing step that control fasle discovery proportion (FDP) below gamma
  ###' P(FDP>gamma)<alpha where gamma = 0.1, and alpha could be 0.05
  p.adj.method = 'holmr' 
  # I wrote an implementation of Holm test that controls FDP
  pval.cut=0.05
  
  ### Report results for fixed (BIC) ------------------------------------------####
  tests.fixed=testMedFix(mod.fixed, e2m.fixed, p.adj.method=p.adj.method)
  # report the results including the PVM calculation for MedFix_BIC
  res.fixed.pcut = medH.L2fixed(mod.fixed, e2m.fixed, tests.fixed, pval.cut=pval.cut)
  outfile <- paste0(outdir, "/res.fixed.pcut_sum_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
  fwrite(as.list(res.fixed.pcut), outfile, sep=",", row.names=FALSE, quote=FALSE)
  
  
  directsnps.fixed = reportDirectSNPs(mod=mod.fixed, tst=tests.fixed)
  if(nrow(directsnps.fixed)>=1){
    outfile <- paste0(outdir, "/res.fixed_dsnps_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(directsnps.fixed, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
  ###' Report the mediators selected by each method
  mediators.fixed = reportMediator(mod=mod.fixed, tst=tests.fixed)
  # id.fixed = mediators.fixed[['pval']]$id
  if(nrow(mediators.fixed) >= 1){
    outfile <- paste0(outdir, "/res.fixed_mediator_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(mediators.fixed, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
  indirect_snps.fixed <- reportINDirectSNPs(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X))
  if(nrow(indirect_snps.fixed)>=1){
    outfile <- paste0(outdir, "/res.fixed_indirect_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(indirect_snps.fixed, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
  
  ### Report results for eq (i.e., 0.5) ------------------------------------------####
  tests.eq=testMedFix(mod.eq,e2m.eq,p.adj.method=p.adj.method)
  # report the results including the PVM calculation for MedFix_{0.5}
  res.eq.pcut = medH.L2fixed(mod.eq, e2m.eq, tests.eq, pval.cut=pval.cut)
  outfile <- paste0(outdir, "/res.eq_sum_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
  # outfile <- paste0("largedata/code2020/res.eq.pcut_gddsilk_0.032_pruned_ld268.csv")
  fwrite(as.list(res.eq.pcut), outfile, sep=",", row.names=FALSE, quote=FALSE)
  
  directsnps.eq = reportDirectSNPs(mod=mod.eq, tst=tests.eq)
  if(nrow(directsnps.eq)>=1){
    outfile <- paste0(outdir, "/res.eq_dsnps_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(directsnps.eq, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
  ###' Report the mediators selected by each method
  mediators.eq = reportMediator(mod=mod.eq, tst=tests.eq)
  if(nrow(mediators.eq) >= 1){
    outfile <- paste0(outdir, "/res.eq_mediator_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(mediators.eq, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
  indirect_snps.eq <- reportINDirectSNPs(mod=mod.eq, e2m=e2m.eq, tst=tests.eq, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X))
  if(nrow(indirect_snps.eq)>=1){
    outfile <- paste0(outdir, "/res.eq_indirect_", names(pheno)[ti], "_ts_", gsub(".*_", "", rnafile))
    fwrite(indirect_snps.eq, outfile, sep=",", row.names=FALSE, quote=FALSE)}
  
}