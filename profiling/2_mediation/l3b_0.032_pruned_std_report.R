setwd('/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020')
library(data.table)
rna <- fread("/common/jyanglab/shared/dbcenter/Kremling_Nature3RNASeq282_March2018/Expression_matrix/df_STAR_HTSeq_counts_B73_match_based_on_genet_dist_DESeq2_normed_fpm_rounded_origNames_and_Altnames.txt", data.table=FALSE)
l3b <- subset(rna, Tissue %in% c("L3Base"))
dim(l3b) # 295 37136
# remove the duplicated ones
l3b0 <- l3b[!duplicated(l3b$HMP32Name), ]
dim(l3b0) #278 37136

#determine the first column of the gene expression
idx1 <- which(names(l3b0) == "AC148152.3_FG001")

X <- l3b0[, c(6,10:ncol(l3b0))]
X$HMP32Name <- gsub("282set_", "", X$HMP32Name)


variances<-apply(X[,2:ncol(X)], 2, var)
X <- X[, -(which(variances<=0.05)+1)] # 278 25038
means<-apply(X[,2:ncol(X)], 2, mean)
X <- X[, -(which(means<=1)+1)] #278 19357


p <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/geno_282/allchr_bisnp_n282_snpid_maf01_geno2_0.032_pruned_NA_0_matrix.txt", data.table=FALSE)
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



fwrite(as.list(colnames(X)), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/colnamesX_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)
fwrite(as.list(colnames(Z)), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/colnamesZ_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)



library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)
source('highmed2019.r')
source('fromSKAT.R')


###' New wrapper functions for each method, they source the function script highmed2019.r

adlassoMixWrapper <- function(Z,X,X0,y,kernel='linear',ncores=1,pmax=length(y)-2){  #This function is a wrapper function that runs the adaptive lasso for linear mixed model for a sequence of penalty parameter lambda, and output the result that minimizes BIC
  if(kernel=='shrink_EJ'){
    K = A.mat(Z,shrink=list(method="EJ"))
  }
  if(kernel=='linear'){
    K = A.mat(Z,shrink=FALSE)
  }
  n = length(y)
  eigK = eigen(K+sqrt(n)*diag(n))
  Qeig = eigK$vectors
  thetaeig = eigK$values-sqrt(n)
  Xt = t(Qeig)%*%cbind(cbind(1,X0),X)
  yt = t(Qeig)%*%y
  p = ncol(X)
  lambda.seq = getLambda(Xt,yt,nlambda=100,intercept=F,penalty.factor = c(rep(0,1+ncol(X0)),rep(1,p)))
  #    getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  results.all = foreach(lambda=lambda.seq) %dopar% {
    library(glmnet)
    library(MASS)
    source('highmed2019.r')
    try(adlassoMix(X,y,K,X0=X0,lambda=lambda,init=list(tau=var(y)/2,ve=var(y)/2),pmax = pmax,err.max=1e-5,iter.max=200,iter.lasso.max=1e4,method.var='REML')) ## 'MLE') ## method.var=MLE or REML
  }
  names(results.all) = lambda.seq
  negll.all = sapply(results.all, getElement,name='negll')
  s0.all = sapply(results.all, getElement,name='s0')
  lambda.all = sapply(results.all, getElement,name='lambda')
  bic.all = sapply(results.all, getElement,name='bic')
  mod.final = results.all[[which.min(bic.all)]]
  return(list(model=mod.final,bic=bic.all,negll=negll.all,s0=s0.all,lambda=lambda.all,results.all=results.all)) 
}

###' 
###' 
adlasso2typeWrapper <- function(y,X0,X,Z,pz=seq(0.01,0.99,length.out=20),pmax=length(y)-2,ncores=16){ ## wrapper of medfix with tuning parameter selection based on BIC
  ## the fixed model should NOT include the intercept in X0
  n = length(y)
  p0 = ncol(X0)
  p = ncol(X)
  q = ncol(Z)
  getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  out.all = foreach(ppzz=pz) %dopar% {
    library(glmnet)
    library(MASS)
    source('highmed2019.r')
    adlasso.2type(y,X,Z,ppzz,X0,pmax=pmax)
  }
  stopCluster(cl)
  bic = sapply(out.all, getElement,name='bic.min')
  id.opt = which.min(bic)[1]
  out.final = out.all[[id.opt]]
  lambda.opt = out.final[['lambda.seq']][which.min(out.final[['bic']])[1]]
  mod.final = out.final[['model']]
  mod.final[['X0']] = X0
  mod.final[['X']] = X
  mod.final[['Z']] = Z
  mod.final[['y']] = y
  mod.eq= adlasso.2type(y,X,Z,pz=0.5,X0,pmax=pmax)[['model']]
  mod.eq[['X0']] = X0
  mod.eq[['X']] = X
  mod.eq[['Z']] = Z
  mod.eq[['y']] = y    
  return(list(model=mod.final,mod.eq=mod.eq,pz.opt=pz[id.opt],lambda.opt=lambda.opt,bic.prop=bic,bic.lambda=out.final[['bic']],pz.seq=pz,lambda.seq=out.final[['lambda.seq']]))
}

e2mFixed <- function(mod,ncores=8){ ## calculate the exposure to mediator effect using the fixed model
  ### mod = mod.fixed
  X = mod[['X']]
  n = nrow(X)
  X0 = mod[['X0']]#[,-1]
  if(is.null(X0)){
    p0=1
  }else{
    p0=ncol(X0)+1
  }
  Z = mod[['Z']]
  b = as.matrix(mod[['b']])
  ixnz = (b!=0)
  ns = sum(ixnz)
  #    print(ns)
  if(ns>0){
    bnz = b[ixnz]
    Xnz = matrix(X[,ixnz],nrow=n)
    #    print('good so far')
    getDoParWorkers()
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    Ball = foreach(ii=1:ns) %dopar% {
      library(glmnet)
      library(MASS)
      source('highmed2019.r')
      ## ii = 1
      vs.lasso = adLasso(X=Z,y=Xnz[,ii],X0=X0)
      all.lasso = bic.glmnet(vs.lasso,cbind(1,cbind(X0,Z)),Xnz[,ii],p0=p0)
      mod.lasso = all.lasso$model
      list(coef=mod.lasso$b,cov=mod.lasso$cov.coef,dv=mod.lasso$dev,s0=mod.lasso$s0)#,bic=all.lasso$bic,df=all.lasso$df,lambda=all.lasso$lambda.seq)
    }
    stopCluster(cl)
  }else{
    Ball = NULL
  }
  return(Ball)
}


ncores = 16

###' run the estimation step for MedFix
vs.fixed = adlasso2typeWrapper(y,X0,X,Z,ncores=ncores) 
mod.fixed = vs.fixed[['model']]  # extract the model that minimizes BIC
mod.eq = vs.fixed[['mod.eq']] # extract the model that assign equal penalty on the two data types.
e2m.fixed= e2mFixed(mod.fixed,ncores=ncores) # for mod.fixed, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
e2m.eq= e2mFixed(mod.eq,ncores=ncores) # for mod.eq, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
###' run the estimation step for MedMix using linear kernel, and extract the model that minimizes BIC 
vs.mix.linear = adlassoMixWrapper(Z,X,X0,y,kernel='linear',ncores=ncores)
mod.mix.linear = vs.mix.linear[['model']] 
###' run the estimation step for MedMix using shrink_EJ kernel, and extract the model that minimizes BIC 
vs.mix.shrink = adlassoMixWrapper(Z,X,X0,y,kernel='shrink_EJ',ncores=ncores)
mod.mix.shrink = vs.mix.shrink[['model']]

###' perform the hypotheses testing step that control fasle discovery proportion (FDP) below gamma
###' P(FDP>gamma)<alpha where gamma = 0.1, and alpha could be 0.05
p.adj.method = 'holmr' # I wrote an implementation of Holm test that controls FDP
pval.cut=0.05
tests.fixed=testMedFix(mod.fixed,e2m.fixed,p.adj.method=p.adj.method)
(res.fixed.pcut = medH.L2fixed(mod.fixed,e2m.fixed,tests.fixed,pval.cut=pval.cut))# report the results including the PVM calculation for MedFix_BIC
fwrite(as.list(res.fixed.pcut), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/res.fixed.pcut_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)

tests.eq=testMedFix(mod.eq,e2m.eq,p.adj.method=p.adj.method)
(res.eq.pcut = medH.L2fixed(mod.eq,e2m.eq,tests.eq,pval.cut=pval.cut))# report the results including the PVM calculation for MedFix_{0.5}
fwrite(as.list(res.eq.pcut), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/res.eq.pcut_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)

tests.mix.linear = testMedMix(mod.mix.linear,p.adj.method=p.adj.method)
(res.mix.linear.pcut = medH.L2(mod.mix.linear,tests.mix.linear,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_linear
fwrite(as.list(res.mix.linear.pcut), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/res.mix.linear.pcut_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)

tests.mix.shrink = testMedMix(mod.mix.shrink,p.adj.method=p.adj.method)
(res.mix.shrink.pcut = medH.L2(mod.mix.shrink,tests.mix.shrink,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_shrink
fwrite(as.list(res.mix.shrink.pcut), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/res.mix.shrink.pcut_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)


###' a function extracting the mediators making the cut

reportMediator <- function(mod,tst,pval.cut){
  coef.outcome = mod[['cov.coef']][['b.nz']]
  med = tst[['med.individual']]
  ix.output = (med$padj<=pval.cut)
  output = list(pval=med[ix.output,],coef.outcome=coef.outcome[ix.output])
  return(output)
}

reportDirectSNPs <- function(mod,tst,pval.cut){
  coef.outcome = mod[['cov.coef']][['u.nz']]
  p.direct = tst$direct
  output = list(pval=p.direct,coef.outcome=coef.outcome)
  return(output)
}

colnames_X<-colnames(X)
reportINDirectSNPs <- function(mod,e2m,tst,pval.cut, Z, colnames_X){
  med = tst[['med.individual']]
  ix.output = (med$padj<=pval.cut)
  idx <- c(1:length(ix.output))[ix.output]
  e2m.sig <- e2m[idx]
  medid <- med$id[ix.output]
  colnames_X <- colnames_X
  mednames <- colnames_X[medid]
  Znames <- colnames(Z)
  #
  output <- data.frame()
  for(i in 1:length(med$id[ix.output])){
    idxnsnps_i <- (e2m.sig[[i]]$coef != 0)
    snpsnames_i <- Znames[idxnsnps_i]
    med_snps <- data.frame(medi=rep(mednames[i], length(snpsnames_i)), snps_for_medi=snpsnames_i)
    output <- rbind(output, med_snps)
  }
  #
  
  return(output)
}


Xcolnames <- fread("/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/colnamesX_l3b.csv", data.table=FALSE)

###' Report the mediators selected by each method
mediators.fixed = reportMediator(mod.fixed,tests.fixed,pval.cut)
id.fixed = mediators.fixed[['pval']]$id
if(length(id.fixed) >= 1){fwrite(Xcolnames[id.fixed], "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.fixed_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)}
directsnps.fixed = reportDirectSNPs(mod.fixed,tests.fixed,pval.cut)
if(length(names(directsnps.fixed$coef.outcome))>=1){fwrite(as.list(names(directsnps.fixed$coef.outcome)), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.fixed_gddsilk_0.032_pruned_l3b_dsnps.csv", sep=",", row.names=FALSE, quote=FALSE)}
indirect_snps.fixed <- reportINDirectSNPs(mod.fixed, e2m.fixed, tests.fixed,pval.cut,Z,colnames_X)
if(nrow(indirect_snps.fixed)>=1){fwrite(indirect_snps.fixed, "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.fixed_gddsilk_0.032_pruned_l3b_indirect_snps.csv", sep=",", row.names=FALSE, quote=FALSE)}




mediators.eq = reportMediator(mod.eq,tests.eq,pval.cut)
id.eq = mediators.eq[['pval']]$id
if(length(id.eq) >= 1){fwrite(Xcolnames[id.eq], "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.eq_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)}
directsnps.eq = reportDirectSNPs(mod.eq,tests.eq,pval.cut)
if(length(names(directsnps.eq$coef.outcome))>=1){fwrite(as.list(names(directsnps.eq$coef.outcome)), "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.eq_gddsilk_0.032_pruned_l3b_dsnps.csv", sep=",", row.names=FALSE, quote=FALSE)}
indirect_snps.eq <- reportINDirectSNPs(mod.eq, e2m.eq, tests.eq,pval.cut,Z,colnames_X)
if(nrow(indirect_snps.eq)>=1){fwrite(indirect_snps.eq, "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.eq_gddsilk_0.032_pruned_l3b_indirect_snps.csv", sep=",", row.names=FALSE, quote=FALSE)}





mediators.mix.linear = reportMediator(mod.mix.linear,tests.mix.linear,pval.cut)
id.mix.linear = mediators.mix.linear[['pval']]$id
if(length(id.mix.linear) >= 1){fwrite(Xcolnames[id.mix.linear], "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.mix.linear_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)} 

mediators.mix.shrink = reportMediator(mod.mix.shrink,tests.mix.shrink,pval.cut)
id.mix.shrink = mediators.mix.shrink[['pval']]$id
if(length(id.mix.shrink) >= 1){fwrite(Xcolnames[id.mix.shrink], "/common/jyanglab/zhikaiyang/projects/mediation-analysis/largedata/code2020/id.mix.shrink_gddsilk_0.032_pruned_l3b.csv", sep=",", row.names=FALSE, quote=FALSE)}
