setwd('/common/jyanglab/zhikaiyang/projects/mediation')

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)



nQTLperM <- as.numeric(as.character(args[1]))
nQTLar <- as.numeric(as.character(args[2]))
h2_as <- as.numeric(as.character(args[3]))
seed <- as.numeric(as.character(args[4]))



library(data.table)



trait = "sim"

library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)
source('lib/highmed2019.r')
source('lib/fromSKAT.R')


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
    source('lib/highmed2019.r')
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
    source('lib/highmed2019.r')
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
      source('lib/highmed2019.r')
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


ncores = 4



###' a function extracting the mediators making the cut
reportDirectSNPs <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['u.nz']]
  p.direct = tst$direct
  output <- data.frame()
  if(!is.null(coef.outcome)){
    output = data.frame(snp=names(coef.outcome), pval=p.direct, coef=coef.outcome)
  }
  return(output)
}

reportMediator <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['b.nz']]
  med <- data.frame()
  
  if(!is.null(coef.outcome)){
    med = tst[['med.individual']]
    # ix.output = (med$padj <= pval.cut)
    # output = list(pval=med[ix.output,],coef.outcome=coef.outcome[ix.output])
    med$coef <- coef.outcome
    med$id <- row.names(med)
  }
  return(med)
}

reportINDirectSNPs <- function(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X)){
  med = tst[['med.individual']]
  ix.output = med$padj <= pval.cut
  idx <- c(1:length(ix.output))[ix.output]
  
  output <- data.frame()
  if(length(idx)>=1){
    e2m.sig <- e2m[idx]
    medid <- med$id[ix.output]
    colnames_X <- colnames_X
    mednames <- colnames_X[medid]
    Znames <- colnames(Z)
    for(i in 1:length(med$id[ix.output])){
      idx <- which(e2m.sig[[i]]$coef != 0)
      snpsnames <- Znames[idx]
      med_snps <- data.frame(medi=mednames[i], snps_for_medi=snpsnames, coef=e2m.sig[[i]]$coef[idx])
      output <- rbind(output, med_snps)
    }
  }
  return(output)
}



mediation_run = function(nQTLperM=10, nQTLar=10, h2_as=0.75, seed=1231){
  
  for (topi in 1:12) {
    
    
    a = fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/gddtosilk_22ksnps_matrix_new.txt",header = T, data.table = F)
    geno = cbind(1:nrow(a), a[,3:ncol(a)])
    
    
    set.seed(seed)
    id_sample = sort(sample(2:ncol(geno),20000))
    geno = geno[,c(1,id_sample)]
    colnames(geno)[1] = "id"
    geno = as.matrix(geno)
    geno[,2:ncol(geno)]=geno[,2:ncol(geno)]+1
    dna = geno
    
    nloci=20000
    inter=4
    pos=which(1:nloci %% inter ==0)
    
    
    
    pathname = paste("largedata/simulation/mediationh2/M_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    rna <- fread(pathname, header=F,data.table=FALSE)
    
    pathname = paste("largedata/simulation/mediationh2/qtlar_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    qtlar = fread(pathname, header=F,data.table=FALSE)
    
    pathname = paste("largedata/simulation/mediationh2/qtl_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    qtl = fread(pathname, header=F,data.table=FALSE)
    
    qtl = as.vector(as.matrix(qtl))
    qtlar = as.vector(as.matrix(qtlar))
    qtlid = unique(c(qtlar, qtl))
    qtlid = sort(qtlid)
    
    print(topi)
    dna = dna[,c(1, c(pos[qtlid]+1))]
    
    pcs = prcomp(as.matrix(dna[,sample(2:ncol(dna),1000)]))
    pcs_m = pcs$x[,1:3]
    
    pathname = paste("largedata/simulation/mediationh2/y_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    y = fread(pathname, header=F, data.table=FALSE)
    
    
    y = as.matrix(y)
    Z = as.matrix(dna[,-1])
    X = as.matrix(rna)
    X0 = as.matrix(pcs_m)
    ### dim(X)  #244 12993
    
    ###remove genes with expression correlation with phenotype lower than 0.05
    cors<-apply(X, 2, function(X) cor(X, y))
    X  <- X[, -(which((cors<=0.05) & (cors>=-0.05) ))] #255 12570
    
    
    pathname = paste("largedata/simulation/mediationh2/resulth2/colnamesX_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    fwrite(as.list(colnames(X)), pathname, sep=",", row.names=FALSE, quote=FALSE)
    pathname = paste("largedata/simulation/mediationh2/resulth2/colnamesZ_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    fwrite(as.list(colnames(Z)), pathname, sep=",", row.names=FALSE, quote=FALSE)
    
    Z=Z-1
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
    
    pathname = paste("largedata/simulation/mediationh2/resulth2/colnamesX_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    Xcolnames <- fread(pathname, data.table=FALSE)
    
    
    tests.mix.linear = testMedMix(mod.mix.linear,p.adj.method=p.adj.method)
    (res.mix.linear.pcut = medH.L2(mod.mix.linear,tests.mix.linear,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_linear
    pathname = paste("largedata/simulation/mediationh2/resulth2/res_mix_linear_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    fwrite(as.list(res.mix.linear.pcut), pathname, sep=",", row.names=FALSE, quote=FALSE)
    
    tests.mix.shrink = testMedMix(mod.mix.shrink,p.adj.method=p.adj.method)
    (res.mix.shrink.pcut = medH.L2(mod.mix.shrink,tests.mix.shrink,pval.cut=pval.cut))# report the results including the PVM calculation for MedMix_shrink
    pathname = paste("largedata/simulation/mediationh2/resulth2/res_mix_shrink_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    fwrite(as.list(res.mix.shrink.pcut), pathname, sep=",", row.names=FALSE, quote=FALSE)
    
    
    
    ###' Report the mediators selected by each method
    
    #############
    
    
    
    
    
    mediators.mix.linear = reportMediator(mod=mod.mix.linear, tst=tests.mix.linear)
    pathname = paste("largedata/simulation/mediationh2/resulth2/mediators_mix_linear_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    if(nrow(mediators.mix.linear) >= 1){fwrite(mediators.mix.linear, pathname, sep=",", row.names=FALSE, quote=FALSE)}
    
    
    
    mediators.mix.shrink = reportMediator(mod=mod.mix.shrink, tst=tests.mix.shrink)
    pathname = paste("largedata/simulation/mediationh2/resulth2/mediators_mix_shrink_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".csv",sep = "")
    if(nrow(mediators.mix.shrink) >= 1){fwrite(mediators.mix.shrink, pathname, sep=",", row.names=FALSE, quote=FALSE)}
    
    
    
  }
  
  
  
}


mediation_run(nQTLperM, nQTLar, h2_as, seed)

