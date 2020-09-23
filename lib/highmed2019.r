require(glmnet)
require(MASS)
require(rrBLUP)
require(parallel)
require(doParallel)
require(CompQuadForm)
### Utility functions
### COPYRIGHT: QI ZHANG

myginv <- function (X, tol = sqrt(.Machine$double.eps)){
            if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
                                  stop("'X' must be a numeric or complex matrix")
                        if (!is.matrix(X))
                                                  X <- as.matrix(X)
                        Xeigen <- eigen(X,symmetric=T)
                        Positive <- Xeigen$values > max(tol * Xeigen$values[1L], 0)
                        if (all(Positive))
                                                  Xeigen$vectors %*% (1/Xeigen$values * t(Xeigen$vectors))
                        else if (!any(Positive))
                                                  array(0, dim(X)[2L:1L])
                        else Xeigen$vectors[, Positive, drop = FALSE] %*% ((1/Xeigen$values[Positive]) *
                                                                                                                         t(Xeigen$vectors[, Positive, drop = FALSE]))
   }

standardize <- function(x){
  out = (x-mean(x,na.rm=T))/sd(x,na.rm=T)
  out[is.na(out)]=0
  out
}

bonfg <- function(p,n=length(p)){
    p[is.na(p)] = 1
    out = p*n
    return(out)
}

holm.k <- function(p,k,n=length(p)){
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)
    i <- seq_len(lp)
    i[1:k] = k
    o <- order(p)
    ro <- order(o)
    p0nna = pmin(1, cummax((n + k - i) * p[o]/k))[ro]
    p0[nna]=p0nna
    return(p0)
}

holm.r <- function(p,r,n=length(p)){
    nm <- names(p)
    p <- as.numeric(p)
    p0 <- setNames(p, nm)
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)
    i <- seq_len(lp)
    ri = floor(i*r)
    o <- order(p)
    ro <- order(o)
    cj = function(j) sum(1/(1:j))
    cri = sapply(ri+1, cj)
    p0nna = pmin(1, cummax(cri*(n + ri+1 - i) * p[o]/(ri+1)))[ro]
    p0[nna]=p0nna
    return(p0)
}


myPadjust <- function(p,method=p.adjust.methods,k=2,r=0.1,n=length(p)){
    if(method=='holmk'){
        output = holm.k(p,k=k,n=n)
        }
    else{
        if(method=='holmr')
          output = holm.r(p,r=r,n=n)
        else{
            if(method=='bonfg')
             output = bonfg(p,n=n)
        else
             output = p.adjust(p,method=p.adjust.methods,n=n)
          }
        }
    return(output)
}



###' functions for medmix
 
getLambda <- function(X,y,nlambda=50,intercept=T,penalty.factor){ ### determine the candidate sequence of the regularization parameter lambda
  ### first run plain lasso by default to get a sequence of prediction errors and a sequence of lambda's   
  mod = glmnet(X,y,intercept=intercept,penalty.factor=penalty.factor,standardize=F)  
  beta = mod$beta
  a0 = mod$a0
  df = mod$df
  n = length(y)
  p = ncol(X)
  m = length(a0)
  ypred = as.matrix( matrix(a0,n,m,byrow=T)+ X%*%beta)
  res = matrix(y,n,m)-ypred
  vres = colSums(res^2)/(n-pmin(n-1,df))
  ### define the range of the candiate lambda in log-scale for medmix, note that the regularization parameters for medmix is comparable with that for the scaled lasso, so the candidate sequence should be scaled.  
  loglambda.range = log(range(pmax(1e-8,mod$lambda))/sqrt(c(max(vres),min(vres))))
  loglambda = seq(loglambda.range[1]-log(2),loglambda.range[2],length.out = nlambda)
  return(exp(loglambda))
}

negll <- function(y,X,gamma,theta,rho,delta){ #neglog-likelihood for medmix times 1/n
    ypred = X%*%gamma
  negll = mean((rho*y-ypred)^2/(delta*theta+1))/2 - log(rho)+mean(log(theta*delta+1))/2   
  return(negll)
}

  
cov.solve <- function(y,X,gamma,theta,gt, init,bounds = c(1e-09,1e+09)){ ### update the variance components by CCCP
    ypred = X%*%gamma
    ypred2 = ypred^2
    yc = y*ypred
    y2=y^2
    n = length(y)
  f.opt <- function(par){ ## the objective function
      rho = par[1]
      delta=par[2]
      sum((rho^2*y2-2*rho*yc+ypred2)/(delta*theta+1))/2 - n*log(rho)+delta*gt 
  }
  soln <- optim(init,f.opt,lower=bounds[1],upper=bounds[2], method="L-BFGS-B") 
  par.opt <- soln$par
  Ve.opt <- 1/par.opt[1]^2
  tau.opt <- Ve.opt*par.opt[2]
  return(list(ve=Ve.opt,tau=tau.opt))  
}

cov.adlasso <- function(X,b,wadj,lambda,s2,p0){ ## computing the covariance of the adaptive lasso estimates.
  n = nrow(X)
  ix1 = (b!=0)
  b1 = b[ix1]
  X1 = X[,ix1]
  w1 = wadj[ix1]
  xtx = t(X1)%*%X1
#  print(summary(diag(xtx)-2*n*lambda*w1/abs(b1)))
  xtxlamb = myginv(xtx+2*n*lambda*diag(w1/abs(b1),nrow=length(w1)))### glmnet penalized on MSE while adaptive lasso paper penalized on SSE, need to scale with n
  covb = s2*xtxlamb%*%xtx%*%xtxlamb
#  covb = s2*ginv(xtx)
  std.err = sqrt(diag(covb))
  rownames(covb) = names(b1)
  colnames(covb) = names(b1)
  names(std.err) = names(b1)
  if(length(b1)>p0){
  output = list(a=b1[1:p0],b.nz = b1[-(1:p0)], std.err.a = std.err[1:p0],std.err.b = std.err[-(1:p0)],s2=s2,cov.ab=covb,cov.bb=covb[-(1:p0),-(1:p0)])
  }else{
  output = list(a=b1,b.nz = NULL, std.err.a = std.err,std.err.b = NULL,s2=s2,cov.ab=covb,cov.bb=NULL)
  }
  return(output)
}


adlassoMix <- function(X,y,K,X0=NULL,lambda,init=list(tau=0.1,ve=var(y)),pmax=length(y)-2,err.max=1e-6,iter.max=100,iter.lasso.max=10,method.var='REML'){##,p.ridge=0.5){ ### The main function of the variable selection for hd lmm
### tau is the variance associated with the random effect 
  n = length(y)
  p= dim(X)[2]
  ## known covariates, including the intercept
  if(is.null(X0)){ 
      X0 = matrix(1,n,1)
  }else{
      X0 = cbind(1,X0)
  }
  p0 = ncol(X0)
  ### eigen decomposition of K for transforming the predictors and the responses so that the noise term has diagonal covariance matrix.     
  offset <- sqrt(n)
  Hb <- K + offset*diag(n)
  Hb.system <- eigen(Hb, symmetric = TRUE)
  theta <- Hb.system$values - offset
  if (min(theta) < -1e-6) {stop("K not positive semi-definite.")}
  Qk <- Hb.system$vectors##[,1:(n-1)]
  ytld0 = t(Qk)%*%y
  Xtld0 = t(Qk)%*%cbind(X0,X)
#  ## fix the weights for all interations

      #### prepare for the algorithm
  q.all = rep(NA,iter.max)
  ### initial values 
  tau = as.numeric(init[['tau']])
  ve = as.numeric(init[['ve']])  
  err = err.max*100
  q = Inf
  iter = 1
  while(err>err.max&iter<=iter.max){
    ## transformed observations normalized by weights
    rho = 1/sqrt(ve)
    delta = tau/ve
    gt = sum(theta/(theta*delta+1))/2  
    wx = (theta*delta+1)^(-0.5)
    ytld = rho*wx*ytld0
    Xtld = diag(wx)%*%Xtld0
    ## scale the lasso tuning parameter. So we are using scaled adaptive lasso
    ## note that tis is lambda0/sqrt(ve) in our case, because we also scaled the data by 1/sqrt(ve)
    ##  use adaptive lasso. first run lasso to get the variable weights
    glmnet.iter0 = glmnet(Xtld,ytld,family='gaussian',standardize=F,intercept=F,penalty.factor=c(rep(0,p0),rep(1,p)),pmax=pmax,lambda=lambda,maxit=iter.lasso.max) #,weights = wobs)
    ## use the lasso output to define variable weights, the zeros have the same weight as the smallest non-zero coefficient 
    gamma.iter0 = as.numeric(glmnet.iter0$beta)
    wad = c(rep(0,p0),1/pmax(abs(gamma.iter0[-(1:p0)]),sd(ytld)/p))
    wad = wad/mean(wad)
    ## adaptive lasso using the scaled tuning parameter
    glmnet.iter = glmnet(Xtld,ytld,family='gaussian',standardize=F,intercept=F,penalty.factor=wad,pmax=pmax,lambda=lambda,maxit=iter.lasso.max) 
    gamma.iter = matrix(glmnet.iter$beta,ncol=1)
    s0 = sum(gamma.iter[-(1:p0)]!=0)
    if(s0==0) gamma.iter[1:p0] = ginv(t(Xtld[,1:p0])%*%Xtld[,1:p0])%*%(t(Xtld[,1:p0])%*%ytld)  ## glmnet may return empty model if the model selected contains more than pmax variables
    ### now prepare for variance component estimation from the residuals
    vars = cov.solve(ytld0,Xtld0,gamma.iter,theta,gt, init=c(rho,delta),bounds = c(1e-09,1e+09))
    tau = vars[['tau']]
    ve = vars[['ve']]
    beta.iter = sqrt(ve)*gamma.iter
    negll.new = negll(ytld0,Xtld0,gamma.iter,theta,1/sqrt(ve),tau/ve)
    ### update the objective function, i.e., the penalized likelihood with scaled adaptive lasso penalty
    q.new = negll.new + lambda*sum(abs(gamma.iter[-(1:p0)])*wad[-(1:p0)])
    err = abs(q-q.new) 
    q.all[iter]=q.new
    q = q.new 
    iter = iter+1
  } 
  s0 = sum(beta.iter[-(1:p0)]!=0)
  bic=2*negll.new+(p0+s0)*log(n)/n
  f.pred= X0%*%beta.iter[1:p0]+X%*%beta.iter[-(1:p0)]
  u.pred = Qk%*%((1-1/((tau/ve)*theta+1))*(t(Qk)%*%(y-f.pred)))
  cov.coef=cov.adlasso(Xtld,beta.iter,wad,lambda,s2=ve*sum((ytld-Xtld%*%gamma.iter)^2)/(n-p0-s0),p0)### pay attention to this
  output = list(a0=beta.iter[1:p0],b=beta.iter[-(1:p0)],tau = tau,ve = ve,cov.coef=cov.coef,K=K,X=X,X0=X0,y=y,pred = list(f=f.pred,u=u.pred),lambda=lambda,s0=s0,negll=negll.new,bic=bic,converged=(err<err.max),iter=iter-1,q.all = q.all[1:(iter-1)])
  return(output)
}


###' functions for calculating the mediation proportions for MedMix

uMix<- function(y,X,Kmtx){ ## output the genetic prediction values of traits
  ymix = mixed.solve(y,X=X,K=Kmtx,method='REML')
  out = as.numeric(ymix$u)
  names(out) = names(y)
  out
}


medH.L2 <- function(mod,res.test=NA,p.adj.method=NA,pval.cut=1){  ## estimate the indirect and indirect effects of the mixed model. 
  ###mod = mod.final
  X = mod[['X']]
  X0 = mod[['X0']]
  y = mod[['y']]
  K = mod[['K']]
  tk = sum(diag(K))
  u.pred = mod[['pred']][['u']]
  b = mod[['b']]
  n = nrow(X)
  ixnz = (b!=0)
  ns = sum(ixnz)
  if(!is.na(res.test)){
     if(is.na(p.adj.method)){
     ixsig = res.test[['med.individual']]$padj<=pval.cut
     }else{
     ixsig = res.test[['med.individual']][[p.adj.method]]<=pval.cut
#     print(res.test[['med.individual']][[p.adj.method]])
     }
    ns = sum(ixsig)
  }
  if(ns>0){
    bnz = b[ixnz]
    Xnz = X[,ixnz]
    Xnz = matrix(Xnz,nrow=n)
    if(!is.na(res.test)){
        bnz=bnz[ixsig]
        Xnz = Xnz[,ixsig]
    }
    bnz = matrix(bnz,nrow=ns)
    Xnz = matrix(Xnz,nrow=n)
    umed = uMix(Xnz%*%bnz,X=X0,Kmtx=K)
  }else{
    umed=0
  }
  uy = u.pred + umed
  v.tot = sum(uy^2)
  v.med = sum(umed^2)
  v.dir = sum(u.pred^2)
  p.med.pure = v.med/v.tot 
  return(c(pmed.pure = p.med.pure,v.tot=v.tot,v.med=v.med,v.dir=v.dir,n.med=ns,pval.cut=pval.cut))
}


###' Functions for hypotheses testing for MedMix model

scoreTestBase<- function(Q,eig){ ## base function for score test as in SKAT 
  out = Get_Davies_PVal.Lambda(Q,eig)
  return(out$p.value)
}

scoreTestMedmix <- function(mod,X0,X){### score tests for the direct exposure effect on the outcome, and the exposure effects on the individual mediators with non-zero mediator-outcome coefficients.
  y = mod[['y']]
#  X = mod[['X']]
#  X0 = mod[['X0']]
#  Z = mod[['Z']]
  K = mod[['K']]
  med.coef = mod[['b']]
    n = length(y)
  idnz = which(med.coef!=0)
  if(length(idnz)>0){
      Xnz = X[,idnz]
      Xall = cbind(X0,Xnz)
  }else{
      Xall = X0
  }
  pall = ncol(Xall)
  ## compute the projection matrices 
  XtXinvall <- myginv(t(Xall)%*%Xall)
#  print(dim(Xall))
#  print(dim(XtXinvall))
#  print(n)
  Sall <- diag(n) - Xall%*%XtXinvall%*%t(Xall)
  ## estimate the noise variances under the null
  s2y = sum((Sall%*%y)^2)/(n-pall)
  ## evaluate the eigen values of the kernel matrices
  SKSy = Sall%*%K%*%Sall
  eigy = eigen(SKSy,symmetric = T)$values
  eigy= eigy[eigy>=sum(eigy)/1e6]  ## discard the negative and near-zero eigen values
  ## compute the score test statistics 
  Qy= sum(y*(SKSy%*%y))/s2y
  ## evaluate the p-values of the score test usin Davis exact method programed as in SKAT pacakge.
  pvaly = scoreTestBase(Qy,eigy) ###  p-values for the direct exposure effect

  if(length(idnz)>0){
      p0 = ncol(X0)
      XtXinv0 <- myginv(t(X0)%*%X0)
      S0 <- diag(n) - X0%*%XtXinv0%*%t(X0)
      s2med = (colSums((S0%*%Xnz)^2))/(n-p0)
      SKSmed = S0%*%K%*%S0
      eigmed = eigen(SKSmed,symmetric = T)$values
      eigmed= eigmed[eigmed>=sum(eigmed)/1e6]
      Qmed = diag(t(Xnz)%*%(SKSmed%*%Xnz))/s2med
      pvalmed = sapply(Qmed,scoreTestBase,eig=eigmed) ## p-values for the exposure effect on the individual mediators
      med=data.frame(id=idnz,e2m=pvalmed)
  }else{
      med=NULL
  }
  return(list(direct=pvaly,med=med))
}

zTestAdlasso <- function(b,std){
    if(is.null(b)|is.null(std)){
        pval=NULL  ### return NULL if no non-zero coefs.
    }else{
    z = b/std
    pval = 2*(1-pnorm(abs(z)))
    }
    return(pval)
}

####

testMedMix <- function(mod,p.adj.method='holm'){
    score.test = scoreTestMedmix(mod,mod[['X0']],mod[['X']])
    ixnz = (mod[['b']]!=0)
    if(sum(ixnz)!=0){
    pval.adal = zTestAdlasso(mod[['cov.coef']][['b.nz']],mod[['cov.coef']][['std.err.b']])
    pval.med = score.test[['med']]
    pval.med$m2y = pval.adal
    pval.med$e2m2y = pmin(1, pmax(pval.med$m2y,pval.med$e2m)) 
#    pval.med$e2m2y = pval.med$m2y+pval.med$e2m 
    nmed = length(pval.adal)
    pval.med$padj = myPadjust(pval.med$e2m2y, method = p.adj.method) #pmin(1, nmed*pval.med$e2m2y)
#    pval.med$padj = sapply(pval.med$bonf, function(xx) min(1,sum(pval.med$bonf[pval.med$bonf<=xx])))
    output = list(direct=score.test[['direct']],indirect=min(pval.med$e2m2y),total=min(score.test[['direct']],min(pval.med$e2m2y)),med.individual=pval.med) 
    }else{
    output = list(direct=score.test[['direct']],indirect=1,total=score.test[['direct']],med.individual=NULL) 
    }
    return(output)
}


###' functions for MedFix


adLasso <- function(X,y,X0=NULL,pmax=length(y)-2){ ## our implementation of adaptive lasso with gamma = 1 and weights based on regular lasso. ### X0 does not include intercept
    n = length(y)
  if(is.null(X0)){ 
      Xall = X
      p0 = 1
  }else{
    Xall = cbind(X0,X)
    p0= ncol(X0)+1
  }
    p=ncol(X)
    ## first run lasso regression to get the variable weights.
    vs.lasso = glmnet(Xall, y,family='gaussian',intercept=T,standardize=F,penalty.factor=c(rep(0,p0-1),rep(1,p)))
    lambda = vs.lasso$lambda
    #lambda = c(1e100*max(vs.lasso$lambda),vs.lasso$lambda)
    coef0 = as.matrix(rbind(vs.lasso$a0,vs.lasso$beta))
    beta0 = coef0[-(1:p0),]
    #beta0 = cbind(0,beta0)
    nlamb = ncol(beta0)
    winv = pmax(abs(beta0), sd(y)/p)
    wad = rbind(matrix(0,p0,nlamb),1/winv)
    wad = wad%*%diag(1/colMeans(wad))
    adlasso.all = lapply(1:nlamb, function(ii) glmnet(Xall, y,family='gaussian',intercept=T,standardize=F,penalty.factor=wad[-1,ii],lambda=lambda[ii]))
    beta = do.call('cbind',lapply(adlasso.all, getElement,name='beta'))
    dimnames(beta) = list(NULL,NULL)
    df = sapply(adlasso.all,getElement,name='df')
    a0 = sapply(adlasso.all,getElement,name='a0')
    dv = sapply(adlasso.all,deviance)
    beta=rbind(a0,beta)
    print(c(p,p0))
    print(dim(beta))
    ### include the null model deviance(lm(y~NULL))
    if(is.null(X0)){
        mod0 = lm(y~NULL)
    }else{
        mod0 = lm(y~X0)
    }
    coef0 = mod0$coefficients
    dv0 = deviance(mod0)
    beta = cbind(c(coef0,rep(0,p)),beta)
    dv = c(dv0,dv)
    df = c(p0-1,df)
    lambda = c(0,lambda) ## let lambda for the null model be 0 to make the cov calculation easier
    wad = cbind(0,wad)   ## let weights to be 0 to make the cov calc easier.
    return(list(beta=beta,df=df,dv=dv,lambda=lambda,p0=p0,wad=wad))
}

bic.glmnet <- function(mod,X,y,p0){ ## model selection with BIC, X and p0 should include the intercept
    n = length(y)
    p = ncol(X)
    d = ncol(mod$beta)
    res = matrix(y,n,d)-X%*%(mod$beta)
    ve = colSums(as.matrix(res)^2)/(n-mod$df)
    bic = log(ve) + (log(n)/n)*mod$df
    id.opt = which.min(bic)
    beta.opt=mod$beta[,id.opt]
    wad.opt=mod$wad[,id.opt]
#    print('so far so good.')
    cov.coef=cov.adlasso(X=X,b=beta.opt,wadj=wad.opt,lambda=mod$lambda[id.opt],s2=ve[id.opt],p0)
    mod.final = list(a0=beta.opt[1:p0],b=beta.opt[-(1:p0)],s0=sum(beta.opt!=0)-p0,dev=mod$dv[1]-mod$dv[id.opt],lambda=mod$lambda[id.opt],cov.coef=cov.coef,ve=ve[id.opt],wad=wad.opt)
    return(list(model=mod.final,dv=mod$dv,bic=bic,df=mod$df, lambda.seq=mod$lambda))
}


cov.adlasso.2type <- function(X,b,wadj,lambda,s2,p0,p){ ## computing the covariance of the adaptive lasso estimates.
  n = nrow(X)
  ix1 = which(b!=0)
  if(length(ix1)>p0){
#  print(length(ix1))
#  print(sum(ix1))
  b1 = b[ix1]
  X1 = X[,ix1]
  w1 = wadj[ix1]
  xtx = t(X1)%*%X1
#  print(summary(diag(xtx)+2*n*lambda*w1/abs(b1)))
#  print(dim(xtx))
#  print(summary(w1))
#  print(summary(b1))
#  print(dim(diag(w1/abs(b1))))
#  print(lambda)
  xtxlamb = myginv(xtx+2*n*lambda*diag(w1/abs(b1),nrow=length(w1)))### glmnet penalized on MSE while adaptive lasso paper penalized on SSE, need to scale with n
#  print(dim(xtxlamb))
  covb = s2*xtxlamb%*%xtx%*%xtxlamb
#  covb = s2*ginv(xtx)
  std.err = sqrt(diag(covb))
  rownames(covb) = names(b1)
  colnames(covb) = names(b1)
  names(std.err) = names(b1)
  ix10 = (ix1<=p0)
  ix11 = (ix1>p0&ix1<=(p+p0))
  ix12 = (ix1>(p+p0))
  if(sum(ix11)>0){
     b.nz=b1[ix11]
     std.err.b = std.err[ix11]
     cov.bb=covb[ix11,ix11]
  }else{
     b.nz=NULL
     std.err.b = NULL
     cov.bb=NULL
  }
  if(sum(ix12)>0){
     u.nz=b1[ix12]
     std.err.u = std.err[ix12]
     cov.uu=covb[ix12,ix12]
  }else{
     u.nz=NULL
     std.err.u = NULL
     cov.uu=NULL
  }  
  output=list(a=b1[ix10],b.nz = b.nz,u.nz=u.nz, std.err.a = std.err[ix10],std.err.b = std.err.b,std.err.u = std.err.u,s2=s2,cov.ab=covb,cov.bb=cov.bb,cov.uu=cov.uu)
  }else{
  output=list(a=NULL,b.nz = NULL,u.nz=NULL, std.err.a = NULL,std.err.b = NULL,std.err.u = NULL,s2=NULL,cov.ab=NULL,cov.bb=NULL,cov.uu=NULL)
  }
  return(output)
}

adlasso.2type <- function(y,X,Z,pz,X0=NULL,pmax=length(y)-2){ ### main function for the fixed effect DEM model
    n = length(y) ##
    p = ncol(X)
    q = ncol(Z)
  if(is.null(X0)){ 
    Xall = cbind(X,Z)
    p0 = 1
  }else{
    Xall = cbind(cbind(X0,X),Z)
    p0 = ncol(X0)+1
  }
#    print(dim(Xall))
    ysd = sd(y)
    cc = ((1-pz)*p+pz*q)/(p+q+p0)
    w0 =c(rep(0,p0),rep(1-pz,p),rep(pz,q))
    w0 = w0/cc
    glmnet0 = glmnet(Xall,y,family='gaussian',standardize=F,intercept=T,penalty.factor=w0,pmax=pmax)
#    print('glmnet0')
    lambdatot = glmnet0$lambda
    lambda = lambdatot/cc
    nmod = length(lambdatot)
    coef0 = rbind(as.vector(glmnet0$a0),as.matrix(glmnet0$beta))
    gamma0 = abs(coef0[(p0+1):(p0+p),])
    beta0 = abs(coef0[(p0+p+1):(p0+p+q),])
    wgamma = 1/pmax(gamma0,ysd/(p+q))
    wbeta = 1/pmax(beta0,ysd/(p+q))
    wgamma = (1-pz)*(wgamma%*%diag(1/colMeans(wgamma)))
    wbeta = pz*(wbeta%*%diag(1/colMeans(wbeta)))
    wad = (1/cc)*rbind(matrix(0,p0,nmod),rbind(wgamma,wbeta))
    #wad = cc*diag(c(rep(0,p0),rep(1-pz,p),rep(pz,q)))%*%wad
#        print('glmnet01')
    modall = lapply(1:nmod, function(ii) glmnet(Xall,y,family='gaussian',standardize=F,intercept=T,penalty.factor=wad[-1,ii],pmax=pmax,lambda=lambdatot[ii]))
#    print('glmnet1')
    a0all = sapply(modall,getElement, name='a0')
    coefall = rbind(a0all,as.matrix(as.data.frame(lapply(lapply(modall,getElement, name='beta'),as.matrix))))
    dfall = sapply(modall,getElement, name='df')   
    res = matrix(y,n,nmod)-cbind(1,Xall)%*%coefall
    ve = colSums(as.matrix(res)^2)/(n-dfall)
    bic = log(ve) + (log(n)/n)*dfall
    id.opt = which.min(bic)[1]
#    mod.final = modall[[id.opt]]
#    id.opt = which.min(bic)[1]
    coef.opt=coefall[,id.opt]
#    print('glmnet11')
#    print(id.opt)
#    print(summary(lambdatot))
#    print(lambdatot[id.opt])
    cov.coef=cov.adlasso.2type(cbind(1,Xall),coef.opt,wad[,id.opt],lambda=lambdatot[id.opt],s2=ve[id.opt],p0,p)
    mod.final = list(a0=coef.opt[1:p0],b=coef.opt[(1+p0):(p0+p)],u=coef.opt[(p0+p+1):(p0+p+q)],wad=wad[,id.opt],lambda=lambdatot[id.opt],cov.coef=cov.coef,ve=ve[id.opt])
    return(list(model=mod.final,bic.min = bic[id.opt],bic=bic,lambda.seq=lambdatot))
}
 
 
adlasso.2type.fixedlambda <- function(lambda,y,X,Z,pz,X0=NULL,pmax=length(y)-2){ ### main function for the fixed effect DEM model ## pz = 0.5
    n = length(y) ##
    p = ncol(X)
    q = ncol(Z)
  if(is.null(X0)){ 
    Xall = cbind(X,Z)
    p0 = 1
  }else{
    Xall = cbind(cbind(X0,X),Z)
    p0 = ncol(X0)+1
  }
    ysd = sd(y)
    cc = ((1-pz)*p+pz*q)/(p+q+p0)
    w0 =c(rep(0,p0),rep(1-pz,p),rep(pz,q))
    w0 = w0/cc
    glmnet0 = glmnet(Xall,y,family='gaussian',standardize=F,intercept=T,penalty.factor=w0,pmax=pmax,lambda=lambda)
    coef0 = c(as.vector(glmnet0$a0),as.vector(glmnet0$beta)) ## class(coef0) ## length(coef0)
    gamma0 = abs(coef0[(p0+1):(p0+p)])
    beta0 = abs(coef0[(p0+p+1):(p0+p+q)])
    wgamma = 1/pmax(gamma0,ysd/(p+q))
    wbeta = 1/pmax(beta0,ysd/(p+q))
    wgamma = (1-pz)*wgamma/mean(wgamma)
    wbeta = pz*wbeta/mean(wbeta)
    wad = (1/cc)*c(rep(0,p0),wgamma,wbeta)
    mod.opt = glmnet(Xall,y,family='gaussian',standardize=F,intercept=T,penalty.factor=wad[-1],pmax=pmax,lambda=lambda)
    coef.opt = c(as.vector(getElement(mod.opt,name='a0')),as.vector(getElement(mod.opt,name='beta')))
    df.opt = as.numeric(getElement(mod.opt, name='df'))   
    res = y-cbind(1,Xall)%*%coef.opt
    ve = sum(res^2)/(n-df.opt)
    cov.coef=cov.adlasso.2type(cbind(1,Xall),coef.opt,wad,lambda=lambda,s2=ve,p0,p)
    mod.final = list(a0=coef.opt[1:p0],b=coef.opt[(1+p0):(p0+p)],u=coef.opt[(p0+p+1):(p0+p+q)],wad=wad,lambda=lambda,cov.coef=cov.coef,ve=ve)
    return(list(model=mod.final,lambda=lambda))
}
 




medH.L2fixed <- function(mod,e2m,res.test=NA,p.adj.method=NA,pval.cut=1){ ## calculate the direct and indirect effects of the fixed model
  y = mod[['y']]
  X = mod[['X']]
  n = nrow(X)
  X0 = mod[['X0']]
  Z = mod[['Z']]
  u = as.matrix(mod[['u']])
  zu = Z%*%u    
  b = as.matrix(mod[['b']])
  ixnz = (b!=0)
  ns = sum(ixnz)
#  print(ns)
  if(!is.na(res.test)){
     if(is.na(p.adj.method)){
     ixsig = res.test[['med.individual']]$padj<=pval.cut
     }else{
     ixsig = res.test[['med.individual']][[p.adj.method]]<=pval.cut
     }
     ns = sum(ixsig)
  }
#  print(ns)
  if(ns>0){
    bnz = b[ixnz]
    Xnz = X[,ixnz]
    Xnz = matrix(Xnz,nrow=n)
    if(!is.na(res.test)){
        bnz=bnz[ixsig]
        Xnz = Xnz[,ixsig]
        e2m=e2m[ixsig]
    }
    bnz = matrix(bnz,nrow=ns)
    Xnz = matrix(Xnz,nrow=n)
    Ball =  do.call('cbind',lapply(e2m, getElement,name='coef'))
    zumed = Z%*%(Ball%*%bnz)##zumtx%*%bnz
  }else{
    zumed = 0
    Ball = NULL
  }
  zuy = zu + zumed
  v.tot = sum(zuy^2)#/n  ### all variance components are estimated up to a constant
  v.med = sum(zumed^2)#/n
  v.dir = sum(zu^2)#/n
  p.med.pure = v.med/v.tot 
  return(c(pmed.pure = p.med.pure,v.tot=v.tot,v.med=v.med,v.dir=v.dir,n.med=ns,pval.cut=pval.cut,n.direct=sum(u!=0)))
}




###' testing the effects for MedFix

devTestAdlasso <- function(dv,df){
    return(pchisq(dv,df,lower.tail=F))
}

chi2TestAdlasso <- function(b,covb){
    if(is.null(b)|is.null(covb)){
        pval=1 ## return 1 if there is no non-zero estimated coefs
    }else{
    dfb = length(b)
    b = matrix(b,ncol=1)
    if(!is.matrix(covb)) covb = matrix(covb,dfb,dfb)
    qb = sum(b*(myginv(covb)%*%b))
    pval = pchisq(qb,dfb,lower.tail=F)
    }
    return(pval)
}


testMedFix <- function(mod,e2m,p.adj.method='holm'){
    pval.direct = chi2TestAdlasso(mod[['cov.coef']][['u.nz']],mod[['cov.coef']][['cov.uu']])    
    idnz = which(mod[['b']]!=0)
    if(length(idnz)>0){
    pval.m2y = zTestAdlasso(mod[['cov.coef']][['b.nz']],mod[['cov.coef']][['std.err.b']])
    pval.e2m = unlist(lapply(e2m, function(xx) devTestAdlasso(xx[['dv']],xx[['s0']])))
#    pval.e2m = unlist(lapply(e2m, function(xx) chi2TestAdlasso(xx[['cov']][['b.nz']],xx[['cov']][['cov.bb']])))
    pval.med = data.frame(id=idnz,e2m=pval.e2m,m2y=pval.m2y)
    pval.med$e2m2y = pmin(1,pmax(pval.med$m2y,pval.med$e2m))
#    pval.med$e2m2y = pval.med$m2y+pval.med$e2m
    nmed = length(idnz)
    pval.med$padj = myPadjust(pval.med$e2m2y, method = p.adj.method) #pmin(1, nmed*pval.med$e2m2y)
#    pval.med$padj = sapply(pval.med$bonf, function(xx) min(1,sum(pval.med$bonf[pval.med$bonf<=xx])))
    output = list(direct=pval.direct,indirect=min(pval.med$e2m2y),total=min(pval.direct,min(pval.med$e2m2y)),med.individual=pval.med)
    }else{
    output = list(direct=pval.direct,indirect=1,total=pval.direct,med.individual=NULL)
    }
    return(output)
}




### functions that use foreach function, and need to source the code file itself

