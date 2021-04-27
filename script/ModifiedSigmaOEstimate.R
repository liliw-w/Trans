# Estimate the SigmaMeta for PCMeta p-value computation
# This is a pre-computation step for PCMeta
# input: the correlation matrix among Z-scores
# output: the SigmaMeta matrix 
ModifiedSigmaOEstimate = function(Sigma,p.method="TruncPCO",simNum=2000,method = "davies"){
  
  X.PCMinP = rep(NA,simNum)
  X.PCFisher = rep(NA,simNum)
  X.PCLC = rep(NA,simNum)
  X.WI = rep(NA,simNum)
  X.Wald =rep(NA,simNum)
  X.VC =rep(NA,simNum)
  
  for(i in 1:simNum){
    Z.vec = t(mvtnorm::rmvnorm(1,rep(0, nrow(Sigma)), Sigma))
    
    if(!is.vector(Z.vec)){
      Z.vec = as.vector(Z.vec)
    }
    
    ##########
    eigen.res = eigen(Sigma)
    lambdas = eigen.res$values
    eigen.vec = eigen.res$vectors
    
    if(p.method != "PCO"){
      lambdas = lambdas[lambdas>0.1]
      eigen.vec = eigen.vec[, which(lambdas>0.1)]
    }
    K = length(lambdas)
    
    z.tmp = as.matrix(Z.vec)
    PC.vec = t(eigen.vec)%*%z.tmp
    
    
    ## 1: PCMinP
    PC.std = (PC.vec)^2/lambdas
    PC.p = pchisq(PC.std,df=1,lower.tail = FALSE)
    p.min = min(PC.p)
    p.PCMinP = 1- (1-p.min)^K
    
    ## 2: PCFisher
    PC.Fisher.stat = -2*sum(log(PC.p))
    p.PCFisher = pchisq(PC.Fisher.stat,df=2*K,lower.tail = FALSE)
    
    ## 3: PCLC
    w.vec = 1/lambdas
    PCLC.stat = sum(w.vec*PC.vec)
    PCLC.stat.std = PCLC.stat^2/sum(w.vec)
    p.PCLC = pchisq(PCLC.stat.std,df=1,lower.tail = FALSE)
    
    ## 4: WI
    T.WI = sum((PC.vec)^2)
    
    if(method=="davies"){
      p.WI = davies(T.WI,lambdas)$Qq
      if(p.WI ==0 | p.WI < 0){
        p.WI = liumod(T.WI,lambdas)
      }
    }else if(method=="liu"){
      p.WI = liu(T.WI,lambdas)
    }else if(method=="liumod"){
      p.WI = liumod(T.WI,lambdas)
    }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
    
    
    ## 5: Wald
    T.Wald = sum(PC.std)
    p.Wald = pchisq(T.Wald,df=K,lower.tail=FALSE)
    
    
    ## 6: VC
    T.VC = sum((w.vec*PC.vec)^2)
    
    if(method=="davies"){
      p.VC = davies(T.VC,sort(1/lambdas, decreasing = TRUE))$Qq
      if(p.VC == 0 | p.VC < 0){
        p.VC = liumod(T.VC,sort(1/lambdas, decreasing = TRUE))
      }
    }else if(method=="liu"){
      p.VC = liu(T.VC,sort(1/lambdas, decreasing = TRUE))
    }else if(method=="liumod"){
      p.VC = liumod(T.VC,sort(1/lambdas, decreasing = TRUE))
    }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
    
    
    
    X.PCMinP[i] = qnorm(p.PCMinP)
    X.PCFisher[i] = qnorm(p.PCFisher)
    X.PCLC[i] = qnorm(p.PCLC)
    X.WI[i] = qnorm(p.WI)
    X.Wald[i] = qnorm(p.Wald)
    X.VC[i] = qnorm(p.VC)
  }
  
  X.mat = cbind(X.PCMinP,X.PCFisher,X.PCLC,X.WI,X.Wald,X.VC)
  X.mat = X.mat[!rowSums(!is.finite(X.mat)),] ## remove rows containing NA, NaN,Inf
  SigmaO = cor(X.mat)
  return(SigmaO)
  
}
