ModifiedPCOMerged = function(Z.mat,Sigma,SigmaO,p.method,method = "davies"){
  ## avoid Z.mat is a row vector which won't work later
  ## Z.mat should be a column vector object, not a 1*K matrix 
  
  ##########
  eigen.res = eigen(Sigma)
  lambdas = eigen.res$values
  eigen.vec = eigen.res$vectors
  
  if(p.method != "PCO"){
    lambdas = lambdas[lambdas>1]
    eigen.vec = eigen.vec[, which(lambdas>1)]
  }
  
  K = length(lambdas)
  w.vec = 1/lambdas
  
  z.tmp = t(Z.mat)
  PC.vec = t(eigen.vec)%*%z.tmp
  
  
  ## 1: PCMinP
  PC.std = (PC.vec)^2 * w.vec
  PC.p = pchisq(PC.std,df=1,lower.tail = FALSE)
  p.min = apply(PC.p, 2, min)
  p.PCMinP = 1- (1-p.min)^K
  
  ## 2: PCFisher
  PC.Fisher.stat = -2*colSums(log(PC.p))
  p.PCFisher = pchisq(PC.Fisher.stat,df=2*K,lower.tail = FALSE)
  
  ## 3: PCLC
  PCLC.stat = colSums(PC.vec * w.vec)
  PCLC.stat.std = PCLC.stat^2/sum(w.vec)
  p.PCLC = pchisq(PCLC.stat.std,df=1,lower.tail = FALSE)
  
  ## 4: WI
  T.WI = colSums((PC.vec)^2)
  
  if(method=="davies"){
    p.WI = sapply(T.WI, function(x) davies(x, lambdas)$Qq)
    if(!all(p.WI > 0)){
      ind.liumod = !(p.WI > 0)
      snp.liumod = names(p.WI[ind.liumod])
      p.WI[snp.liumod] = liumod(T.WI[snp.liumod], lambdas)
    }
  }else if(method=="liu"){
    p.WI = liu(T.WI,lambdas)
  }else if(method=="liumod"){
    p.WI = liumod(T.WI,lambdas)
  }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
  
  ## 5: Wald
  T.Wald = colSums(PC.std)
  p.Wald = pchisq(T.Wald,df=K,lower.tail=FALSE)
  
  
  ## 6: VC
  T.VC = colSums((PC.vec * w.vec)^2)
  wt.VC = sort(w.vec, decreasing = TRUE)
  
  if(method=="davies"){
    p.VC = sapply(T.VC, function(x) davies(x, wt.VC)$Qq)
    if(!all(p.VC > 0)){
      ind.liumod = !(p.VC > 0)
      snp.liumod = names(p.VC[ind.liumod])
      p.VC[snp.liumod] = liumod(T.VC[snp.liumod], wt.VC)
    }
  }else if(method=="liu"){
    p.VC = liu(T.VC, wt.VC)
  }else if(method=="liumod"){
    p.VC = liumod(T.VC, wt.VC)
  }else {stop("Please specify a valid method: davies, liu, liumod.\n")}
  
  
  ## PCO
  T.minp = apply(rbind(p.PCMinP,p.PCFisher,p.PCLC, p.WI,p.Wald,p.VC), 2, min)
  
  K.X = dim(SigmaO)[1]
  if(K.X!=6){
    stop("The dimension of SigmaO is not 6!!")
  }
  
  upper_bound = rep(Inf,K.X); mv.mean = rep(0, K.X)
  pval = sapply(qnorm(T.minp),
                 function(x) 1-mvtnorm::pmvnorm(lower = rep(x, K.X), upper = upper_bound,
                                                mean = mv.mean, sigma = SigmaO)[1])
  return(pval)
}
