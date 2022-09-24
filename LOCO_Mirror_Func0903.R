######Select-Estimate-Predict#######
SEP <- function(Xtr, Ytr, Xte, Yte){
  #####Step 1: Use training set to find selected S#####
  M = cv.glmnet(Xtr,Ytr)
  
  # ind0 = which(M$nzero>1 & M$nzero<(nrow(Xtr)*0.8))
  # ind = which.min(M$cvm[ind0]) #remove the case: no variable is significant
  # beta = M$glmnet.fit$beta[,ind0[ind]]
  # nonbeta = which(beta!=0)
  
  ind = which.min(M$cvm)
  beta = M$glmnet.fit$beta[,ind]
  nonbeta = which(beta!=0)
  
  Yfit = Xtr%*%beta
  
  Ypred = Xte%*%beta
  #####Step 2: Use traing set to do LSE on S######
  # if(length(Ytr)>length(nonbeta)){
  # lmMod = lm(Ytr~Xtr[,nonbeta])
  # 
  # #####Step 3: Use test set to predict#####
  # Ypred = lmMod$coefficients[1]+Xte[,nonbeta]%*%lmMod$coefficients[-1]
  # # Ypred = predict(lmMod,as.data.frame(Xte[,nonbeta]))
  # }
  
  return(list(S=nonbeta,predEr = abs(Yte-Ypred), fitEr = abs(Ytr-Yfit), beta=beta))
  
}

###### One Split: LOCO-Mirror #######
LOCO_Mirror <- function(X,Y, q, nonzero, n, p){
  
  #####Step 0: Split data into two equal parts
  Ind1 = sample(1:n,floor(n/2))
  Ind0 = (1:n)[-Ind1]
  
  X0 = X[Ind0,];Y0 = Y[Ind0];
  X1 = X[Ind1,];Y1 = Y[Ind1];
  
  #####Step 1: Find S and do LSE on S with D0, predict error with D1
  Shat = SEP(X0, Y0, X1, Y1)
  
  # #####Step 2: Leave-one-covariate-out to repeat Step 1
  # deltaJ = sapply(Shat$S, function(t){
  #   Sjhat = SEP(X0[,-t], Y0, X1[,-t], Y1)
  #   d1j = Sjhat$fitEr - Shat$fitEr
  #   d2j = Sjhat$predEr - Shat$predEr
  #   # return(Sjhat$predEr - Shat$predEr)
  #   return(mean(d1j)*mean(d2j))
  # })
  # W = rep(0,p)
  # W[Shat$S] = deltaJ
  

  deltaJ = sapply(Shat$S, function(t){
    Sjhat = SEP(X0[,-t], Y0, X1[,-t], Y1)
    return(Sjhat$predEr - Shat$predEr)
  })
  #
  #####Step 3: Compute hatGammaS with formula (27) ######
  # gammaS = colMeans(deltaJ)
  W = rep(0,p)
  W[Shat$S] = colMeans(deltaJ)
  
  W = W*abs(Shat$beta)
  
  ####Step 4: FDR control#####
  
  t = sort(abs(W))
  
  Tt = sapply(t,function(x){(1+length(W[W<=(-x)]))/max(1,length(W[W>=x]))})
  thre = min(t[Tt<=q])
  det = which(W >= thre)
  #det = which(Z0>=thre)
  FDP = length(setdiff(det,nonzero))/max(1,length(det))
  TPR = length(intersect(det,nonzero))/length(nonzero)
  
  return(list(W=W, FDP=FDP, TPR=TPR, det = det,deltaJ = deltaJ,nonbeta=Shat$S))
}

###### LOCO BH #######
LOCO_BH <- function(deltaJ,q,nonzero,nonbeta){
  
  nd = nrow(deltaJ)
  pd = ncol(deltaJ)
  
  ttest = sapply(1:pd, function(t){
    t.test(deltaJ[,t])$statistic*sqrt(nd/(nd-1))
    })
  
  pvalue = 1-pnorm(ttest)
  
  detBH = BH(pvalue, pd, q)
  detBH = nonbeta[detBH]
  
  FDPBH = length(setdiff(detBH,nonzero))/max(1,length(detBH))
  TPRBH = length(intersect(detBH,nonzero))/max(1,length(nonzero))
  
  return(c(length(detBH),FDPBH,TPRBH))
}


###### Multiple Split: LOCO-Mirror with SDA Bagging strategy #######
# BaggingSplit <- function(X,Y, q, nonzero, n, p, B=10){
  
  # #### Repeat one-split mirror B times ####
  # OneSplitRes <- lapply(1:B,function(t){
  #   #### One-split results, can adjust the Mirror function ####
  #   LOCO_Mirror(X,Y, q, nonzero, n, p) #LOCO 
  # })
BaggingSplit <- function(OneSplitRes, q, n, p, B=10){  
  
  #### Union all detected mediators ####
  detAll = lapply(1:B, function(t){OneSplitRes[[t]]$det})
  detBag = as.data.frame(table(unlist(detAll)))
  
  #### Find those that are selected more than B/2 times ####
  detStar = as.numeric(detBag[detBag$Freq > (B/2), 1])
  
  #### Find which replication has the most intersection ####
  Len.Inter = sapply(1:B, function(t){
    length(intersect(detStar,detAll[[t]]))
  })
  k = which.max(Len.Inter)
  
  ##### output the refined results #####
  Res = c(length(OneSplitRes[[k]]$det),OneSplitRes[[k]]$FDP, OneSplitRes[[k]]$TPR)
  # rownames(Res) = c("NumDet","FDR","TPR")
  # Res = as.data.frame(t(Res))
  # Res$Method = c("SDA-Refined")
  
  return(Res)
}


###### Multiple Split: LOCO-Mirror with meanW strategy #######
# BaggingSplit2 <- function(X,Y, q, nonzero, n, p, B=10){
#   
#   #### Repeat one-split mirror B times ####
#   OneSplitRes <- lapply(1:B,function(t){
#     #### One-split results, can adjust the Mirror function ####
#     LOCO_Mirror(X,Y, q, nonzero, n, p) #LOCO 
#   })
BaggingSplit2 <- function(OneSplitRes, nonzero, q, n, p, B=10){  
  
  #### Summarize the W ########
  Wall = sapply(1:B, function(t){OneSplitRes[[t]]$W})
  W = rowMeans(Wall)
  t = sort(abs(W))
  
  Tt = sapply(t,function(x){(1+length(W[W<=(-x)]))/max(1,length(W[W>=x]))})
  thre = min(t[Tt<=q])
  det = which(W >= thre)
  #det = which(Z0>=thre)
  FDP = length(setdiff(det,nonzero))/max(1,length(det))
  TPR = length(intersect(det,nonzero))/length(nonzero)
  
  ##### output the refined results #####
  Res = c(length(det),FDP, TPR)
  # rownames(Res) = c("NumDet","FDR","TPR")
  # Res = as.data.frame(t(Res))
  # Res$Method = c("OneSplit", "Refined")
  
  return(Res)
}


###### Multiple Split: LOCO-Mirror-BH #######
MultipleBH <- function(OneSplitRes, nonzero, q, n, p, B=10){  

  pvalAll <- sapply(1:B, function(t){
    deltaJ = OneSplitRes[[t]]$deltaJ
    nd = nrow(deltaJ)
    pd = ncol(deltaJ)
    ttest = sapply(1:pd, function(t){
      t.test(deltaJ[,t])$statistic*sqrt(nd/(nd-1))
    })
    
    pvalue = rep(1,p)
    pvalue[OneSplitRes[[t]]$nonbeta] = 1-pnorm(ttest)
    
    return(pvalue)
  })

  gammin = 0.05
  gam = seq(gammin,1-0.01,0.01)
  pvalQstar <- sapply(1:p, function(t){
    Qgam = sapply(gam, function(s){
      min(1,quantile(pvalAll[t,]/s,s))
    }) 
    return(min(1,(1-log(gammin))*min(Qgam)))
  })
  detBH = BH(pvalQstar, p, q)
  
  AvaP = which(rowSums(pvalAll)!=B)
  detBH = BH(pvalQstar[AvaP], length(AvaP), q)
  detBH = AvaP[detBH]
  
  FDPBH = length(setdiff(detBH,nonzero))/max(1,length(detBH))
  TPRBH = length(intersect(detBH,nonzero))/max(1,length(nonzero))
  
  return(c(length(detBH),FDPBH,TPRBH))

}
