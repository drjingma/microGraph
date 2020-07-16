##----------------------------------------------
##
## Evaluation methods
##
##----------------------------------------------

estNet = Sigma + sample(c(0.01, 1, 0), size=p*p, replace=T)
trueNet = Sigma + sample(c(0.01, 1, 0), size=p*p, replace=T)

eval_MSE = function(estNet, trueNet){
  mean(((estNet - trueNet)[lower.tri(estNet, diag=F)])^2)
}


eval_presence = function(estNet, trueNet, 
                         threshold = 0.3  # in SparCC paper, correlation is thresholded at 0.3 as final network
){ 
  estNet[estNet<threshold] <-0
  p=dim(estNet)[1]
  total_edge = p*(p-1)/2
  estsignpos = estNet[lower.tri(estNet, diag=F)]>0
  truesignpos = trueNet[lower.tri(estNet, diag=F)]>0
  estsignneg = estNet[lower.tri(estNet, diag=F)]<0
  truesignneg = trueNet[lower.tri(estNet, diag=F)]<0
  
  TP = sum(estsignpos*truesignpos) + sum(estsignneg*truesignneg) 
  TN = sum(estNet[lower.tri(estNet, diag=F)]==0 & trueNet[lower.tri(estNet, diag=F)]==0)
  FP = sum(estNet[lower.tri(estNet, diag=F)]!=0 & trueNet[lower.tri(estNet, diag=F)]==0) + 1/2 * (sum(estsignpos*truesignneg) + sum(truesignpos * estsignneg))
  FN = sum(estNet[lower.tri(estNet, diag=F)]==0 & trueNet[lower.tri(estNet, diag=F)]!=0) + 1/2 * (sum(estsignpos*truesignneg) + sum(truesignpos * estsignneg)) 
  
  true_edge = sum(trueNet[lower.tri(estNet, diag=F)]!=0)
  return(list(TP=TP, TN=TN, FP=FP, FN=FN, true_edge = true_edge))
}

eval_mean_absolute_error =function(estNet, trueNet){
  mean((abs(estNet - trueNet)[lower.tri(estNet, diag=F)])^2)
}


eval_AUC = function(estNet, trueNet # only says recover nonzero entries, so do not distinguish signs
){
  library(pROC)
  tmp=roc(as.vector(as.factor(trueNet[lower.tri(estNet, diag=F)]!=0)),as.vector(estNet[lower.tri(estNet, diag=F)]))
  plot(tmp)
  auc(tmp)[1]
}


tmp=eval_AUC(estNet, trueNet)
