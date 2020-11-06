###-----------------------------------------------------------
###
### Implementation of available methods for inferring microbiome network
###
### Kun Yue, yuek@uw.edu
### Date: 2020/11/02
### 
### (modify the file paths before using on your computer)
###-----------------------------------------------------------

library(glasso)
library(DescTools)
library(ccrepe)


##--------------------
## utility functions
##--------------------
# convert count matrix to compositional matrix
count_to_comp = function(W){
  X = sweep(W,1,STATS = rowSums(W), FUN='/') 
  return(X)
}

# compute FP and TP (adapted from COAT source code)
calTprFpr <- function(sigmaTrue, sigmaHat, eps = 1e-11){
  # modify to only account for off-diagonal entries (restrict to lower triangular matrix)
  index = lower.tri(sigmaTrue, diag=F)
  indTrueZero <- abs(sigmaTrue) < eps
  indTrueNonzero <- abs(sigmaTrue) >= eps
  indTestNonzero <- abs(sigmaHat) >= eps
  nTrueSparse <- sum((indTrueNonzero & indTestNonzero)[index])
  tpr <- nTrueSparse/sum(indTrueNonzero[index])
  nFalseSparse <- sum(indTrueZero[index] & indTestNonzero[index])
  fpr <- nFalseSparse/sum(indTrueZero[index])
  return(list(tpr = tpr, fpr = fpr))
}

# obtain ROC of correlation graph recovery, based on thresholding correlation matrix
get_cov_ROC = function(Cor, target){
  if(any(diag(Cor)!=1)){
    Cor = diag(1/sqrt(diag(Cor))) %*% Cor %*%diag(1/sqrt(diag(Cor)))
  }
  
  cov_tp <- cov_fp <- NULL
  for(thresh in seq(0, 1, length.out = 40)){
    tmp = calTprFpr(sigmaTrue = target, sigmaHat = (abs(Cor)>thresh))
    cov_tp = c(cov_tp,tmp$tpr)
    cov_fp = c(cov_fp, tmp$fpr)
    
  }
  
  return(list(fp = cov_fp, tp = cov_tp))
  
}

# obtain ROC of inverse covariance graph recovery, based on estimated covariance matrix and glasso  
find_glasso_ROC = function(covariance, lambda_seq, target_graph_inv){
  tmp = huge(covariance, lambda=lambda_seq, method='glasso')
  ROC = huge.roc(tmp$path, target_graph_inv)
  
  return(ROC) 
}

##--------------
## main function to compute ROC for the seven methods
##--------------

compare_methods = function(data_rep, # the collection of data matrixs, data[[1,k]] is composition, data[[2,k]] is count matrix, run for one repetition at a time
                           est_mat = c('covariance', 'precision')[1],
                           method = c('CoNet', 'SparCC','CCLasso', 'COAT', 'SpiecEasi', 'gCoDa','Spring')[1],
                           target_graph_cov, target_graph_inv, # input the targeted graph
                           option, 
                           lambda_seq # penalty parameter sequence pre-specified, in order to make ROC under different replications comparable
){
  n = nrow(data_rep[[1,1]])
  p = ncol(data_rep[[2,1]])
  
  diag(target_graph_cov) <- diag(target_graph_inv)<-0
  
  if (grepl(est_mat,'covariance',ignore.case=TRUE)){
    
    if (sum(grepl(method,c('CoNet', 'SparCC', 'CCLasso', 'COAT'),ignore.case=TRUE))==0) stop('target is Covariance graph, specified method not targeting the correct graph')
    
    if(grepl(method,'CoNet',ignore.case=TRUE)){
      
      #--------------
      # CoNet (also named as ReBoot)
      #--------------
      ## (input data should be p by n)
      
      ## use ccrepe package for single correlation with permulation and renormalization
      
      time = sum(system.time({
        conet_res = ccrepe(x = data_rep[[1,1]], sim.score = cor) # default 1000 iterations for p values
      })[2:3])
      
      pvals = conet_res$p.values 
      pvals_FDR = conet_res$q.values
      
      fp_null = list(fp_null = mean(pvals<0.05, na.rm=T),
                       fp_null_FDR = mean(pvals_FDR<0.05, na.rm=T))
        
      
      # constructing the ROC under alternative by thresholding the p values, for covariance graph only
      
      # pvals have NA columns/rows due to insufficient data/too many zeros, for those entries we supply pvals=2 to denote non-edge
      pvals[is.na(pvals)] <-2
      cov_tp <- cov_fp <- NULL
      for(p_thresh in seq(0, 1, length.out = 40)){
        tmp = calTprFpr(sigmaTrue = target_graph_cov, sigmaHat = (pvals<p_thresh))
        cov_tp = c(cov_tp,tmp$tpr)
        cov_fp = c(cov_fp, tmp$fpr)
      }
      
      ROC_cov = list(tp = cov_tp, fp = cov_fp)
      ROC_inv = NULL

      ROC = list(ROC_cov = ROC_cov, ROC_inv = ROC_inv, fp_null = fp_null, time = time)
    
    }else if(grepl(method,'SparCC',ignore.case=TRUE)){
      #--------------
      # SparCC 
      #--------------
      library(SpiecEasi)
      ## take in count data matrix, n by p
      
      fp_null<-fp_null_FDR <- NULL
      ROC_cov <- ROC_inv <- NULL
      
      if(option$hypothesis == 'null'){
        time = sum(system.time({
          sparcc_bootstrap_res = sparccboot(data = data_rep[[2,1]],
                                            R=1000,                                                      # number of bootstraps, maybe this will make p value noisy? with 200 reps, p value of 0.05 should vary by sqrt(0.05*0.95/200)*1.96
                                            ncpus = 4,
                                            sparcc.params = list(iter = 20, inner_iter = 10, th = 0.1)) # this is very slow
          
          # obtain p value for sparCC, use it for null model evaluation
          pval_cov =  pval.sparccboot(sparcc_bootstrap_res)$pvals
          
        })[2:3])
        
        fp_null = mean(pval_cov<0.05, na.rm=T)
        fp_null_FDR = mean(p.adjust(pval_cov, method='fdr')<0.05, na.rm=T)
        
      }else{
        time = sum(system.time({
          # use glasso for alternative model estimation and evaluation
          sparcc_res =  sparcc(data = data_rep[[2,1]],iter = 20, inner_iter = 10, th = 0.1)
        })[2:3])
        
        ROC_cov = get_cov_ROC(sparcc_res$Cor, target_graph_cov)
        ROC_inv = find_glasso_ROC(sparcc_res$Cor, lambda_seq, target_graph_inv)
        
      }

      ROC = list(ROC_cov = ROC_cov, ROC_inv = ROC_inv, 
                 fp_null = list(fp_null = fp_null, fp_null_FDR = fp_null_FDR) ,
                 time = time)
      
    }else if(grepl(method,'CCLasso',ignore.case=TRUE)){
      #--------------
      # CClasso
      #--------------
      ## (obtained from Github: https://github.com/huayingfang/CCLasso)
      source("../CCLasso-master/R/cclasso.R");
      
      # for CCLasso, the input need to be zero-corrected. so add peudocount 1 if there are zero counts
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          print('zero counts corrected by adding pseudocount=1')
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      time= sum(system.time({
        res_ccl_count <- cclasso(x = data_rep[[1,1]],
                                 counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                                 k_cv = 3,               # cv folds, min=2
                                 lam_int = c(1e-6, 3),   # tuning parameter value range;
                                 k_max=20, n_boot =2)    # input format n by p
        
      })[2:3])
      
      ROC_inv = find_glasso_ROC(res_ccl_count$cor_w, lambda_seq, target_graph_inv)
      ROC_cov = get_cov_ROC(Cor =res_ccl_count$cor_w, target = target_graph_cov )
      
      # separately compute the false positive rate under null model, for Covariance
      fp_null = calTprFpr(sigmaHat = res_ccl_count$cor_w, sigmaTrue = matrix(0, p, p))$fp
      
      ROC = list(ROC_cov = ROC_cov, ROC_inv = ROC_inv, fp_null = fp_null, time = time)
      
    }else if (grepl(method,'COAT',ignore.case=TRUE)){
      #--------------
      # COAT 
      #--------------
      ## yuanpeicao/COAT
      source("../COAT-master/COAT-master/coat.R")
      
      # peudocounts need to be added in advance for zero counts, and composition matrix adjusted accordingly
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      
      time = sum(system.time({
        coat_res = coat(data_rep[[1, 1]], nFolder=5, soft=1) # x is n by p data matrix, need to be compositional and zero adjusted
      })[2:3])
      
      
      ROC_inv = find_glasso_ROC(coat_res$corr, lambda_seq, target_graph_inv)
      ROC_cov = get_cov_ROC(Cor = coat_res$corr, target = target_graph_cov)
      fp_null = calTprFpr(sigmaHat = coat_res$corr, sigmaTrue = matrix(0, p, p))$fp
      
      ROC = list(ROC_cov = ROC_cov, ROC_inv = ROC_inv, fp_null = fp_null, time=time)
      
      
    }else{
      stop('no method matched')
    }
    
    
  }else if(grepl(est_mat,'precision',ignore.case=TRUE)){
    if(grepl(method,'SpiecEasi',ignore.case=TRUE)){
      #--------------
      # SPIEC-EASI
      #--------------
      # for inv-covariance graph recovery
      library(pulsar)
      if(option$hypothesis=='null'){
        time = sum(system.time({
          # input non-normalized, count data, without adding pesudocounts
          Spiec_network <- spiec.easi(data_rep[[2, 1]], method='glasso', # choose from 'mb' for neighbourhood selection, or 'glasso'
                                      # if to perform model selection
                                      pulsar.select = T,
                                      pulsar.params = list(
                                        thresh=0.05,                     # Threshold for StARS criterion.
                                        subsample.ratio=0.8,             # Subsample size for StARS.
                                        rep.num = 20),                   # Number of subsamples for StARS.
                                      lambda.min.ratio= min(lambda_seq)/max(lambda_seq), 
                                      lambda.max = max(lambda_seq), nlambda=length(lambda_seq)
          ) 
        })[2:3])
        
        # the result should be a solution path over lambda; the path is adjacency matrix with diagonal being zero
        Spiec_ROC_res = huge::huge.roc(path = Spiec_network$est$path, 
                                       theta = target_graph_inv, verbose=FALSE)
        
        # also need to get fp at optimal tuning value for null model
        precision = getOptNet(Spiec_network)
        fp_inv_null = calTprFpr(sigmaHat = precision, sigmaTrue = matrix(0, p, p))$fp
        
        
      }else{
        # if for ROC computation, we can skip model selection and this is faster
        
        # here to force the lambda sequence we supplied, and to remove the adding 1 effect under unnecessary cases (alt2), we run clr+glasso
        
        for( k in 1:ncol(data_rep)){
          if (any(data_rep[[2,k]]==0)){
            data_rep[[2,k]] = data_rep[[2,k]]+1
            data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
          }
        }
        
        input = data_rep[[1,1]]
        tmp = clr(input, mar=1)
        time = sum(system.time({
          use_clr = huge(t(tmp), method='glasso', lambda = lambda_seq)
        })[2:3]
        )
        
        Spiec_ROC_res = huge.roc(use_clr$path, theta = target_graph_inv) 
        
        fp_inv_null= NULL
        
      }


      ROC_inv = Spiec_ROC_res



      ROC_cov <- ROC_cov_pathbased <- NULL
      
      
      # covariance graph recovery depreciated (encounter error with some simulation settings)
      if(F){
                    
        # thresholding the optimal cov matrix for ROC
      
        # Spiec_cov = as.matrix(getOptCov(Spiec_network))ROC_cov = get_cov_ROC(Spiec_cov, target_graph_cov)
      
              # or use the cov path corresponding to the lambda path, and compute ROC. to eliminate very small values, make Cov into Cor and threhsold at 1e-11
              ROC_cov_pathbased = 
                huge::huge.roc(path = 
                                 lapply(Spiec_network$est$cov, function(x){
                                   x = diag(1/sqrt(diag(x)))%*% x %*% diag(1/sqrt(diag(x)))# transform into correlation matrix
                                   diag(x) = 0
                                   tmp = (abs(x)>1e-11)*1
                                   tmp}), 
                               theta = target_graph_cov)
      }
      


      ROC = list(ROC_cov = ROC_cov, ROC_cov_pathbased=ROC_cov_pathbased,
                 ROC_inv = ROC_inv, fp_inv_null = fp_inv_null, time=time)

      
    }else if(grepl(method,'gcoda',ignore.case=TRUE)){
      #--------------
      # gCoda
      #--------------
      source('../gCoda-master/R/gcoda.R')
      
      # input zero corrected compositional data
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      time = sum(system.time({
        gcoda_network <- gcoda(data_rep[[1,1]], counts = F, pseudo = 1, 
                             lambda.min.ratio=1e-3,       # these will be ignored
                             nlambda= length(lambda_seq), # this has to be the length of provided lambda_seq
                             lambda = lambda_seq);        # it seems gcoda sometimes automatically include more lambda values
      })[2:3])

      # covariance based graph recovery depreciated
      ROC_cov <- ROC_cov_pathbased <- NULL
      if(F){
        # based on thresholding the optimal cov (inverse of optimal icov)
        gcoda_cov = solve(gcoda_network$opt.icov)
        ROC_cov = get_cov_ROC(Cor = gcoda_cov, target_graph_cov)
        
        # based on the glasso path generated cov
        ROC_cov_pathbased = 
          huge::huge.roc(path = 
                           lapply(gcoda_network$icov, function(x){
                             x = solve(x)
                             x = diag(1/sqrt(diag(x)))%*% x %*% diag(1/sqrt(diag(x)))# transform into correlation matrix
                             diag(x) = 0
                             tmp = (abs(x)>1e-11)*1
                             tmp}), 
                         theta = target_graph_cov)
        
        
      }
      
      gcoda_ROC_res = huge::huge.roc(path=gcoda_network$path, theta = option$Sigma_list$A_inv, verbose=F)
      ROC_inv = gcoda_ROC_res
      
      
      precision = gcoda_network$opt.icov
      fp_inv_null = calTprFpr(sigmaHat = precision, sigmaTrue = matrix(0, p, p))$fp
      
      ROC = list(ROC_cov = ROC_cov, ROC_cov_pathbased = ROC_cov_pathbased,
                 ROC_inv = ROC_inv, fp_inv_null = fp_inv_null, lambda_seq = gcoda_network$lambda,
                 time = time)
      
    }else if(grepl(method,'Spring',ignore.case=TRUE)){
      #---------------
      # SPRING
      #---------------
      # devtools::install_github("irinagain/mixedCCA")
      # devtools::install_github("GraceYoon/SPRING")
      library(SPRING)

      # supply uncorrected compositional data with n by p data matrix. It will be mclr transformed inside their function
      if(option$hypothesis == 'null'){
        # only do lambda selection for null hypothesis
        time = sum(system.time({
          fit.spring <- SPRING(data_rep[[1,1]], 
                               quantitative = F, # F means input is compositional
                               lambdaseq = lambda_seq, 
                               # nlambda = 20, 
                               ncores = 1,
                               subsample.ratio = 0.8, rep.num = 20 # these are for tuning parameter selection
          )
        })[2:3])
      }else{
        time = sum(system.time({
          fit.spring <- SPRING(data_rep[[1,1]], 
                               quantitative = F, 
                               lambdaseq = lambda_seq, 
                               # nlambda = 20, 
                               ncores = 1,
                               subsample.ratio = 0.1, rep.num = 2 
          )
        })[2:3])
      }

      spring_ROC_res = huge::huge.roc(fit.spring$fit$est$path, theta = option$Sigma_list$A_inv, verbose = F)
      ROC_inv = spring_ROC_res
      
      
      opt.K <- fit.spring$output$stars$opt.index
      adj.icov <- as.matrix(fit.spring$fit$est$path[[opt.K]]) # adjacency matrix
      fp_inv_null = calTprFpr(sigmaHat = adj.icov, sigmaTrue = matrix(0, p, p))$fp
      
      
      # threshold the optimal cov (depreciated)
      # Estimated partial correlation coefficient, same as negative precision matrix.
      precision <- -as.matrix(SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode = 'maxabs'))
      diag(precision)<-1
      ROC_cov <- ROC_cov_pathbased <-NULL
      
       
      if(F){
         ROC_cov = get_cov_ROC(solve(precision), target_graph_cov)
         # use the glasso path generated cov path
         ROC_cov_pathbased = 
           huge::huge.roc(path = lapply( fit.spring$output$est$beta,
                                         function(x){
                                           x = -as.matrix(SpiecEasi::symBeta(x, mode = 'maxabs'))  # this is the coefficient matrix from neighbourhood selection
                                           diag(x) = 1
                                           x = solve(x) # compute covariance
                                           x = diag(1/sqrt(diag(x)))%*% x %*% diag(1/sqrt(diag(x)))# transform into correlation matrix
                                           diag(x) = 0
                                           tmp = (abs(x)>1e-11)*1
                                           tmp}), 
                          theta = target_graph_cov)
       
      }
      
      
      ROC = list(ROC_cov = ROC_cov, ROC_cov_pathbased = ROC_cov_pathbased,
                 ROC_inv = ROC_inv, fp_inv_null = fp_inv_null, time=time)
      

    }else{
      stop('no method matched')
    }
    
  }
  
  
  return(ROC)
}
