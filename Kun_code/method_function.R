###-----------------------------------------------------------
###
### Implementation of available methods for computing microbiome network
###
###-----------------------------------------------------------
library(glasso)
library(DescTools)
library(ccrepe)

# filepath = 'C:\\Users\\yuek\\Dropbox\\Microbial_Networks\\microGraph' #BOX
# filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' #mac
# filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC
filepath = '~/Desktop/micro_net' #bayes

count_to_comp = function(W){
  X = sweep(W,1,STATS = rowSums(W), FUN='/') # the compositions from log-normal; different Xi for different cells
  return(X)
}

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




find_glasso_ROC = function(covariance, lambda_seq, target_graph_inv){

  
  # tmp = glassopath(s = covariance, rholist = lambda_seq, penalize.diagonal = F, trace=0)
  # 
  # tpfp = sapply(1:length(lambda_seq), function(i){
  #   prec = tmp$wi[,,i] # inverse covariance
  #   diag(prec) = 0
  #   prec = (abs(prec)>1e-11)*1
  #   calTprFpr(sigmaTrue = target_graph_inv, sigmaHat = prec)
  # })
  # 
  # tp = sort(c(0, do.call(c, tpfp[1,]), 1)) # always add this (0,0) and (1,1) point since the start/end point of ROC from (0,0) to (1,1)
  # fp = sort(c(0, do.call(c,tpfp[2,]), 1))
  # auc = AUC(x=fp , y=tp)
  # ROC = list(tp = tp, fp = fp, AUC = auc)
  
  
  
  # tmp = huge(covariance, nlambda=40, lambda.min.ratio = 0.001, method='glasso')
  
  tmp = huge(covariance, lambda=lambda_seq, method='glasso')
  ROC = huge.roc(tmp$path, target_graph_inv)
  
  return(ROC) 
}


compare_methods = function(data_rep, # the collection of data matrixs, data[[1,k]] is composition, data[[2,k]] is count matrix, run for one repetition at a time
                           est_mat = c('covariance', 'precision'),
                           method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                      'SpiecEasi', 'gCoDa','Spring'),
                           target_graph_cov, target_graph_inv, # input the targeted graph; if provide var/inverse cov, must set the diagonal to zero first
                           option,
                           lambda_seq
){
  n = nrow(data_rep[[1,1]])
  p = ncol(data_rep[[2,1]])
  
  diag(target_graph_cov) <- diag(target_graph_inv)<-0
  
  if (grepl(est_mat,'covariance',ignore.case=TRUE)){
    
    if ( sum(grepl(method,c('CoNet', 'SparCC', 'CCLasso', 'COAT'),ignore.case=TRUE))==0) stop('target is Covariance graph, specified method not targeting the correct graph')
    
    if(grepl(method,'CoNet',ignore.case=TRUE)){
      #--------------
      # CoNet (also named as ReBoot)
      #--------------
      ## (input data should be p by n)
      
      ### option 1: use the CoNet software (not in R)
      {
      # filepath_conet = 'CoNet_app_output\\Conet_output2.txt'
      # read_tab = function(filepath){
      #   
      #   processFile = function(filepath) {
      #     full = list()
      #     count=1
      #     con = file(filepath, "r")
      #     while ( TRUE ) {
      #       line = readLines(con, n = 1)
      #       if ( length(line) == 0 ) {
      #         break
      #       }
      #       full[count]=line
      #       count=count+1
      #     }
      #     close(con)
      #     return(full)
      #   }
      #   
      #   see = processFile(filepath)
      #   
      #   # search for nodes information
      #   i=1
      #   while (i <=length(see)){
      #     if(substr(see[[i]],1,6) == ';NODES'){nodehead = i}
      #     if(substr(see[[i]],1,5) == ';ARCS'){edgehead=i}
      #     i=i+1
      #   }
      #   node_info = read.table(filepath, skip=nodehead-1, nrow=edgehead-nodehead, sep='\t')
      #   edge_info = read.table(filepath, skip=edgehead-1, nrow=length(see) - edgehead+1, sep='\t', header=T)
      #   colnames(edge_info)[1]<-'Node2'
      #   edge_info = cbind(Node1 = rownames(edge_info), edge_info)
      #   rownames(edge_info)<-NULL
      #   
      #   return(list(node_info = node_info, edge_info = edge_info))
      # }
      # read_tab(filepath_conet) # still need to determine how to use the edge informations.
      
      
      
      
      ### option 2: use its R version without some features
      
      # install.packages('gsl')
      # library(gsl)
      # file.path(R.home("gsl"), "Makeconf")
      # 
      # 
      # if (!requireNamespace("BiocManager", quietly = TRUE))
      #   install.packages("BiocManager")
      # 
      # BiocManager::install("DirichletMultinomial")
      # 
      # library(devtools)
      # install_github("hallucigenia-sparsa/seqgroup")
      
      
      # library(seqgroup)
      # #reference: https://hallucigenia-sparsa.github.io/seqgroup/reference/barebonesCoNet.html
      # #input data needs p by n; need to name the taxa
      #
      # # not working
      # net_CoNet = barebonesCoNet(abundances = t(data_rep[[2,1]]), 
      #                            methods = c("spearman",'pearson', 'bray',"kld"), 
      #                method.num.T = 4, pval.T = 0.05,
      #                #init.edge.num is using default value (sqrt(p)); this is the number of top and bottom edges to initially keep for later testing
      #                init.edge.num = 50,
      #                pval.cor = FALSE, # do not use usual correlation test
      #                permut = T, renorm = T,permutandboot = T, # their permutation/bootstrap and renormalization will remove spurious correlation
      #                iters = 100, bh = TRUE,
      #                pseudocount = 1, # counts added to zeros when taking log
      #                plot = F, verbose = F) 
      # 
      # # not working
      # net_CoNet = barebonesCoNet(abundances = t(data_rep[[2,1]]), methods = c("spearman", 'kld'), 
      #                            method.num.T = 2, pval.T = 0.05,
      #                            pval.cor = T, # do not use usual correlation test
      #                            init.edge.num = 20,
      #                            permut = F, renorm = F,permutandboot = F, # their permutation/bootstrap and renormalization will remove spurious correlation
      #                            iters = 100, bh = TRUE,
      #                            pseudocount = 1e-11, # counts added to zeros when taking log
      #                            plot = T, verbose = F) 
      # 
      # # this one works
      # net_CoNet = barebonesCoNet(abundances = t(data_rep[[2,1]]), methods = c("spearman"),
      #                            method.num.T = 2, pval.T = 0.05,
      #                            pval.cor = F, # do not use usual correlation test
      #                            # init.edge.num = 40,
      #                            permut = F, renorm = T,permutandboot = T, # their permutation/bootstrap and renormalization will remove spurious correlation
      #                            iters = 100, bh = TRUE,
      #                            pseudocount = 1e-11, # counts added to zeros when taking log
      #                            plot = T, verbose = F)
      # 
      # net_CoNet = barebonesCoNet(abundances = t(data_rep[[2,1]]), methods = c("spearman"),
      #                            method.num.T = 2, pval.T = 0.05,
      #                            pval.cor = T, # do not use usual correlation test
      #                            init.edge.num = 300,
      #                            permut = F, renorm = F,permutandboot = F, # their permutation/bootstrap and renormalization will remove spurious correlation
      #                            iters = 100, bh = TRUE,
      #                            pseudocount = 1e-11, # counts added to zeros when taking log
      #                            plot = T, verbose = F)
      # 
      # 
      # # the output is am igraph object with edge weights;
      # # here consider tuning parameter init.edge.num ( methods, method.num are just set fixed, though can also be changed)
      # # init.
      # 10^(-edge_attr(net_CoNet)$weight) # this is the p value for each edge output
      # network_matrix_CoNet = as_adjacency_matrix(net_CoNet,type='both', attr = "weight", sparse=F)
      # 
      # 
      # plot(net_CoNet)
      }
      
      
      ### working option: use ccrepe package for single correlation with permulation and renormalization
      time = sum(system.time({
        conet_res = ccrepe(x = data_rep[[1,1]], sim.score = cor) # default 1000 iterations for p values
      })[2:3])
      
      pvals = conet_res$p.values 
      pvals_FDR = conet_res$q.values
      
      fp_null = list(fp_null = mean(pvals<0.05, na.rm=T),
                       fp_null_FDR = mean(pvals_FDR<0.05, na.rm=T))
        
      # constructing the RIC under alternative by thresholding the p values, for covariance graph only
      
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
      # library(devtools)
      # install_github("zdk123/SpiecEasi")
      library(SpiecEasi)
      ## take in count data matrix, n by p
      # sparcc_res = sparcc(data_rep[[2,1]], iter = 20, inner_iter = 10, th = 0.1) # need input format n by p; varying threshold (th) has minial impact on the result
      # sparcc_res$Cov # the estimated log counts covariance, the Sigma
      # sparcc_res = res$Cor # the estimated log counts correlation, the corr_mat
      
      
      {
        
        # sparcc_bootstrap_res = sparccboot(data = data_rep[[2,1]], 
        #                                   R=200,  # number of bootstraps
        #                                   ncpus = 4, 
        #                                   sparcc.params = list(iter = 20, inner_iter = 10, th = 0.1)) # this is very slow
        # 
        # pval_cov =  pval.sparccboot(sparcc_bootstrap_res)
        # sparcc_bootstrap_pval <- sparcc_bootstrap_cor <- matrix(0, p, p)
        # sparcc_bootstrap_pval[upper.tri(sparcc_bootstrap_pval)] = pval_cov$pvals
        # sparcc_bootstrap_pval = sparcc_bootstrap_pval+t(sparcc_bootstrap_pval) 
        # # sparcc_bootstrap_pval[1:10, 1:10]
        # diag(sparcc_bootstrap_pval) = 2 # the diagonal set to 2, so never detect an edge on the diagonal
        # 
        # sparcc_bootstrap_cor[upper.tri(sparcc_bootstrap_pval)] = pval_cov$cors
        # sparcc_bootstrap_cor = sparcc_bootstrap_cor+t(sparcc_bootstrap_cor)
        # diag(sparcc_bootstrap_cor) = 1
        # 
        # sparcc_bootstrap_pval_FDR = matrix(p.adjust(p = as.vector(sparcc_bootstrap_pval), method = 'fdr'), p, p, byrow=F)
        # 
        # # will changing the threshold for the matrix by hand to construct ROC
        # nlam = 20
        # get_path_raw <- get_path_FDR <-  list()
        # ind = 0
        # for (alpha in seq(0, 1, length = nlam)){
        #   ind = ind+1
        #   get_path_raw[[ind]] = (sparcc_bootstrap_pval<=alpha)*1
        #   get_path_FDR[[ind]] = (sparcc_bootstrap_pval_FDR<=alpha)*1
        #   
        # }
        # 
        # ROC_raw = huge::huge.roc(path = get_path_raw, theta = target_graph_cov)
        # ROC_FDR = huge::huge.roc(path = get_path_FDR, theta = target_graph_cov)
        # 
        # 
        # # problem with ROC_inv: do not have p values to threshold
        # ROC = list(ROC_cov = list(ROC_raw = ROC_raw, ROC_FDR = ROC_FDR),
        #            ROC_inv = NULL)
        
      }
      
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
        
        # mean(pval_cov < 0.05)
        # mean(pval_cov < 0.05-sqrt(0.05*0.95/200)*1.96) # this is still showing inflated false positive
        
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
      source(paste0(filepath, "/Kun_code/CCLasso-master/R/cclasso.R"));
      
      # for CCLasso, the input need to be zero-corrected. so add peudocount 1 if there are zero counts
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          print('zero counts corrected by adding pseudocount=1')
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      
      #       
      # 
      #       lambda_vec = sort(lambda_vec)
      #       ccl_cor_path = list()
      #       ccl_cor = list()
      #       ccl_inv = list()
      #       ind = 0
      #       for(lambda in lambda_vec){
      #         ind=ind+1
      #         res_ccl_count <- cclasso(x = data_rep[[1,1]], 
      #                                  counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
      #                                  k_cv = 2, # cv folds, min=2
      #                                  lam_int = rep(lambda), #tuning parameter value range; need to vary this for ROC curve
      #                                  k_max=200, n_boot =1);  # input format n by p
      #         
      #         # res_ccl_count$info_cv$lams # the lambda path 
      #         ccl_cor[[ind]] = res_ccl_count$cor_w
      #         # res_ccl_count$p_vals # from bootstrap, ignored in our simulation
      #       }
      #       
      #       
      #       ROC_cov = huge::huge.roc(path = lapply(ccl_cor, function(x){
      #                                               tmp = (abs(x)>1e-11)*1
      #                                               diag(tmp)=0
      #                                               tmp}), 
      #                                theta = target_graph_cov)
      #       
      #       ROC_inv = huge::huge.roc(path = lapply(ccl_cor, function(x){
      #                                               x = solve(x)
      #                                               tmp = (abs(x)>1e-11)*1
      #                                               diag(tmp)=0
      #                                               tmp}), 
      #                                theta = target_graph_inv)
      #       
      #       ROC = list(ROC_cov = ROC_cov, ROC_inv = ROC_inv)
      
      
      # it also has an implementation of SparCC, but seems not giving the same results
      #source("E:\\Dropbox\\Microbial_Networks\\codes\\CCLasso-master\\R/SparCC.R");
      #res_spa_count <- SparCC.count(x = data);
      #res$Cov[1:5, 1:5]
      #res_spa_count$cov.w[1:5, 1:5]
      
      time= sum(system.time({
        res_ccl_count <- cclasso(x = data_rep[[1,1]],
                                 counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                                 k_cv = 3, # cv folds, min=2
                                 lam_int = c(1e-6, 3), #tuning parameter value range; need to vary this for ROC curve
                                 k_max=20, n_boot =20)  # input format n by p
        
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
      source(paste0(filepath, "/Kun_code/COAT-master/COAT-master/simulation.R")) # this contains all different data generating models
      source(paste0(filepath, "/Kun_code/COAT-master/COAT-master/coat.R"))
      
      
      # this function runs to provide the ROC; it seems peudocounts need to be added in advance for zero counts, and composition matrix adjusted accordingly
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      # # look at this function under simulation.R, it seems to be implemented for coat ROC
      # Coat_ROC_res = calCoatROC(dataCell = data_rep, 
      #                           sigmaTrue = option$Sigma_list$Sigma,
      #                           precTrue = option$Sigma_list$Omega,
      #                           nPlotPoint = 31, nGrid = 30, soft = 1) 
      # 
      # # plot(y=Coat_ROC_res$tp,x=Coat_ROC_res$fp, type='o' ) # the output is averaged over different repetitions. Here to distribute the workload just input one repetition
      # 
      # library(DescTools)
      # ROC_cov = list(tp = Coat_ROC_res$tprGrid, fp = Coat_ROC_res$fprGrid)
      # ROC_cov$AUC = AUC(x=ROC_cov$fp , y=ROC_cov$tp)
      # 
      # ROC_inv = list(tp = Coat_ROC_res$tprGrid_inv, fp = Coat_ROC_res$fprGrid_inv)
      # ROC_inv$AUC = AUC(x=ROC_inv$fp , y=ROC_inv$tp)
      # 
      # ROC = list(ROC_cov = ROC, ROC_inv = ROC_inv)
      
      
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

      time = sum(system.time({
        # input non-normalized, count data, without adding pesudocounts
        Spiec_network <- spiec.easi(data_rep[[2, 1]], method='glasso', # choose from 'mb' for neighbourhood selection, or 'glasso'
                                    # if to perform model selection
                                    pulsar.select = T,
                                    pulsar.params = list(
                                      thresh=0.05,# Threshold for StARS criterion.
                                      subsample.ratio=0.8, # Subsample size for StARS.
                                      rep.num = 20), # Number of subsamples for StARS.
                                    lambda.min.ratio= min(lambda_seq)/max(lambda_seq), 
                                    lambda.max = max(lambda_seq), nlambda=length(lambda_seq)
                                    # , lambda = lambda_seq
        ) # lambda is for penalty parameter, is tuning parameter; override to set same as lambda_seq
        
        
      })[2:3])
      
      # Spiec_network$lambda
      # lambda_seq
      cat('done1')
      # the result should be a solution path over lambda; the path is adjacency matrix with diagonal being zero
      Spiec_ROC_res = huge::huge.roc(path = Spiec_network$est$path, 
                                     theta = target_graph_inv, verbose=FALSE)
      
      
      # solve(Spiec_network$est$icov[[1]])[1:10, 1:10]
      # Spiec_network$est$cov[[1]][1:10, 1:10]
      ROC_inv = Spiec_ROC_res
      cat('done2')
            
      # thresholding the optimal cov matrix for ROC
      Spiec_cov = as.matrix(getOptCov(Spiec_network))
      # print(class(Spiec_cov))
      ROC_cov <- ROC_cov_pathbased <- NULL
      
      if(F){
              ROC_cov = get_cov_ROC(Spiec_cov, target_graph_cov)
      
      cat('done3')
      # or use the cov path corresponding to the lambda path, and compute ROC. to eliminate very small values, make Cov into Cor and threhsold at 1e-11
      ROC_cov_pathbased = 
        huge::huge.roc(path = 
                         lapply(Spiec_network$est$cov, function(x){
                            x = diag(1/sqrt(diag(x)))%*% x %*% diag(1/sqrt(diag(x)))# transform into correlation matrix
                            diag(x) = 0
                            tmp = (abs(x)>1e-11)*1
                            tmp}), 
                       theta = target_graph_cov)

      
      cat('done4')

      }
      
      # also need to get fp at optimal tuning value for null model
      precision = getOptNet(Spiec_network)
      fp_inv_null = calTprFpr(sigmaHat = precision, sigmaTrue = matrix(0, p, p))$fp
      
      ROC = list(ROC_cov = ROC_cov, ROC_cov_pathbased=ROC_cov_pathbased,
                 ROC_inv = ROC_inv, fp_inv_null = fp_inv_null, time=time)
      cat('done5')
      # getOptMerge(Spiec_network) # symmetric matrix with edge-wise stability; used for getting 0/1 network
      # getOptInd(Spiec_network) # index of the selected lambda from provided lambda path
      # getOptiCov(Spiec_network) # the optimal inverse covariance matrix (glasso only)
      # getOptCov(Spiec_network) # the optimal covariance matrix associated with the selected network (glasso only)
      # getOptBeta(Spiec_network) # the optimal coefficient matrix (mb only) (should be corresponding to inv-cov level, but may not the same)
      # getOptNet(Spiec_network) # the optimal (StARS-refit) network, 0/1 adjacency matrix
      
    }else if(grepl(method,'gcoda',ignore.case=TRUE)){
      #--------------
      # gCoda
      #--------------
      source(paste0(filepath, '/Kun_code/gCoda-master/R/gcoda.R'))
      # input zero corrected compositional data
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      time = sum(system.time({
        gcoda_network <- gcoda(data_rep[[1,1]], counts = F, pseudo = 1, 
                             lambda.min.ratio=1e-3, # these will be ignored
                             nlambda= length(lambda_seq), # this has to be the length of provided lambda_seq
                             lambda = lambda_seq);  # it seems gcoda will automatically include more lambda values. need to pay attention to this.
      })[2:3])
      # gcoda_network$lambda
      # gcoda_network$opt.index  # if at the boundary may need to change lambda.min.ratio; however maximum cannot be changed... any reason how they choose the lambda.max?
      # gcoda_network$opt.icov
      
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
      source(paste0(filepath, '/Kun_code/Jing_lib/func_libs.R'))
      
      # supply uncorrected compositional data with n by p data matrix. It will be mclr transformed inside their function
      
      time = sum(system.time({
        fit.spring <- SPRING(data_rep[[1,1]], 
                           quantitative = F, # F means input is compositional
                           lambdaseq = lambda_seq, 
                          # nlambda = 20, 
                           ncores = 1,
                           subsample.ratio = 0.8, rep.num = 20 # these are for tuning parameter selection
      )
      })[2:3])
      
      
      spring_ROC_res = huge::huge.roc(fit.spring$fit$est$path, theta = option$Sigma_list$A_inv, verbose = F)
      ROC_inv = spring_ROC_res
      
      
      opt.K <- fit.spring$output$stars$opt.index
      adj.icov <- as.matrix(fit.spring$fit$est$path[[opt.K]]) # adjacency matrix
      fp_inv_null = calTprFpr(sigmaHat = adj.icov, sigmaTrue = matrix(0, p, p))$fp
      
      
      # threshold the optimal cov
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
      
      
      # # StARS-selected lambda index based on the threshold (default = 0.01)
      # opt.K <- fit.spring$output$stars$opt.index
      # # Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
      # adj.K <- as.matrix(fit.spring$fit$est$path[[opt.K]])
      # # Estimated partial correlation coefficient, same as negative precision matrix.
      # pcor.K <- as.matrix(SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode = 'maxabs'))
      
      
    }else{
      stop('no method matched')
    }
    
  }
  
  
  return(ROC)
}
