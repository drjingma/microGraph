###-----------------------------------------------------------
###
### Implementation of available methods for computing microbiome network
###
###-----------------------------------------------------------
filepath = 'C:\\Users\\yuek\\Dropbox\\Microbial_Networks\\microGraph\\Kun_code' #BOX

filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph/Kun_code'
filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph\\Kun_code'

count_to_comp = function(W){
  X = sweep(W,1,STATS = rowSums(W), FUN='/') # the compositions from log-normal; different Xi for different cells
  return(X)
}

#--------------
# CoNet 
#--------------
## (input data should be p by n)

### option 1: use the CoNet software (not in R)

filepath_conet = 'CoNet_app_output\\Conet_output2.txt'
read_tab = function(filepath){
  
  processFile = function(filepath) {
    full = list()
    count=1
    con = file(filepath, "r")
    while ( TRUE ) {
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      full[count]=line
      count=count+1
    }
    close(con)
    return(full)
  }
  
  see = processFile(filepath)
  
  # search for nodes information
  i=1
  while (i <=length(see)){
    if(substr(see[[i]],1,6) == ';NODES'){nodehead = i}
    if(substr(see[[i]],1,5) == ';ARCS'){edgehead=i}
    i=i+1
  }
  node_info = read.table(filepath, skip=nodehead-1, nrow=edgehead-nodehead, sep='\t')
  edge_info = read.table(filepath, skip=edgehead-1, nrow=length(see) - edgehead+1, sep='\t', header=T)
  colnames(edge_info)[1]<-'Node2'
  edge_info = cbind(Node1 = rownames(edge_info), edge_info)
  rownames(edge_info)<-NULL
  
  return(list(node_info = node_info, edge_info = edge_info))
}
# read_tab(filepath_conet) # still need toi determine how to use the edge informations.

### option 2: use its R version without some features
# install.packages('gsl')
# library(gsl)
# file.path(R.home("gsl"), "Makeconf")
# 
# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DirichletMultinomial")

library(devtools)
install_github("hallucigenia-sparsa/seqgroup")
library(seqgroup)
# reference: https://hallucigenia-sparsa.github.io/seqgroup/reference/barebonesCoNet.html
# input data needs p by n; need to name the taxa
net_CoNet = barebonesCoNet(abundances = t(data), 
                           methods = c("spearman",'pearson', 'bray',"kld"), 
               method.num.T = 4, pval.T = 0.05,
               #init.edge.num is using default value (sqrt(p)); this is the number of top and bottom edges to initially keep for later testing
               init.edge.num = 50,
               pval.cor = FALSE, # do not use usual correlation test
               permut = T, renorm = T,permutandboot = T, # their permutation/bootstrap and renormalization will remove spurious correlation
               iters = 100, bh = TRUE,
               pseudocount = 1e-11, # counts added to zeros when taking log
               plot = F, verbose = F) 


# the output is am igraph object with edge weights; 
# here consider tuning parameter init.edge.num ( methods, method.num are just set fixed, though can also be changed)
# init.
edge_attr(net_CoNet)
network_matrix_CoNet = as_adjacency_matrix(net_CoNet,type='both', attr = "weight", sparse=F)

net_CoNet = barebonesCoNet(abundances = t(data), methods = c("spearman"), 
               method.num.T = 2, pval.T = 0.05,
               pval.cor = F, # do not use usual correlation test
               init.edge.num = 200,
               
               permut = F, renorm = T,permutandboot = T, # their permutation/bootstrap and renormalization will remove spurious correlation
               iters = 100, bh = TRUE,
               pseudocount = 1e-11, # counts added to zeros when taking log
               plot = T, verbose = F) 
plot(net_CoNet)









compare_methods = function(data_rep, # the collection of data matrixs, data[[1,k]] is composition, data[[2,k]] is count matrix, run for one repetition at a time
                           target = c('covariance', 'precision'),
                           method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                      'SpiecEasi', 'gCoDa','Spring'),
                           target_graph_cov, target_graph_inv # input the targeted graph; if provide var/inverse cov, must set the diagonal to zero first
         
){
  
  diag(target_graph) <-0
  
  if (grepl(target,'covariance',ignore.case=TRUE)){
    if ( sum(grepl(method,c('CoNet', 'SparCC', 'CCLasso', 'COAT'),ignore.case=TRUE))==0) stop('target is Covariance graph, specified method not targeting the correct graph')
    
    if(grepl(method,'CoNet',ignore.case=TRUE)){
      ...
      
      
      
      
      
      
      
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
      
      sparcc_bootstrap_res = sparccboot(data = data_rep[[2,1]], 
                                        R=200,  # number of bootstraps
                                        ncpus = 4, 
                                        sparcc.params = list(iter = 20, inner_iter = 10, th = 0.1)) # this is very slow
      
      pval_cov =  pval.sparccboot(sparcc_bootstrap_res)
      sparcc_bootstrap_pval <- sparcc_bootstrap_cor <- matrix(0, p, p)
      sparcc_bootstrap_pval[upper.tri(sparcc_bootstrap_pval)] = pval_cov$pvals
      sparcc_bootstrap_pval = sparcc_bootstrap_pval+t(sparcc_bootstrap_pval) 
      diag(sparcc_bootstrap_pval) = 2 # the diagonal set to 2, so never detect an edge on the diagonal
      
      sparcc_bootstrap_cor[upper.tri(sparcc_bootstrap_pval)] = pval_cov$cors
      sparcc_bootstrap_cor = sparcc_bootstrap_cor+t(sparcc_bootstrap_cor)
      diag(sparcc_bootstrap_cor) = 1
      
      sparcc_bootstrap_pval_FDR = matrix(p.adjust(p = as.vector(sparcc_bootstrap_pval), method = 'fdr'), p, p, byrow=F)
      
      # will changing the threshold for the matrix by hand to construct ROC
      nlam = 20
      get_path_raw <- get_path_FDR <-  list()
      ind = 0
      for (alpha in seq(0, 1, length = nlam)){
        ind = ind+1
        get_path_raw[[ind]] = (sparcc_bootstrap_pval<=alpha)*1
        get_path_FDR[[ind]] = (sparcc_bootstrap_pval_FDR<=alpha)*1
        
      }
      
      ROC_raw = huge::huge.roc(path = get_path_raw, theta = option$Sigma_list$A_cov)
      ROC_FDR = huge::huge.roc(path = get_path_FDR, theta = option$Sigma_list$A_cov)
      
      ROC = list(ROC_raw = ROC_raw, ROC_FDR = ROC_FDR)

    }else if(grepl(method,'CCLasso',ignore.case=TRUE)){
      #--------------
      # CClasso
      #--------------
      ## (obtained from Github: https://github.com/huayingfang/CCLasso)
      source(paste0(filepath, "\\CCLasso-master\\R\\cclasso.R"));
      
      lambda_vec = c(seq(0.01,0.1,0.01),seq(0.1,3,0.1)) * sqrt(log(p)/n) # Jing's setting for specifying the lambda list
      
      # for CCLasso, the input need to be zero-corrected. so add peudocount 1 if there are zero counts
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      lambda_vec = sort(lambda_vec)
      ccl_cor_path = list()
      ind = 0
      for(lambda in lambda_vec){
        ind=ind+1
        res_ccl_count <- cclasso(x = data_rep[[1,1]], 
                                 counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                                 k_cv = 2, # cv folds, min=2
                                 lam_int = rep(lambda), #tuning parameter value range; need to vary this for ROC curve
                                 k_max=200, n_boot =1);  # input format n by p
        
        # res_ccl_count$info_cv$lams # the lambda path 
        ccl_cor_path[[ind]] = res_ccl_count$cor_w
        # res_ccl_count$p_vals # from bootstrap, ignored in our simulation
      }
      
      ROC = huge::huge.roc(path = ccl_cor_path, theta = option$Sigma_list$A_cov)

      

      
      # it also has an implementation of SparCC, but seems not giving the same results
      #source("E:\\Dropbox\\Microbial_Networks\\codes\\CCLasso-master\\R/SparCC.R");
      #res_spa_count <- SparCC.count(x = data);
      #res$Cov[1:5, 1:5]
      #res_spa_count$cov.w[1:5, 1:5]
      
      
    }else if (grepl(method,'COAT',ignore.case=TRUE)){
      #--------------
      # COAT 
      #--------------
      source(paste0(filepath, "\\COAT-master\\COAT-master\\simulation.R")) # this contains all different data generating models
      source(paste0(filepath, "\\COAT-master\\COAT-master\\coat.R"))
      
      # x = data_rep[[2,1]] # counts
      # if(any(x==0)) x = x+1
      # x = sweep(data+1,1,STATS = rowSums(data+1), FUN='/')
      # coat(x, nFolder=5, soft=1) # x is n by p data matrix, need to be compositional and zero adjusted
      
      # this function runs to provide the ROC; it seems peudocounts need to be added in advance for zero counts, and composition matrix adjusted accordingly
      for( k in 1:ncol(data_rep)){
        if (any(data_rep[[2,k]]==0)){
          data_rep[[2,k]] = data_rep[[2,k]]+1
          data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
        }
      }
      
      # look at this function under simulation.R, it seems to be implemented for coat ROC
      Coat_ROC_res = calCoatROC(dataCell = data_rep, sigmaTrue = option$Sigma_list$Sigma,
                                nPlotPoint = 31, nGrid = 30, soft = 1) 
      
      # plot(y=Coat_ROC_res$tp,x=Coat_ROC_res$fp, type='o' ) # the output is averaged over different repetitions. Here to distribute the workload just input one repetition
      
      library(DescTools)
      ROC = list(tp = Coat_ROC_res$tprGrid, fp = Coat_ROC_res$fprGrid)
      ROC$AUC = AUC(x=ROC$fp , y=ROC$tp)

      
    }else{
      stop('no method matched')
    }
    
    
    
    
    
  }else if(grepl(target,'precision',ignore.case=TRUE)){
    if(grepl(method,'SpiecEasi',ignore.case=TRUE)){
      #--------------
      # SPIEC-EASI
      #--------------
      # for inv-covariance graph recovery
      library(pulsar)
      # input count data
      Spiec_network <- spiec.easi(data_rep[[2, 1]], method='glasso', # choose from 'mb' for neighbourhood selection, or 'glasso'
                                  # if to perform model selection
                                  pulsar.select = F,
                                  pulsar.params = list(
                                    thresh=0.05,# Threshold for StARS criterion.
                                    subsample.ratio=0.8, # Subsample size for StARS.
                                    rep.num = 20), # Number of subsamples for StARS.
                                  lambda.min.ratio=1e-2, nlambda=20) # lambda is for penalty parameter, is tuning parameter
      
      # the result should be a solution path over lambda; the path is adjacency matrix with diagonal being zero
      Spiec_ROC_res = huge::huge.roc(Spiec_network$est$path, option$Sigma_list$A_inv, verbose=FALSE)

      ROC = Spiec_ROC_res
      
      # getOptMerge(Spiec_network) # symmetric matrix with edge-wise stability; used for getting 0/1 network
      # getOptInd(Spiec_network) # index of the selected lambda from provided lambda path
      # getOptiCov(Spiec_network) # the optimal inverse covariance matrix (glasso only)
      # getOptCov(Spiec_network) # the optimal covariance matrix associated with the selected network (glasso only)
      # getOptBeta(Spiec_network) # the optimal coefficient matrix (mb only) (should be corresponding to inv-cov level, but may not the same)
      # getOptNet(Spiec_network) # the optimal (StARS-refit) network, 0/1 adjacency matrix
      
    }else if(grepl(method,'SpiecEasi',ignore.case=TRUE)){
      #--------------
      # gCoda
      #--------------
      source(paste0(filepath, '\\gCoda-master\\R\\gcoda.R'))
      gcoda_network <- gcoda(data, counts = T, pseudo = 1, lambda.min.ratio=1e-3, nlambda=20);
      
      gcoda_network$lambda
      gcoda_network$opt.index  # if at the boundary may need to change lambda.min.ratio; however maximum cannot be changed... any reason how they choose the lambda.max?
      gcoda_network$opt.icov
      
      gcode_ROC_res = huge::huge.roc(path=gcoda_network$path, theta = option$graph_Sigma$A, verbose=F)
      
    }else if(grepl(method,'Spring',ignore.case=TRUE)){
      #---------------
      # SPRING
      #---------------
      # devtools::install_github("irinagain/mixedCCA")
      # devtools::install_github("GraceYoon/SPRING")
      library(SPRING)
      source(paste0(filepath, '\\Jing_lib\\func_libs.R'))
      # supply directly the compositional data with n by p data matrix. It will be mclr transformed inside their function
      fit.spring <- SPRING(data_rep[[1,1]], quantitative = F, 
                           lambdaseq = "data-specific", 
                           nlambda = 20, 
                           subsample.ratio = 0.1, rep.num = 1 # these are for tuning parameter selection, for faster result we can set these small
      )
      
      
      spring_ROC_res = huge::huge.roc(fit.spring$fit$est$path, theta = option$Sigma_list$A, verbose = F)
      
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
