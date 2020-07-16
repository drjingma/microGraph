###-----------------------------------------------------------
###
### Implementation of available methods for computing microbiome network
###
###-----------------------------------------------------------

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
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DirichletMultinomial")

# library(devtools)
# install_github("hallucigenia-sparsa/seqgroup") 
library(seqgroup)
# reference: https://hallucigenia-sparsa.github.io/seqgroup/reference/barebonesCoNet.html
# input data needs p by n
barebonesCoNet(abundances = t(data_list$data), methods = c("spearman",'pearson', 'bray',"kld"), 
               method.num.T = 4, pval.T = 0.05,
               pval.cor = FALSE, # do not use usual correlation test
               permut = T, renorm = T,permutandboot = T, # their permutation/bootstrap and renormalization will remove spurious correlation
               iters = 100, bh = TRUE,
               pseudocount = 1e-11, # counts added to zeros when taking log
               plot = T, verbose = T) 


barebonesCoNet(abundances = t(data_list$data), methods = c("spearman",'pearson', 'bray',"kld"), 
               method.num.T = 4, pval.T = 0.05,
               pval.cor = FALSE, # do not use usual correlation test
               permut = F, renorm = F,permutandboot = F, # their permutation/bootstrap and renormalization will remove spurious correlation
               iters = 100, bh = TRUE,
               pseudocount = 1e-11, # counts added to zeros when taking log
               plot = T, verbose = T) 

#--------------
# SparCC 
#--------------
#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
res = sparcc(data, iter = 20, inner_iter = 10, th = 0.1) # need input format n by p
res$Cov # the estimated log counts covariance, the Sigma
res$Cor # the estimated log counts correlation, the corr_mat


#--------------
# CClasso
#--------------
## (obtained from Github: https://github.com/huayingfang/CCLasso)
source("E:\\Dropbox\\Microbial_Networks\\codes\\CCLasso-master\\R\\cclasso.R");

res_ccl_count <- cclasso(x = data, counts = T, pseudo = 0.5);  # input format n by p
res_ccl_count$cor_w
res_ccl_count$p_vals # not sure how this is computed

# it also has an implementation of SparCC, but seems not giving the same results
#source("E:\\Dropbox\\Microbial_Networks\\codes\\CCLasso-master\\R/SparCC.R");
#res_spa_count <- SparCC.count(x = data);
#res$Cov[1:5, 1:5]
#res_spa_count$cov.w[1:5, 1:5]



#--------------
# COAT 
#--------------
source('E:\\Dropbox\\Microbial_Networks\\codes\\COAT-master\\COAT-master\\simulation.R') # this contains all different data generating models
source('E:\\Dropbox\\Microbial_Networks\\codes\\COAT-master\\COAT-master\\coat.R')
coat(x, nFolder=5, soft=1) # x is n by p data matrix



#--------------
# SPIEC-EASI
#--------------
se <- spiec.easi(X, method='glasso', # choose from 'mb' for neighbourhood selection, or 'glasso'
                 pulsar.params = list(
                   thresh=0.05,# Threshold for StARS criterion.
                   subsample.ratio=0.8, # Subsample size for StARS.
                   rep.num = 20), # Number of subsamples for StARS.
                 lambda.min.ratio=1e-2, nlambda=5) # lambda is for penalty parameter, is tuning parameter
# the result should be a solution path over lambda. Select based on STARS functions
se$select
plot(huge::huge.roc(se$est$path, graph, verbose=FALSE))
stars.pr(getOptMerge(se), graph, verbose=FALSE)

getOptMerge(se) # symmetric matrix with edge-wise stability; used for getting 0/1 network
getOptInd(se) # index of the selected lambda from provided lambda path
getOptiCov(se) # the optimal inverse covariance matrix (glasso only)
getOptCov(se) # the optimal covariance matrix associated with the selected network (glasso only)
getOptBeta(se) # the optimal coefficient matrix (mb only) (should be corresponding to inv-cov level, but may not the same)
getOptNet(se) # the optimal (StARS-refit) network, 0/1 adjacency matrix


#--------------
# gCoda
#--------------
source('E:\\Dropbox\\Microbial_Networks\\codes\\gCoda-master\\R\\gcoda.R')
res_gcoda_count <- gcoda(x, counts = T, lambda.min.ratio=1e-3, nlambda=20);
res_gcoda_count$lambda
res_gcoda_count$opt.index  # if at the boundary may need to change lambda.min.ratio; however maximum cannot be changed... any reason how they choose the lambda.max?
res_gcoda_count$opt.icov

#---------------
# SPRING
#---------------
devtools::install_github("irinagain/mixedCCA")
devtools::install_github("GraceYoon/SPRING")
library(SPRING)
data("QMP") # load the data available from this package, containing 106 samples and 91 OTUs.

# Apply SPRING on QMP data.
fit.spring <- SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific", nlambda = 50, rep.num = 50)
fit.spring <- SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific", nlambda = 5, rep.num = 5)
# This takes around 23 minutes. We are working on reducing the computation time (10/25/2019).

# StARS-selected lambda index based on the threshold (default = 0.01)
opt.K <- fit.spring$output$stars$opt.index
# Estimated adjacency matrix from sparse graphical modeling technique ("mb" method) (1 = edge, 0 = no edge)
adj.K <- as.matrix(fit.spring$fit$est$path[[opt.K]])
# Estimated partial correlation coefficient, same as negative precision matrix.
pcor.K <- as.matrix(SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K]], mode = 'maxabs'))

