## simulation study for compositional network learning
library(entropy)
library(gtools)
library(SpiecEasi)
library(seqgroup)

data(amgut1.filt)
p = 127
n = 289
set.seed(1)
library_scale = 20000


###------------------------
### Data generating model
###------------------------

##-------------------
## Null model: Type I error evaluation
##-------------------


#---------------
# option one: generate each OUT independently from marginal distributions, with parameter values referencing some empirical data
# can shuffle an available dataset
#---------------
null_data_generate_1 = function(x){ # x is n by p real data matrix
  # shuffle within taxa, so keep its marginal distribution
  n = nrow(x)
  p = ncol(x)
  cat('Shuffle reference dataset, nsub ', n, ', ntaxa ', p, '\n')
  index = replicate(p, sample(n, replace=F))
  x_shuffle = sapply(1:p, function(i)x[index[,i], i])
  return(x_shuffle)
}

#---------------
# option two: generate composition from Dirichlet distribution, then generate counts from multinomial. Equivalent to independent Gamma distributions
#---------------
alpha = runif(p) 
null_data_generate_2 = function(n, p, library_scale, alpha){ # generate total counts from Poisson(library_scale), generate Dirichelt(alpha) with a fixed alpha
  prob_vec = rdirichlet(n=n, alpha = alpha) 
  data = t(sapply(1:n, function(x)rmultinom(n=1, size=rpois(n=1,library_scale), prob=prob_vec[x,]))) 
  return(data)
}

#---------------
# option three: based on option two model, subtract all entries by the data matrix mean, set negative values to zero -> increase sparsity level 
#---------------
null_data_generate_3 = function(n, p, library_scale){
  data = null_data_generate_2(n, p, library_scale)
  data = data-mean(data)
  data[data<0]<-0
  return(data)
}

##### Parametric network model: discovery rate evaluation
##-----------------
## option one: logistic normal multinomial:
##-----------------
para_data_generate_1 = function(n, p, library_scale){
  mu_vec = rnorm(p-1, sd=1) # mu and Sigma need some reference
  Sigma_mat_decomp = matrix(rnorm((p-1)^2, sd=0.1), p-1, p-1)
  Sigma_mat = Sigma_mat_decomp%*%t(Sigma_mat_decomp)
  z = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*(p-1), mean = 0), n, p-1, byrow=T)) + mu_vec
  z = cbind(t(z), rep(0, n))
  x = exp(z)/rowSums(exp(z))
  data = t(sapply(1:n, function(i)rmultinom(n=1, size=rpois(n=1,library_scale), prob=x[i,])))
  return(data)
}
##-----------------
## option two: log normal distribution for generating abudance data; takes in mu and Sigma
##-----------------
para_data_generate_2 = function(n, p, mu, Sigma){
  Sigma_mat_decomp = t(chol(Sigma))
  # Sigma_mat = Sigma_mat_decomp%*%t(Sigma_mat_decomp)
  x = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*p, mean = 0), n, p, byrow=T)) + mu
  y = exp(x)
  data=y  
  return(data)
}

##-----------------
## Conet naive generative model
##-----------------
CoNet_data_generate = function(n, p, library_scale){
  noise = sample(c(-1, 1), replace=T, size=n*p)*rpois(n*p, library_scale/p/20)
  t1 = round(seq(library_scale/4/3/3, library_scale/4/3*3, length=n))
  t2 = rep(library_scale/4/2, n)
  t3 = round(seq(library_scale/4/3*3, library_scale/4/3/3, length=n))
  t4 = rep(library_scale/4/2,n)
  t_rest = round(library_scale/p/(p-4))
  data = cbind(t1, t2, t3, t4, matrix(t_rest, nrow=n, ncol=p-4)) + matrix(noise, nrow=n, ncol=p)
  data[data<0] = 0
  return(data)
}


##-----------------
## SparCC generative model (parametric, log-normal) (simulate aubdance, not necesarily integer counts)
##-----------------
library(entropy)

SparCC_data_generate = function(n, p, mu0, mu1, prob){
  mu = c(mu1, rep(mu0, p-1))
  
  Sigma = matrix(0, p, p)
  rand_corr = as.numeric(runif(p*(p-1)/2)<prob)
  rand_corr[rand_corr!=0] = sample(c(-1, 1), replace=T, size=length(rand_corr[rand_corr!=0]))*rand_corr[rand_corr!=0]
  Sigma[lower.tri(Sigma, diag=F)] <- rand_corr*0.1*0.1 # the sd of each taxon is set to be 0.1
  Sigma = Sigma+ t(Sigma)
  diag(Sigma)<-0.01
  # make sure it is psd
  tmp = eigen(Sigma)
  eigenvalue = tmp$values
  eigenvalue[eigenvalue<=0]<-.Machine$double.eps
  Sigma = tmp$vector %*% diag(eigenvalue)%*%t(tmp$vectors)
  corr_mat = diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
  
  data = exp(tmp$vector%*%diag(sqrt(eigenvalue))%*%matrix(rnorm(n*p),nrow=p, ncol=n) + mu)
  n_eff = mean(sapply(1:n, function(i)entropy.ChaoShen(data[,i])))
  return(list(n_eff=n_eff, prob=prob, mu=mu, Sigma=Sigma, corr_mat = corr_mat, data=data))
}

x=seq(-5, 2, by=0.1)
nn = sapply(x, function(s)SparCC_data_generate(n, p, mu0=2, mu1=s, prob=0.1)$n_eff)
plot(x, nn, type='l')
nn = sapply(x, function(s)SparCC_data_generate(n, p, mu0=-1, mu1=s, prob=0.1)$n_eff) # with mu0 smaller and mu1 smaller, generally can have higher entropy...?
plot(x, nn, type='l', col=2)


#-----------------------
# CClasso generative model (parametric, log-normal for abundance)
#-----------------------
CClasso_mu = function(n,p){
  runif(n=p, -0.5, 0.5)
}
#option: random, neighbor, AR, hub, block
CClasso_Sigma = function(n, p, option){
  Sigma = matrix(0, p, p)
  if(option=='random'){
    Sigma[lower.tri(Sigma, diag=F)] <- sample(c(-0.15, 0.15), replace=T, size=p*(p-1)/2) * as.numeric(runif(p*(p-1)/2)<0.3)
    Sigma = Sigma+t(Sigma) + diag(1, p)
  }
  
  if(option=='neighbor'){
    library(FNN)
    pos = matrix(runif(p*2), p, 2)
    for (i in 1:p){
      neighbor=((1:p)[-i])[knnx.index(pos[-i,], pos[i,,drop=F], k=10)] # gives index of 10 nearest neighbours
      Sigma[i, neighbor] = 0.5
    }
    
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  if(option == 'AR'){
    library(stats)
    # create a toeplitz matrix
    Sigma = toeplitz(c(1, 0.4, 0.2, 0.2, 0.1, rep(0, p-5)))
  }
  
  if(option == 'hub'){
    hub = sample(1:p, size=3, replace = F)
    Sigma[hub, 1:p] = 0.2 * as.numeric(runif(3*p)<0.7)
    Sigma[(1:p)[-hub], 1:p] = 0.2 * as.numeric(runif((p-3)*p)<0.2)
    Sigma[lower.tri(Sigma, diag=T)] <-0
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  if(option=='block'){ # assume we just have p that can be divided by 5
    #library(caret)
    #flds <- createFolds(1:p, k = 5, list = T, returnTrain = FALSE)
    #flds = sapply(flds, sort)
    groupsize = p/5
    flds = lapply(1:5, function(i)((i-1)*groupsize+1):(i*groupsize))
    for(i in 1:5){
      Sigma[as.matrix(expand.grid(flds[[i]], flds[[i]]))] <- 0.4*as.numeric(runif(length(flds[[i]])^2)<0.6)
      Sigma[as.matrix(expand.grid(flds[[i]], (1:p)[-flds[[i]]] ))] <- 0.2*as.numeric(runif(length(flds[[i]])*(p-length(flds[[i]]))   )<0.2)
    }
    Sigma[lower.tri(Sigma, diag=T)] <-0
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  # make Sigma psd
  if(!all(eigen(Sigma)$value>0))
    stop('error, not PSD, not sure how to correct it for now')
  
  return(Sigma)
}


#-----------------------
# COAT generative model (parametric, log-normal / log-gamma)
#-----------------------
source('E:\\Dropbox\\Microbial_Networks\\codes\\COAT-master\\COAT-master\\simulation.R') # this contains all different data generating models
data_list = generateDataCell(n, modelCov, p1 = 50, p2 = 100, p3 = 200, nRep = 100) #modelCov = 1: hub cov; modelCov = 2: block cov; modelCov = 3: sparse cov


#----------------------
# SpiecEasi generative model (from real data)
#----------------------
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)


data(amgut1.filt)
range(rowSums(amgut1.filt))

SpiecEasi_graph = function(amgut1.filt, type='erdos_renyi'){
  depths <- rowSums(amgut1.filt) # raw counts data, rows are obs, n by p
  amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
  amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
  
  d <- ncol(amgut1.filt.cs)
  n <- nrow(amgut1.filt.cs)
  e <- d
  

  graph <- make_graph(type, d, e) # choose from band, cluster, erdos_renyi, hub, scale_free, block; can also write my own graph (essentially an adjacency matrix)
  return(graph)
}

SpiecEasi_generate = function(amgut1.filt, graph){
  depths <- rowSums(amgut1.filt) # raw counts data, rows are obs, n by p
  amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
  amgut1.filt.cs <- round(amgut1.filt.n * min(depths))
  
  d <- ncol(amgut1.filt.cs)
  n <- nrow(amgut1.filt.cs)
  e <- d
  
Prec  <- graph2prec(graph, posThetaLims = c(2, 3), # this is precision strength range
                    targetCondition = 100, # this is condition number kappa
                    epsBin = 0.01, numBinSearch = 100) # here can set different generating parameters
Cor   <- cov2cor(prec2cov(Prec))

X <- synth_comm_from_counts(amgut1.filt.cs, mar=2, 
                            distr='zinegbin', # choose from zipois, zinegbin, nbinom, pois, lnorm, negbin
                            Sigma=Cor, n=n) # output data n by p

return(list(data=X, graph=graph, Correlation = Cor, Precision=Prec ))
}



#####################
data = null_data_generate_2(n, p, library_scale)
data = para_data_generate_1(n, p, library_scale)
data = CoNet_data_generate(n, p, library_scale)
data_list = SparCC_data_generate(n, p, mu0=-1, mu1=-3, prob=0.1); data = t(data_list$data)
data = para_data_generate_2(n, p, mu=CClasso_mu(n,p), Sigma = CClasso_Sigma(n, p, 'random'))



dim(data) # n by p
data = data.frame(t(data))
rownames(data) = paste0('taxa', 1:p)
write.table(data, file='E:\\Dropbox\\Microbial_Networks\\codes\\data\\data.txt', sep='\t', col.names = FALSE, quote=F)


###-----------------------------------------------------------
###
### Implementation of available methods
###
###-----------------------------------------------------------

#--------------
# CoNet 
#--------------
## (input data should be p by n)
## option 1: use the CoNet software (not in R)
filepath_conet = 'E:\\Dropbox\\Microbial_Networks\\codes\\data\\Conet_output2.txt'
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

## another option: use its R version without most of the features
# install.packages('gsl')
# library(gsl)
# file.path(R.home("gsl"), "Makeconf")
# 
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DirichletMultinomial")


library(devtools)
install_github("hallucigenia-sparsa/seqgroup") 
library(seqgroup)
# reference: https://hallucigenia-sparsa.github.io/seqgroup/reference/barebonesCoNet.html
barebonesCoNet(abundances, metadata = NULL, methods = c("spearman",
                                                        "kld"), T.up = NA, T.down = NA, method.num.T = 2, pval.T = 0.05,
               init.edge.num = max(2, round(sqrt(nrow(abundances)))), min.occ = 0,
               keep.filtered = TRUE, norm = FALSE, stand.rows = FALSE,
               pval.cor = FALSE, permut = FALSE, renorm = FALSE,
               permutandboot = FALSE, iters = 100, bh = TRUE,
               pseudocount = 1e-11, plot = FALSE, verbose = FALSE) 

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




