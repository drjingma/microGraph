#######################
##
##
##  Data Examples for book chapter
##
##
########################
 
# filepath = 'C:\\Users\\yuek\\Dropbox\\Microbial_Networks\\microGraph' #BOX
# 
# filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' # MAC
# 
# filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC

filepath = '~/Desktop/micro_net' #bayes


setwd(filepath)
source('lib/func_libs.R')
source('Kun_code/generation_function.R')
source('Kun_code/method_function.R')
library(SpiecEasi)
library(MASS)

data(amgut1.filt)
reference_data = amgut1.filt
# dim(reference_data) #data format n by p
# n = c(100, 200, 500)[1]
# p = 200

args = commandArgs(trailingOnly = T) # (n p null1 2),  n and p override by reference data set, e.g. 
# args = c(50, 50, 'alt1', 2)
# args = c(20, 30, 'null2', 2)
# args = c(100, 127, 'null1.1', 100)
n = as.integer(args[1])
p = as.integer(args[2])
choose_model = as.character(args[3])
nreps = as.integer(args[4])
distr = ifelse(is.na(args[5]), NA, as.character(args[5]))
if(is.na(distr)) distr = NULL

set.seed(102)

#---------------
# part I: null hypothesis
#---------------

if(choose_model == 'null1'){
  ## Null 1
  # shuffle reference data set, empty graph for both cov and inv-cov
  reference_data = amgut1.filt[sample(1:nrow(amgut1.filt), size=n, replace = F), ]
  p = ncol(reference_data)
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data, Sigma_list = Sigma_list)
}

if(choose_model == 'null1.1'){
  ## Null1 1.1
  # Use reference data set, marginal distributions are NB, and the network is set empty
  reference_data = amgut1.filt[sample(1:nrow(amgut1.filt), size=n, replace = F), ]
  p = ncol(reference_data)
  # generate a null graph object
  option = list(hypothesis = 'null', model = 'copula', reference_data = reference_data, Sigma_list = NULL, distr = distr)
  option$Sigma_list$A_inv <-  option$Sigma_list$A_cov <- matrix(0, p, p)
  option$Sigma_list$Sigma <- option$Sigma_list$Omega <- diag(1, p)
}

if(choose_model == 'null2'){
  ## Null 2
  # generate from Dirichlet distribution, so will assume both inv and cov based graph being empty?
  mu = runif(p, 0, 4) # the mu is for underlying log-normal of Dirichlet parameters
  alpha=100
  library_scale = rnegbin(n,mu=3e4, theta=9e8/(3e6-3e4)) # use library_scale to directly specify the total counts per sample
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu, Sigma_list = Sigma_list)
}



#---------------
# part II: Alternative hypothesis
#---------------


if(choose_model=='alt1'){
  ## Alternative 1
  # use Copula generative model, graph is based on inverse covariance matrix.
  reference_data = amgut1.filt[1:n, ]
  p = ncol(reference_data)
  Sigma_list = SpiecEasi_graph_Sigma(data = reference_data, type = c('band', 'cluster', 'erdos_renyi', 'hub', 'scale_free', 'block')[3])
  option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, Sigma_list = Sigma_list, distr = distr)
}


if(choose_model == 'alt2'){
  ## Alternative 2
  # generate from random graph, p nodes and m edges; const for correcting non-positive definite Sigma
  # this graph corresponds to inverse covariance
  Sigma_list  = model_gnm(p=p, m=p, const = 0.1)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  Sigma_list$A_inv = Sigma_list$Adj
  eps = 1e-11
  Sigma_list$A_cov = (abs(Sigma_list$Sigma)>eps)*1; diag(Sigma_list$A_cov) = 0
  mu = runif(p, 0, 4)
  option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list)
  
}

# data_list = data_generation(n, p, option)
# dim(data_list$data)
# 
# 
# data = data_list$data
# colnames(data) <- paste('taxa', 1:ncol(data))
# rownames(data) <- paste('cell', 1:nrow(data))


## wrap nrep replications for running the methods; data_rep[[1,k]] is composition matrix; data_rep[[2,k]] is count matrix
## the data matrix is not zero-corrected. 

data_rep = matrix(list(), 2, nreps)
for(k in 1:nreps){
  set.seed(k*100)
  data_list = data_generation(n, p, option)
  dim(data_list$data)
  
  data = data_list$data
  colnames(data) <- paste('taxa', 1:ncol(data))
  rownames(data) <- paste('cell', 1:nrow(data))
  
  data_rep[[1, k]] = count_to_comp(data) # compositions
  data_rep[[2, k]] = data # counts
}

save(data_rep, file = paste0('dist_data/', option$distr, '/n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))
save(list=c('data_rep', 'option', 'n', 'p'), file = paste0('dist_data/', option$distr, '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))






