#######################
##
##
##  Data Examples for book chapter
##
##
########################
filepath = 'C:\\Users\\yuek\\Dropbox\\Microbial_Networks\\microGraph' #BOX


filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' # MAC

filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph'
setwd(filepath)
source('lib/func_libs.R')
source('Kun_code/generation_function.R')
source('Kun_code/method_function.R')
library(SpiecEasi)
library(MASS)

data(amgut1.filt)
reference_data = amgut1.filt
dim(reference_data) #data format n by p
#---------------
# part I: null hypothesis
#---------------

# parameters that need to discuss:
n = c(100, 200, 500)[1]
p = 200
reference_data = amgut1.filt
library_scale = rnegbin(n,mu=3e4, theta=9e8/(3e6-3e4)) # use library_scale to directly specify the total counts per sample

set.seed(102)

## Null 1
# shuffle reference data set, empty graph for both cov and inv-cov
Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data, Sigma_list = Sigma_list)

## Null 2
# generate from Dirichlet distribution, so will assume both inv and cov based graph being empty?
mu = runif(p, 0, 4) # the mu is for underlying log-normal of Dirichlet parameters
alpha=100
Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu, Sigma_list = Sigma_list)

## Alternative 1
# generate from random graph, p nodes and m edges; const for correcting non-positive definite Sigma
# this graph corresponds to inverse covariance
Sigma_list  = model_gnm(p=p, m=p, const = 0.1)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
Sigma_list$A_inv = Sigma_list$Adj
eps = 1e-11
Sigma_list$A_cov = (abs(Sigma_list$Sigma)>eps)*1; diag(Sigma_list$A_cov) = 0
mu = runif(p, 0, 4)
option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list)

## Alternative 2
# use Copula generative model, graph is based on inverse covariance matrix.
Sigma_list = SpiecEasi_graph_Sigma(data = reference_data, type = c('band', 'cluster', 'erdos_renyi', 'hub', 'scale_free', 'block')[3])
option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, Sigma_list = Sigma_list)


data_list = data_generation(n, p, option)
dim(data_list$data)


data = data_list$data
colnames(data) <- paste('taxa', 1:ncol(data))
rownames(data) <- paste('cell', 1:nrow(data))

# wrap nrep replications for running the methods; data_rep[[1,k]] is composition matrix; data_rep[[2,k]] is count matrix
# the data matrix is not zero-corrected. 
nreps=2
data_rep = matrix(list(), 2, nreps)
for(k in 1:nreps){
  set.seed(k)
  data_list = data_generation(n, p, option)
  dim(data_list$data)
  
  data = data_list$data
  colnames(data) <- paste('taxa', 1:ncol(data))
  rownames(data) <- paste('cell', 1:nrow(data))
  
  data_rep[[1, k]] = count_to_comp(data) # compositions
  data_rep[[2, k]] = data # counts
}

n = nrow(data_rep[[1,1]])
p = ncol(data_rep[[1,1]])

save(data_rep, file = 'data_rep.RData')

# compare methods for their ROC
# input the collection of data matrixs, data[[1,k]] is composition, data[[2,k]] is count matrix, for total of nrep repetitions
cov_conet_roc = compare_methods(data_rep[,1, drop=F], 
                target = c('covariance', 'precision')[1],
                method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                           'SpiecEasi', 'gCoDa','Spring')[1])

cov_sparcc_roc = compare_methods(data_rep[,1, drop=F],  
                             target = c('covariance', 'precision')[1],
                             method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                        'SpiecEasi', 'gCoDa','Spring')[2])

cov_cclasso_roc = compare_methods(data_rep[,1, drop=F], 
                              target = c('covariance', 'precision')[1],
                              method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                         'SpiecEasi', 'gCoDa','Spring')[3])

cov_coat_roc = compare_methods(data_rep[,1, drop=F],  
                           target = c('covariance', 'precision')[1],
                           method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                      'SpiecEasi', 'gCoDa','Spring')[4])






