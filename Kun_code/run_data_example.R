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
option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data)

## Null 2
mu = runif(p, 0, 4) # the mu is for underlying log-normal of Dirichlet parameters
alpha=100
option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu)

## Alternative 1
# generate from random graph, p nodes and m edges; const for correcting non-positive definite Sigma
Sigma_list  = model_gnm(p=p, m=p, const = 0.1)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
mu = runif(p, 0, 4)
option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list)

## Alternative 2
graph_Sigma = SpiecEasi_graph_Sigma(data = reference_data, type = c('band', 'cluster', 'erdos_renyi', 'hub', 'scale_free', 'block')[3])
option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, graph_Sigma = graph_Sigma)


data_list = data_generation(n, p, option)
dim(data_list$data)

data = data_list$data
colnames(data) <- paste('taxa', 1:ncol(data))
rownames(data) <- paste('cell', 1:nrow(data))

# wrap for COAT method
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




