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


args = commandArgs(trailingOnly = T) # (n p null1 2 1),  n and p override by reference data set, e.g. args = c(50, 50, 'alt1', 2, 1)
n = as.integer(args[1])
p = as.integer(args[2])
choose_model = as.character(args[3])
nreps = as.integer(args[4])
run_rep = as.integer(args[5])


load(paste0('data/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))


# use gcoda implemented lambda path computation, based on the first rep data
x = data_rep[[2,1]]
if (any(x==0)){
  x = x+1
}
x = sweep(x,1,STATS = rowSums(x), FUN='/')

# Log transformation for compositional data
S <- var(log(x) - rowMeans(log(x)));
# Generate lambda via lambda.min.ratio and nlambda
lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
lambda.min.ratio = 1e-4
lambda.min <- lambda.min.ratio * lambda.max;


lambda_seq = exp(c(seq(log(lambda.min*0.9),log(lambda.max*1.1),length = 40)))
lambda_seq = sort(lambda_seq, decreasing = T)


# compare methods for their ROC

cov_conet_roc = compare_methods(data_rep[,run_rep, drop=F],
                                est_mat = c('covariance', 'precision')[1],
                                method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                           'SpiecEasi', 'gCoDa','Spring')[1],
                                target_graph_cov = option$Sigma_list$A_cov,
                                target_graph_inv = option$Sigma_list$A_inv,
                                option = option,
                                lambda_seq = lambda_seq)

cov_sparcc_roc = compare_methods(data_rep[,run_rep, drop=F],  
                                 est_mat = c('covariance', 'precision')[1],
                                 method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                            'SpiecEasi', 'gCoDa','Spring')[2],
                                 target_graph_cov = option$Sigma_list$A_cov,
                                 target_graph_inv = option$Sigma_list$A_inv,
                                 option = option,
                                 lambda_seq = lambda_seq)

cov_cclasso_roc = compare_methods(data_rep[,run_rep, drop=F], 
                                  est_mat = c('covariance', 'precision')[1],
                                  method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                             'SpiecEasi', 'gCoDa','Spring')[3],
                                  target_graph_cov = option$Sigma_list$A_cov,
                                  target_graph_inv = option$Sigma_list$A_inv,
                                  option = option)

cov_coat_roc = compare_methods(data_rep[,run_rep, drop=F],  
                               est_mat = c('covariance', 'precision')[1],
                               method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                          'SpiecEasi', 'gCoDa','Spring')[4],
                               target_graph_cov = option$Sigma_list$A_cov,
                               target_graph_inv = option$Sigma_list$A_inv,
                               option = option)

cov_spieceasi_roc = compare_methods(data_rep[,run_rep, drop=F],  
                               est_mat = c('covariance', 'precision')[2],
                               method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                          'SpiecEasi', 'gCoDa','Spring')[5],
                               target_graph_cov = option$Sigma_list$A_cov,
                               target_graph_inv = option$Sigma_list$A_inv,
                               option = option)

cov_gcoda_roc = compare_methods(data_rep[,run_rep, drop=F],  
                                    est_mat = c('covariance', 'precision')[2],
                                    method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                               'SpiecEasi', 'gCoDa','Spring')[6],
                                    target_graph_cov = option$Sigma_list$A_cov,
                                    target_graph_inv = option$Sigma_list$A_inv,
                                    option = option)

cov_spring_roc = compare_methods(data_rep[,run_rep, drop=F],  
                                est_mat = c('covariance', 'precision')[2],
                                method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                           'SpiecEasi', 'gCoDa','Spring')[7],
                                target_graph_cov = option$Sigma_list$A_cov,
                                target_graph_inv = option$Sigma_list$A_inv,
                                option = option)


rm('data_rep')

save.image(paste0('data/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_run_rep', run_rep, '.RData'))
