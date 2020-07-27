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


lambda_seq = c(seq(0.01,0.1,0.01),seq(0.1,3,0.1)) * sqrt(log(p)/n) # Jing's setting for specifying the lambda list
lambda_seq = sort(lambda_seq, decreasing = T)


# compare methods for their ROC

cov_conet_roc = compare_methods(data_rep[,1, drop=F],
                                est_mat = c('covariance', 'precision')[1],
                                method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                           'SpiecEasi', 'gCoDa','Spring')[1],
                                target_graph_cov = option$Sigma_list$A_cov,
                                target_graph_inv = option$Sigma_list$A_inv,
                                option = option)

cov_sparcc_roc = compare_methods(data_rep[,run_rep, drop=F],  
                                 est_mat = c('covariance', 'precision')[1],
                                 method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                            'SpiecEasi', 'gCoDa','Spring')[2],
                                 target_graph_cov = option$Sigma_list$A_cov,
                                 target_graph_inv = option$Sigma_list$A_inv,
                                 option = option)

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
