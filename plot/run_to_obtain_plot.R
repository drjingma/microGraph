###--------------------------
### run functions in this file to obtain plots in the manuscript
###---------------------------

## modify the file path and the data folder name/data file prefix name before using this script.


library(ggplot2)
library(kableExtra)
library(cowplot)
theme_set(theme_cowplot())

filepath = '......'
setwd(filepath)
source('plot_function.R')


#########
## results for null models
#########

###### 
## null1.1
get_plt_null1.1 = 
  plt_null1.1(n_list= c(100, 200, 300, 500),
              p=127,
              folder_name= 'data/submission',
              nrep=200,
              include=100,
              copula_distr = 'zinegbin',
              choose_model = 'null1.1',
              network_option='none',
              network_condition_number = 0,
              fix = 'vary_' # data file prefix name, if any
              )


ggsave(get_plt_null1.1$res, filename = 'type_I_null1.1.png',width =11.7*1.7, height= 11.7*1.7, units = 'cm')

ggsave(get_plt_null1.1$time, filename = 'time_null1.1.png',width =11.7*1.7, height= 11.7*1.7, units = 'cm')


##### 
## null2
get_plt_null2 = 
  plt_null2(
    n_list = c(100, 200, 300, 500),
    p = 200,
    folder_name= 'data/submission',
    nreps=200,
    include=100,
    copula_distr = 'none',
    choose_model = 'null2',
    network_option='none',
    network_condition_number = 0,
    fix = 'fix_lambda_')

ggsave(get_plt_null2$res, filename = 'type_I_null2.png', width = 11.7*1.7, height = 11.7*1.7, units = 'cm')
ggsave(get_plt_null2$time, filename = 'time_null2.png',width =11.7*1.7, height= 11.7*1.7, units = 'cm')


library(ggpubr)
null_plot = ggarrange(get_plt_null1.1$res,get_plt_null2$res,labels = c('A', 'B'), common.legend = T, legend = 'right' )
ggsave(null_plot, filename = paste0('type_I_null1.1_null2.png'), width = 11.7*1.7, height = 7*1.7, units = 'cm')


#####################
## results for alternative models
####################

##########
## alt1
get_plt_alt1_cond100 = 
  plt_alt(
    n_list = c(100, 200, 300, 500),
    p = 127,
    folder_name = 'data/submission',
    nreps = 200,
    include = 20,
    copula_distr = 'zinegbin',
    choose_model = 'alt1',
    network_option= 'erdos_renyi',
    network_condition_number = 100,
    fix = 'vary_',
    target='inv')

ggsave(get_plt_alt1_cond100$res, filename = 'ROC_alt1_cond100.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')
ggsave(get_plt_alt1_cond100$time, filename = 'time_alt1_cond100.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')

get_plt_alt1_cond1000 = 
  plt_alt(
    n_list = c(100, 200, 300, 500),
    p = 127,
    folder_name = 'data/submission',
    nreps = 200,
    include = 20,
    copula_distr = 'zinegbin',
    choose_model = 'alt1',
    network_option= 'erdos_renyi',
    network_condition_number = 1000,
    fix = 'vary_',
    target='inv')

ggsave(get_plt_alt1_cond1000$res, filename = 'ROC_alt1_cond1000.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')
ggsave(get_plt_alt1_cond1000$time, filename = 'time_alt1_cond1000.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')


###########
##alt 2

get_plt_alt2_cond100 = 
  plt_alt(
    n_list = c(100, 200, 300, 500),
    p = 200,
    folder_name = 'data/submission',
    nreps = 200,
    include = 20,
    copula_distr = 'none',
    choose_model = 'alt2',
    network_option= 'erdos_renyi',
    network_condition_number = 100,
    fix = 'fix_lambda_',
    target='inv')

ggsave(get_plt_alt2_cond100$res, filename = 'ROC_alt2_cond100.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')
ggsave(get_plt_alt2_cond100$time, filename = 'time_alt2_cond100.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')

get_plt_alt2_cond1000 = 
  plt_alt(
    n_list = c(100, 200, 300, 500),
    p = 200,
    folder_name = 'data/submission',
    nreps = 200,
    include = 20,
    copula_distr = 'none',
    choose_model = 'alt2',
    network_option= 'erdos_renyi',
    network_condition_number = 1000,
    fix = 'fix_lambda_',
    target='inv')

ggsave(get_plt_alt2_cond1000$res, filename = 'ROC_alt2_cond1000.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')
ggsave(get_plt_alt2_cond1000$time, filename = 'time_alt2_cond1000.png', width = 11.7*1.7, height = 8*1.7, units = 'cm')



save.image('all_plts.RData')
