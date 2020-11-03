### for submission, summarize files to obtain plots
## useful data objects:
## null1.1_none, null2_none, alt1_data_100, alt1_data_1000, alt2_data_100, alt2_data_1000

## useful graph objects
## pp ( null1.1), pp1 (null2), alt1_plot_all_100, alt1_plot_all_1000, alt2_plot_all_100, alt2_plot_all_1000









library(ggplot2)
library(kableExtra)
library(cowplot)
theme_set(theme_cowplot())

filepath = '\\\\fs2-vip\\students\\yuek\\Desktop\\micro_net' #BOX
setwd(filepath)
source('summary_code_function.R')
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7","#F0E442")
part_name = list('1'='ReBoot', '2' = 'SparCC', '3' = 'CCLasso', '4' = 'COAT','5' = 'SPIEC-EASI', '6' = 'gCoda', '7' = 'SPRING')

#########
## results for null models
#########

###### 
## null1.1
folder_name= 'data/submission'
include=100

copula_distr = 'zinegbin'
choose_model = 'null1.1'
network_option='none'
network_condition_number = 0

fix = 'fix_lambda_'
fix = 'vary_'

null1.1= lapply(c(100, 200, 300, 500), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
  get_summary_fp(fix, folder_name, choose_model,n=n, p=127, nreps=200, part=part, part_name, copula_distr=copula_distr,
                 network_option = network_option, network_condition_number = network_condition_number))),
  dist = copula_distr, network_option = network_option, network_condition_number = network_condition_number))
null1.1 = do.call(rbind, null1.1)
head(null1.1)


null1.1_none = null1.1[null1.1$adjust=='none', ]
pp = ggplot(data = null1.1_none, aes(x=n, y=x, group=method, color=method, 
                                     linetype = method))+
  coord_cartesian(ylim = c(0, min(max(null1.1_none$x), 0.4)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  ylab('Type I Error')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[1:7])+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  scale_linetype_manual(values=1:7)+
  # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  scale_x_continuous(breaks=c(100,  200, 300,  400,  500))
ggsave(pp, filename = paste0('plot/(', include, ' reps)', fix, 'submission_', copula_distr, '_type_I_null1.1.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')

null1.1_none[null1.1_none$method==part_name[3],]


# summarize time information for each method under varying sample size

pp_time = ggplot(data = null1.1_none, aes(x=n, y=log10(time), group=method, color=method, 
                                          linetype = method))+
  # coord_cartesian(ylim = c(0, min(max(null1.1_none$time), 0.2)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  ylab('log10(time) (seconds)')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[1:7])+
  theme(legend.title = element_blank())+
  scale_linetype_manual(values=1:7)+
  # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  scale_x_continuous(breaks=c(100, 150, 200, 250, 300, 350, 400, 450, 500))
ggsave(pp_time, filename = paste0('plot/', fix, 'submission_time_null1.1.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')

##################
##### 
## null2

folder_name= 'data/submission'
include=100
copula_distr = 'none'
choose_model = 'null2'
network_option='none'
network_condition_number = 0
add_prex = 'fix_lambda_'

null2 = lapply(c(200), function(p){
  null2= lapply(c(100, 200, 300, 500), function(n)
    cbind(do.call(rbind, lapply(1:7,  function(part)
      get_summary_fp(add_prex, folder_name, choose_model,n=n, p=p, nreps=200, part=part, part_name, 
                     copula_distr=copula_distr,network_option = network_option, 
                     network_condition_number = network_condition_number))),
      dist = copula_distr, 
      network_option = network_option, 
      network_condition_number = network_condition_number))
  
  null2 = do.call(rbind, null2)
  null2
})

null2 = do.call(rbind, null2)
head(null2)

null2_none = null2[null2$adjust=='none' & null2$p==200,]
pp1 =ggplot(data = null2_none, aes(x=n, y=x, group=method, color=method, 
                                   linetype=method))+
  coord_cartesian(ylim = c(0, min(max(null2_none$x), 0.2)))+
  geom_line(size=1)+
  scale_linetype_manual(values=1:7)+
  geom_point(size=2, shape=18)+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ylab('Type I Error')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[c(1:7)])+
  scale_x_continuous(breaks=c(100,  200, 300,400,500))+
  scale_y_continuous(breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06))

ggsave(pp1, filename = paste0('plot/(', include, ' reps)fix_lambda_submission_type_I_null2.png'), width = 11.7*1.7, height = 11.7*1.7, units = 'cm')

pp_time = ggplot(data = null2_none, aes(x=n, y=log10(time), group=method, color=method, 
                                        linetype = method))+
  # coord_cartesian(ylim = c(0, min(max(null1.1_none$time), 0.2)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  ylab('log10(time) (seconds)')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[1:7])+
  theme(legend.title = element_blank())+
  scale_linetype_manual(values=1:7)+
  # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  scale_x_continuous(breaks=c(100,  200, 300,400,500))
ggsave(pp_time, filename = paste0('plot/fix_lambda_submission_time_null2.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')




library(ggpubr)
null_plot = ggarrange(pp,pp1,labels = c('A', 'B'), common.legend = T, legend = 'right' )
ggsave(null_plot, filename = paste0('plot/submission_type1error_null1.1_null2.png'), width = 11.7*1.7, height = 7*1.7, units = 'cm')


#####################
## results for alternative models
####################

###################
## alt1
#alt1
folder_name = 'data/submission'
# fix = 'fix_lambda_'
fix = 'vary_'
include = 20

{
  target = 'inv' # or 'cov'
  copula_distr = 'zinegbin'
  choose_model = 'alt1'
  network_option='erdos_renyi'
  network_condition_number = 100
}


{
  target = 'inv' # or 'cov'
  copula_distr = 'zinegbin'
  choose_model = 'alt1'
  network_option='erdos_renyi'
  network_condition_number = 1000 
}


if(T){
  data = lapply(c(100,  200, 300, 500), function(n){
    output_plot_all(fix,folder_name, choose_model,n=n, p=127, nreps=200, part_name, copula_distr,network_option, network_condition_number, target)$data
  })
  data = do.call(rbind, data)
  data$n = as.factor(data$n)
  
  tmpdata = data[data$group=='mean' & data$method != 'ReBoot',] # for
  plot_all = ggplot(tmpdata, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    facet_wrap(~n, labeller = label_both, nrow = 2)+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))
  ggsave(plot_all, filename = paste0('plot/(', include, ' reps)', fix, 'submission_ROC_alt1_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
  tmpdata2 = with(tmpdata, aggregate(time, by = list(n=n, p=p, method=method), FUN=mean, ))
  for(i in c(1))tmpdata2[,i] = as.numeric(levels(tmpdata2[,i]))[tmpdata2[,i]]
  plot_time = ggplot(tmpdata2, aes(x=n, y=log10(x), group=method, color = method,linetype=method))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    scale_shape_manual(values=2:7)+
    xlab('n')+
    ylab('log10(time) (seconds)')+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))+
    scale_x_continuous(breaks=c(100, 150, 200, 250))
  ggsave(plot_time, filename = paste0('plot/', fix, 'submission_time_alt1_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
}

# for condition number 100:
alt1_data_100 = data
alt1_plot_all_100 = plot_all
alt1_plot_time_100 = plot_time


# for condition number 1000:
alt1_data_1000 = data
alt1_plot_all_1000 = plot_all
alt1_plot_time_1000 = plot_time

###########
##alt 2
folder_name = 'data/submission'
fix = 'fix_lambda_'
include=20

{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='erdos_renyi'
  network_condition_number = 100
}


{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='erdos_renyi'
  network_condition_number = 1000 
  
} 

if(T){
  data = lapply(c(100, 200,300, 500), function(n){
    output_plot_all(fix, folder_name,choose_model,n=n, p=200, nreps=200, part_name, copula_distr,network_option, network_condition_number, target)$data
  })
  data = do.call(rbind, data)
  data$n = as.factor(data$n)
  
  tmpdata = data[data$group=='mean' & data$method != 'ReBoot',]
  plot_all = ggplot(tmpdata, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    facet_wrap(~n, labeller = label_both, nrow = 2)+
    scale_color_manual(values = cbPalette[c(2:5,6, 7)])+
    scale_linetype_manual(values = c(2:5,6, 7))
  ggsave(plot_all, filename = paste0('plot/(', include, ' reps)', fix, 'submission_ROC_alt2_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
  tmpdata2 = with(tmpdata, aggregate(time, by = list(n=n, p=p, method=method), FUN=mean, ))
  for(i in c(1))tmpdata2[,i] = as.numeric(levels(tmpdata2[,i]))[tmpdata2[,i]]
  plot_time = ggplot(tmpdata2, aes(x=n, y=log10(x), group=method, color = method,linetype=method))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    scale_shape_manual(values=c(2:5,6,7))+
    xlab('n')+
    ylab('log10(time) (seconds)')+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    scale_color_manual(values = cbPalette[c(2:5,6, 7)])+
    scale_linetype_manual(values = c(2:5,6,7))+
    scale_x_continuous(breaks=c(100, 200,300, 400, 500))
  ggsave(plot_time, filename = paste0('plot/', fix, 'submission_time_alt2_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
}

# for condition number 100:
alt2_data_100 = data
alt2_plot_all_100 = plot_all
alt2_plot_time_100 = plot_time


# for condition number 1000:
alt2_data_1000 = data
alt2_plot_all_1000 = plot_all
alt2_plot_time_1000 = plot_time


save.image('for_plts.RData')
