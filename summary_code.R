filepath = '\\\\fs2-vip\\students\\yuek\\Desktop\\micro_net' #BOX
setwd(filepath)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7","#F0E442")

# choose_model = c('null1', 'null2', 'alt1', 'alt2')[2]
# n = c(100, 200, 500, 289)[1]
# p = c(200, 127)[1]
# nreps=100
# run_rep = 1:10 # just compute for 10 repetitions
# part = (1:7)[1] # for the 7 methods
part_name = list('1'='CoNet', '2' = 'SparCC', '3' = 'CCLasso', '4' = 'COAT','5' = 'SPIEC-EASI', '6' = 'gCoda', '7' = 'SPRING')

copula_distr = 'pois'             # poisson generating model, for nul1.1 and alt1
copula_distr = 'negbin'           # for null1.1 and alt1
copula_distr = 'zinegbin'         # for null1.1 and alt1
copula_distr = 'none'             # for varying library scale and mu~uniform(0,4), alt2; also for null1
#copula_distr = 'fixN'
#copula_distr = 'fixNandrmvnorm'
#copula_distr = 'fixNandrmvnorm_mu01'
copula_distr = 'rmvnorm_mu04'     # for varying library scale and mu~uniform(0,4), null2

# compute fp summary for null models, dist_data files
get_summary_fp = function(choose_model, n, p, nreps, part, part_name, copula_distr){
  # conet and Sparcc have p values, others only fp_null
  res_for_20 = do.call(rbind, lapply(1:20, function(run_rep){
    load(paste0(filepath,'/dist_data/',copula_distr, '/', choose_model,'/res_n_', n, '_p_', p, '_', choose_model, 
                '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData'))
    if(part %in% c(1, 2)){
      rbind(data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                       fp=roc$fp_null$fp_null, adjust = 'none', rep = run_rep),
            data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                       fp=roc$fp_null$fp_null_FDR, adjust = 'FDR', rep=run_rep))
    }else if(part %in% c(3, 4)){
      data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                 fp=roc$fp_null, adjust = 'none', rep=run_rep)
    }else{
      data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                 fp=roc$fp_inv_null, adjust = 'none', rep=run_rep)
    }
  }))
  res_for_20$method = factor(res_for_20$method, levels=do.call(c,part_name))
  
  with(res_for_20, aggregate(fp, by = list(n=n, p=p, model= model, method=method, adjust=adjust), function(x)mean(x,na.rm=T) ))
  
}

library(ggplot2)
library(kableExtra)
library(cowplot)
theme_set(theme_cowplot())


copula_distr = 'none'

choose_model = 'null1'

null1= lapply(c(100, 289), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
  get_summary_fp(choose_model,n=n, p=127, nreps=100, part=part, part_name, copula_distr=copula_distr))),
  dist = copula_distr))
null1 = do.call(rbind, null1)
head(null1)


pp = ggplot(data = null1, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
  coord_cartesian(ylim = c(0, min(max(null1$x), 0.2)))+
  geom_line(data=null1, aes(linetype = interaction(method, adjust)))+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))
ggsave(pp, filename = paste0('plot/', copula_distr, '_', choose_model,'.png'))



# null1.1
copula_distr = 'zinegbin'

choose_model = 'null1.1'

null1.1= lapply(c(100, 120, 150, 200, 289), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
  get_summary_fp(choose_model,n=n, p=127, nreps=100, part=part, part_name, copula_distr=copula_distr))),
  dist = copula_distr))
null1.1 = do.call(rbind, null1.1)
head(null1.1)



pp = ggplot(data = null1.1, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
  coord_cartesian(ylim = c(0, min(max(null1.1$x), 0.2)))+
  geom_line(data=null1.1, aes(linetype = interaction(method, adjust)))+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))
ggsave(pp, filename = paste0('plot/', copula_distr, '_', choose_model,'.png'))

null1.1_none = null1.1[null1.1$adjust=='none', ]
pp = ggplot(data = null1.1_none, aes(x=n, y=x, group=method, color=method, 
                                linetype = method))+
  coord_cartesian(ylim = c(0, min(max(null1.1_none$x), 0.2)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  ylab('Type I Error')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[1:7])+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  scale_x_continuous(breaks=c(100, 120, 150, 200, 289))
ggsave(pp, filename = paste0('plot/type_I_null1.1.png'),width = 7, height=7)


# # null1 models
# null1.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null1',n=289, p=127, nreps=100, part=part, part_name)))
# null1.2$x = round(null1.2$x, 4)
# kable(null1.2, 'latex', booktab=T, col.names = c(names(null1.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# 
# null1.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null1',n=100, p=127, nreps=100, part=part, part_name)))
# null1.2$x = round(null1.2$x, 4)
# kable(null1.2, 'latex', booktab=T, col.names = c(names(null1.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()

# null2
copula_distr = 'rmvnorm_mu04'

choose_model = 'null2'
null2 = lapply(c(20, 50, 100, 150, 200), function(p){
  null2= lapply(c(100, 200, 500), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
  get_summary_fp(choose_model,n=n, p=p, nreps=100, part=part, part_name, copula_distr=copula_distr))),
  dist = copula_distr))
null2 = do.call(rbind, null2)
null2
})

null2 = do.call(rbind, null2)
head(null2)

pp1 = ggplot(data = null2, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
  coord_cartesian(ylim = c(0, min(max(null2$x), 0.2)))+
  geom_line(data=null2, aes(linetype = interaction(method, adjust)))+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  facet_wrap(~p)
pp2 = ggplot(data = null2, aes(x=p, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
  coord_cartesian(ylim = c(0, min(max(null2$x), 0.2)))+
  geom_line(data=null2, aes(linetype = interaction(method, adjust)))+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  facet_wrap(~n)

ggsave(pp1, filename = paste0('plot/', copula_distr, '_', choose_model,'n.png'))
ggsave(pp2, filename = paste0('plot/', copula_distr, '_', choose_model,'p.png'))

null2_none = null2[null2$adjust=='none' & null2$p==200,]
pp1 =ggplot(data = null2_none, aes(x=n, y=x, group=method, color=method, 
                              linetype=method))+
  coord_cartesian(ylim = c(0, min(max(null2_none$x), 0.2)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  ylab('Type I Error')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[c(1:7)])+
  scale_x_continuous(breaks=c(100,  200, 300,400,500))+
  scale_y_continuous(breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06))

ggsave(pp1, filename = paste0('plot/type_I_null2.png'), width = 7, height = 7)

# null2.1 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=100, p=200, nreps=100, part=part, part_name)))
# null2.1$x = round(null2.1$x, 4)
# kable(null2.1, 'latex', booktab=T, col.names = c(names(null2.1)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# null2.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=200, p=200, nreps=100, part=part, part_name)))
# null2.2$x = round(null2.2$x, 4)
# kable(null2.2, 'latex', booktab=T, col.names = c(names(null2.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# 
# null2.3= do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=500, p=200, nreps=100, part=part, part_name)))
# null2.3$x = round(null2.3$x, 4)
# kable(null2.3, 'latex', booktab=T, col.names = c(names(null2.3)[1:5], 'False Positive Rate'))%>%
#   kable_styling()

#######################################################################
# get ROC summary for alternative models, dist_data files
get_summary_roc = function(choose_model, n, p, nreps, part, part_name, copula_distr){
  # all models we compare the inverse cov graph ROC
  # to have an 'averaged' ROC curve, we may take average of fp and tp at each grid? (since we are using the same sequence of lambda)
  res_for_20 = lapply(1:20, function(run_rep){
    load(paste0(filepath,'/dist_data/', copula_distr, '/', choose_model,'/res_n_', n, '_p_', p, '_', choose_model, 
                '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData'))
    list(n=n, p=p, model = choose_model, method = part_name[[part]],rep=run_rep,
         tp = roc$ROC_inv$tp, fp = roc$ROC_inv$fp, auc = roc$ROC_inv$AUC)

  })
  
  ave_tp = rowMeans(do.call(cbind, sapply(res_for_20, `[`, 'tp')))
  ave_fp = rowMeans(do.call(cbind, sapply(res_for_20, `[`, 'fp')))
  ave_auc = mean(do.call(cbind, sapply(res_for_20, `[`, 'auc')))
    
  ave_res = list(n=n, p=p, model = choose_model, method = part_name[[part]],
                 ave_tp = ave_tp, ave_fp = ave_fp, ave_auc = ave_auc)
  return(list(ave_res = ave_res, res_for_20 = res_for_20))
}


plot_roc = function(roc_res){
  tmp1=data.frame(tp = roc_res$ave_res$ave_tp, fp = roc_res$ave_res$ave_fp, group=rep('mean', length(roc_res$ave_res$ave_fp)))
  tmp2= do.call(rbind, lapply(1:20, function(i)
    data.frame(tp = roc_res$res_for_20[[i]]$tp, fp = roc_res$res_for_20[[i]]$fp, 
               group=rep(as.character(i), length(roc_res$res_for_20[[i]]$fp)))))
  tmp = rbind(tmp1, tmp2)  
  tmp$n = roc_res$ave_res$n[1]
  tmp$p = roc_res$ave_res$p[1]
  tmp$model = roc_res$ave_res$model[1]
  tmp$method = roc_res$ave_res$method[1]

  
  gplot = ggplot(tmp, aes(x=fp, y=tp, group=group, color=group))+
    geom_line()+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    ggtitle(paste('n',tmp$n[1], 'p', tmp$p[1], tmp$model[1], tmp$method[1]))+
    xlab('False Positive')+
    ylab('True Positive')+
    scale_color_manual(values = c('black', rep('grey',10)))+
    scale_linetype_manual(values = c(1, rep(4, 10)))+
    geom_line(data = tmp[tmp$group=='mean',], aes(x=fp, y=tp))+
    theme(legend.position = 'none')
  
  return(list(plot = gplot, data = tmp))
}

# if to plot all methods together:
output_plot_all = function(choose_model, n, p, nreps, part_name, copula_distr){
  data = do.call(rbind, lapply(2:7, function(part){
    roc_res = get_summary_roc(choose_model, n, p, nreps, part, part_name, copula_distr)
    plot_roc(roc_res)$data
  }))
  data$group = factor(data$group, levels = c(1:20, 'mean'))
  data$method = factor(data$method, levels=do.call(c,part_name))
  
  plot_all = ggplot(data, aes(x=fp, y=tp, group=group, color = group))+
    geom_line()+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('dist', copula_distr, 'n',data$n[1], 'p', data$p[1], data$model[1]))+
    theme(legend.position = 'none')+
    scale_color_manual(values=c(rep('grey', 20), 1))+
    facet_wrap(~method)
  
  # ggsave(plot_all, filename = paste0('plot/', copula_distr, '_n',data$n[1], 'p', data$p[1], data$model[1],'.png'))
  return(list(data= data, plots = plot_all))
}



#alt1

copula_distr = 'zinegbin'

data = lapply(c(100, 120,150,200,289 ), function(n){
    output_plot_all('alt1',n=n, p=127, nreps=100, part_name, copula_distr)$data
})
data = do.call(rbind, data)
data$n = as.factor(data$n)


tmpdata = data[data$group=='mean',]
plot_all = ggplot(tmpdata, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
  geom_line(size=1.2)+
  geom_abline(slope=1, intercept=0,linetype = 3)+
  xlab('False Positive')+
  ylab('True Positive')+
  #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
  theme(legend.position = 'right')+
  facet_wrap(~n, labeller = label_both)+
  scale_color_manual(values = cbPalette[2:7])
ggsave(plot_all, filename = paste0('plot/ROC_alt1.png'), width = 14, height = 9)

# output_plot_all('alt1',n=100, p=127, nreps=100, part_name, copula_distr)
# output_plot_all('alt1',n=120, p=127, nreps=100, part_name, copula_distr)
# output_plot_all('alt1',n=150, p=127, nreps=100, part_name, copula_distr)
# output_plot_all('alt1',n=200, p=127, nreps=100, part_name, copula_distr)
# output_plot_all('alt1',n=289, p=127, nreps=100, part_name, copula_distr)


# alt2
copula_distr = 'edge_2_3'
# copula_distr = '2pnetwork'

data2 = lapply(c(100, 200, 500), function(n){
  tmp = lapply(c(20, 50, 100, 150, 200), function(p){
    output_plot_all('alt2',n=n, p=p, nreps=100, part_name, copula_distr)$data
  })
  do.call(rbind, tmp)
})
data2 = do.call(rbind, data2)
data2$p = as.factor(data2$p)

for(nn in c(100, 200, 500)){
  tmpdata2 = data2[data2$group=='mean'& data2$n==nn,]
  plot_all = ggplot(tmpdata2, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('sample size n =', nn))+
    theme(legend.position = 'right')+
    scale_color_manual(values=cbPalette[2:7])+
    facet_wrap(~p, labeller = label_both)
  ggsave(plot_all, filename = paste0('plot/ROC_alt2_n_',nn,'.png'), width = 14, height = 9)
  
}

for(pp in c(20,50,100,150,200)){
  tmpdata2 = data2[data2$group=='mean'& data2$p==pp,]
  plot_all = ggplot(tmpdata2, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('number of genes p =', pp))+
    theme(legend.position = 'right')+
    scale_color_manual(values=cbPalette[2:7])+
    facet_wrap(~n, labeller = label_both)
  ggsave(plot_all, filename = paste0('plot/ROC_alt2_p_',pp,'.png'), width = 14, height = 4.5)
  
}




data3 = lapply(c(100, 200, 500), function(n){
  tmp = lapply(c(20, 50, 100, 150, 200), function(p){
    output_plot_all('alt3',n=n, p=p, nreps=100, part_name, copula_distr)$data
  })
  do.call(rbind, tmp)
})
data3 = do.call(rbind, data3)
data3$p = as.factor(data3$p)

for(nn in c(100, 200, 500)){
  tmpdata3 = data3[data3$group=='mean'& data3$n==nn,]
  plot_all = ggplot(tmpdata3, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('sample size n =', nn))+
    theme(legend.position = 'right')+
    scale_color_manual(values=cbPalette[2:7])+
    facet_wrap(~p, labeller = label_both)
  ggsave(plot_all, filename = paste0('plot/ROC_alt3_n_',nn,'.png'), width = 14, height = 9)
  
}

for(pp in c(20,50,100,150,200)){
  tmpdata3 = data3[data3$group=='mean'& data3$p==pp,]
  plot_all = ggplot(tmpdata3, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('number of genes p =', pp))+
    theme(legend.position = 'right')+
    scale_color_manual(values=cbPalette[2:7])+
    facet_wrap(~n, labeller = label_both)
  ggsave(plot_all, filename = paste0('plot/ROC_alt3_p_',pp,'.png'), width = 14, height = 4.5)
  
}


# output_plot_all('alt2',n=200, p=200, nreps=100, part_name, copula_distr)
# output_plot_all('alt2',n=500, p=200, nreps=100, part_name, copula_distr)

# inspect CClasso problem
nn=100
for(pp in c(20, 50, 100, 150, 200)){
  plot_all =output_plot_all('alt2',n=nn, p=pp, nreps=100, part_name, copula_distr)$plots
  ggsave(plot_all, filename = paste0('plot/', copula_distr, '_n_',nn,'_p_',pp,'alt2.png'))
}


# output_plot_all('alt2',n=200, p=100, nreps=100, part_name, copula_distr)$plots
# output_plot_all('alt2',n=200, p=150, nreps=100, part_name, copula_distr)$plots
# output_plot_all('alt2',n=200, p=200, nreps=100, part_name, copula_distr)$plots