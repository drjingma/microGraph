#######################
##
##
##  Data Examples for book chapter
##
##
########################
filepath = 'E:\\Dropbox\\Microbial_Networks\\codes'
source('codes.R')
library(SpiecEasi)

data(amgut1.filt)
reference_data = amgut1.filt
dim(reference_data) #data format n by p
#---------------
# part I: null hypothesis
#---------------

# parameters that need to discuss:
n = c(100, 200, 500)
p = 200
reference_data = amgut1.filt
library_scale
alpha # for Dirichlet
mu = runif()
; Sigma # for pamametric



option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data)
option = list(hypothesis = 'null', model='Dirichelt', library_scale = library_scale, alpha = alpha)
option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma = Sigma)
option = 
function(n, p, 
         reference_data,
         option,
         ){
  if(option$hypothesis == 'null'){
    if(option$model =='shuffle'){
      #input real data set of n by p, shuffle real data
      model_data = null_data_generate_1(option$reference_data) 
    }
    if(option$model == 'Dirichlet'){
      model_data = null_data_generate_2(n, p, library_scale = option$library_scale, alpha = option$alpha) 
    }
  }
  
  if(option$hypothesis == 'alternative'){
    if(option$model == '')
    if(option$model == 'log-normal'){
      model_data = para_data_generate_2(n, p, mu = option$mu, Sigma = option$Sigma)
    }
  }

}
## Null 1: shuffle available data set x



alpha = runif(p)
null_data = null_data_generate_2(n, p, library_scale, alpha) # generate null from Dirichlet





#---------------
#part II: alternative hypothesis
#----------------

# use SpiecEasi generative model

graph = SpiecEasi_graph(amgut1.filt,'hub')

data_list= SpiecEasi_generate(amgut1.filt, graph)

graph = SpiecEasi_graph(amgut1.filt, 'erdos_renyi')

data_list= SpiecEasi_generate(amgut1.filt, graph)


# use log-normal, specify the precision matrix level, generate with covariance
Sigma = data_list$Correlation
mu = runif(n=p, min=0, max=10)
para_data_generate_2(n, p, mu, Sigma)
