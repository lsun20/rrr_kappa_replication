########
# set up
########

rm(list=ls())
install.packages(c('glmnet','grplass','gglasso','plotrix'))
install.packages(c('glmnet'))

library("foreign")
library("dplyr")
library("ggplot2")
# library("quantreg")

library("MASS")
library("glmnet")	#lasso, group lasso, and ridge, for outcome models, LPM, mlogit. Also, unpenalized mlogit via glmnet is far faster than the mlogit package.
library("grplasso")	#glmnet group lasso requires the same n per group (ie multi-task learning), which is perfect for mlogit but wrong for the outcome model.
# library("mlogit")	#slower than unpenalized estimation in glmnet, but glmnet won't fit only an intercept
library("nnet")	#quicker multinomial logit
library("gglasso")
library("plotrix")
library("gridExtra")

# setwd("~/Desktop/research/rrr_lasso_401k")
# loop the simulations

#######################
# clean and format data
#######################
N_sim <- 500; N <- 1000; # simulation spec
out <- vector() # initialize an empty vector to store
source('specifications_sim.R')
for (s in 1:N_sim){
  data<-get_sim_data(N,s) 
  #SS: need to directly edit the script 'specifications_sim.R' to supplement grid; first column of Y is participation
  Y=data[[1]]
  T=data[[2]] #eligibility
  X=data[[3]] #no intercept
  perc=data[[4]] #grid pts
  
  ##################
  # helper functions
  ##################
  
  source('primitives.R')
  source('stage1_kappa.R')
  
  # dictionary
  dict=b2 # b for partially linear model, b2 for interacted model
  p=length(b(T[1],X[1,]))
  
  #p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
  p0=ceiling(p/4) #p/2 for quintile=0. 4 otw
  if (p>60){
    p0=ceiling(p/40) #p/20 for quintile=0. 40 otw
  }
  
  
  D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
  D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
  max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights
  
  ###########
  # algorithm
  ###########
  
  set.seed(1) # for sample splitting
  lasso=1 #else dantzig
  
  source('stage2_kappa.R')
  results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,lasso) 
  
  theta_hat <- results[[1]]
  out <- rbind( out, t(theta_hat))

}

# dir <- './results/dml_rrr_sim_cntr_200105/'
dir <- './results/dml_rrr_sim_cntr_200105/'
robustness <- '_lambda_L5_keep'
spec <- 'untreated'
filename <- paste(dir,spec,robustness,'.csv',sep = "")
write.csv(out, filename,row.names = FALSE)


# parallelize the simulations
###########
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# load auxilliary functions
source('primitives.R')
source('stage1_kappa.R')
N_sim <- 500; N <- 1000; # simulation spec
source('specifications_sim.R')
data<-get_sim_data(N,1) # get one draw of data to set up parameter values
  #SS: need to directly edit the script 'specifications_sim.R' to supplement grid; first column of Y is participation
T=data[[2]] #eligibility
X=data[[3]] #no intercept
perc=data[[4]] #grid pts
# dictionary
dict=b2 # b for partially linear model, b2 for interacted model
p=length(b(T[1],X[1,]))

#p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
p0=ceiling(p/4) #p/2 for quintile=0. 4 otw
if (p>60){
  p0=ceiling(p/40) #p/20 for quintile=0. 40 otw
}


D_LB=0 #each diagonal entry of \hat{D} lower bounded by D_LB
D_add=.2 #each diagonal entry of \hat{D} increased by D_add. 0.1 for 0, 0,.2 otw
max_iter=10 #max number iterations in Dantzig selector iteration over estimation and weights

lasso=1 #else dantzig

source('stage2_kappa.R')

finalMatrix <- foreach(i=1:N_sim, .combine=rbind) %dopar% {
  data<-get_sim_data(N,i) 
  Y=data[[1]]
  T=data[[2]] #eligibility
  X=data[[3]] #no intercept
  
  set.seed(1) # for sample splitting
  results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,lasso) 
  
  theta_hat <- results[[1]]
  t(theta_hat)
  # tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
}
#stop cluster
stopCluster(cl)
dir <- ''
robustness <- '_lambda_L5_keep'
spec <- 'untreated'
filename <- paste(dir,spec,robustness,'.csv',sep = "")
write.csv(finalMatrix, filename,row.names = FALSE)

