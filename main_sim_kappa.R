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
library("randomForest")

# setwd("~/Desktop/research/rrr_lasso_401k")

# parallelize the simulations
###########
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
cl <- makeCluster(12)
registerDoParallel(cl)

# load auxilliary functions
source('../primitives.R')
source('stage1_kappa.R')
N_sim <- 1000; N <- 1000; # simulation spec
source('../specifications_sim.R')
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

alpha_estimator=1
gamma_estimator=1

#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn

#source('stage2_kappa.R') # rrr() only calculates auto-dml
source('stage2_kappa_dml.R') # rrr() also calculates dml 
source('naive_kappa.R') # naive_kappa() kappa weights
rm(finalMatrix)
finalMatrix <- foreach(i=1:N_sim, .combine=rbind, .packages=c('glmnet','MASS')) %dopar% {
  data<-get_sim_data(N,i) 
  Y=data[[1]]
  T=data[[2]] #eligibility
  X=data[[3]] #no intercept
  
  set.seed(1) # for sample splitting
  results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator) 
  
  theta_hat <- results[[1]]
  # t(theta_hat)
  # tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)

  ###########
  # joint confidence bands
  ###########
  j_hat <- results[[2]]
  omega_hat <- results[[3]]
  V <- j_hat^-1 *omega_hat * j_hat^-1
  theta_hat_var = diag(V)/N
  diagV <- diag(V)^-0.5; diagV[diagV==Inf] <- 0;
  sigma <- diag(diagV ) %*% V %*% diag(diagV )
  R <- abs(mvrnorm(10000, rep(0,length(theta_hat)), sigma))
  R_max <- apply(R, 1,FUN=max) # take row-wise max
  crit <- quantile(R_max,0.95)
  theta_l <- theta_hat - crit * theta_hat_var^0.5
  theta_u <- theta_hat + crit * theta_hat_var^0.5
  #cbind(t(theta_hat),t(theta_l),t(theta_u))
  theta_kappa_hat <- naive_kappa(Y,T,X)[[1]]
  theta_dml_hat <- results[[4]]
  
  cbind(t(theta_hat),t(theta_l),t(theta_u),t(theta_hat_var),t(theta_dml_hat),t(theta_kappa_hat))
}

dir <- '../results/dml_rrr_sim_cntr_211004/'
robustness <- '_lambda_L5_keep_full'
spec <- 'untreated'
filename <- paste(dir,spec,robustness,'.csv',sep = "")
write.csv(finalMatrix, filename,row.names = FALSE)

#stop cluster
stopCluster(cl)
