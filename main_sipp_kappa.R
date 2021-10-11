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

#######################
# clean and format data
#######################

# source('specifications.R')
# for (quintile in 1:5){
source('../specifications_kappa.R')
quintile=0 #quintile in (1-5). 0 means all quintiles

  print(paste0('quintile: '))
  print(paste0(quintile))
  
df  <- read.dta("../sipp1991.dta")
spec=3 #spec in (1-3)
data<-get_data(df,spec,quintile) #trimming like Farrell; different than Chernozhukov et al. 
# SS: need to directly edit the script 'specifications_kappa.R' to supplement grid; first column of Y is participation

Y=data[[1]]
T=data[[2]] #eligibility
X=data[[3]] #no intercept
perc=data[[4]] #grid pts

##################
# helper functions
##################

source('../primitives.R')
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

alpha_estimator=1
gamma_estimator=1

#alpha_estimator: 0 dantzig, 1 lasso
#gamma_estimator: 0 dantzig, 1 lasso, 2 rf, 3 nn

###########
# algorithm
###########

set.seed(1) # for sample splitting

source('stage2_kappa.R')
results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator) 

theta_hat <- results[[1]]
j_hat <- results[[2]]
omega_hat <- results[[3]]
###########
# joint confidence bands
###########
n <- nrow(X)
V <- j_hat^-1 *omega_hat * j_hat^-1
theta_hat_var = diag(V)/n
sigma <- diag(diag(V)^-0.5 ) %*% V %*% diag(diag(V)^-0.5 )
R <- abs(mvrnorm(10000, rep(0,length(theta_hat)), sigma))
R_max <- apply(R, 1,FUN=max) # take row-wise max
crit <- quantile(R_max,0.95)
theta_l <- theta_hat - crit * theta_hat_var^0.5
theta_u <- theta_hat + crit * theta_hat_var^0.5

# printer(results) #SS: haven't updated these to show vectors
# for_tex(results)
out <- cbind( perc, theta_hat, theta_l, theta_u)
dir <- '../results/sipp_cntr_210923/'
spec <- 'untreated'
filename <- paste(dir,spec,'_CB','.csv',sep = "")
write.csv(out, filename,row.names = FALSE)


# }