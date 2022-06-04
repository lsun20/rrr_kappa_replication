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
# library("grplasso")	#glmnet group lasso requires the same n per group (ie multi-task learning), which is perfect for mlogit but wrong for the outcome model.
# library("mlogit")	#slower than unpenalized estimation in glmnet, but glmnet won't fit only an intercept
library("nnet")	#quicker multinomial logit
library("gglasso")
#library("plotrix")
#library("gridExtra")
#library("randomForest")
# setwd("~/Desktop/research/rrr_lasso_401k")

#######################
# clean and format data
#######################

# source('specifications.R')
# for (quintile in 1:5){
source('../specifications_kappa_AE98.R')
quintile=0 #quintile in (1-5). 0 means all quintiles

print(paste0('quintile: '))
print(paste0(quintile))

df  <- read.dta("m_d_806_processed.dta")
#df  <- read.dta("m_d_806_processed_subsample.dta")
spec=1 #spec in (1-3) # 3 is making mlogit slow
data<-get_data(df,spec,quintile) #trimming like Farrell; different than Chernozhukov et al. 
# SS: need to directly edit the script 'specifications_kappa_AE98.R' to change the covariate and instrument

Y=data[[1]]
T=data[[2]] #instrument
X=data[[3]] #no intercept

##################
# helper functions
##################

source('../primitives.R')
source('../stage1_kappa.R')

# dictionary
dict=b2 # b for partially linear model, b2 for interacted model
p=length(b(T[1],X[1,]))

#p0=dim(X0) used in low-dim dictionary in the stage 1 tuning procedure
p0=ceiling(p/4) #p/2 for quintile=0. 4 otw
if (p>60){
  p0=ceiling(p/40) #p/20 for quintile=0. 40 otw
}
p0=2#the minimum number of dict

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

source('../stage2_kappa.R')
results<-rrr(Y,T,X,p0,D_LB,D_add,max_iter,dict,alpha_estimator,gamma_estimator) 

theta_hat <- results[[1]]
j_hat <- results[[2]]
omega_hat <- results[[3]]
#moments_full <- results[[4]] 
psi_hat <- results[[5]] 

# multi2nd instruments
theta1_hat <- results[[1]]
j1_hat <- results[[2]]
omega1_hat <- results[[3]]
#moments_full <- results[[4]] 
psi1_hat <- results[[5]] 
result1 <- results

# samesex instruments
theta2_hat <- results[[1]]
j2_hat <- results[[2]]
omega2_hat <- results[[3]]
#moments_full <- results[[4]] 
psi2_hat <- results[[5]] 
result2 <- results



# Save multiple objects
#save(result1, result2, file = "educm_mean_spec1_full.RData")
save(result1, result2, file = "age2nd_mean_spec1_full.RData")

# To load the data again
#load("educm_mean_spec1_full.RData")
load("age2nd_mean_spec1_full.RData")

###########
# joint confidence bands
###########
n <- nrow(X)
psi_hat <- cbind(psi1_hat,psi2_hat)
omega_hat <- t(psi_hat) %*% psi_hat  / n

V <- omega_hat/n; V[1,1] <- V[1,1]/j1_hat^2; V[2,2] <- V[2,2]/j2_hat^2
V[1,2] <- V[1,2]/j1_hat/j2_hat; V[2,1] <- V[2,1]/j1_hat/j2_hat;


crit <- 1.96
theta1_l <- theta1_hat - crit * V[1,1]^0.5
theta1_u <- theta1_hat + crit * V[1,1]^0.5
theta2_l <- theta2_hat - crit * V[2,2]^0.5
theta2_u <- theta2_hat + crit * V[2,2]^0.5
tstat <- (theta1_hat - theta2_hat)/(V[1,1]-V[1,2]-V[2,1]+V[2,2])^0.5

cbind(theta1_l,theta1_u, V[1,1]^0.5)
cbind(theta2_l,theta2_u, V[2,2]^0.5)
pt(tstat, df = n-1, lower.tail = T) # one-sided test p-value
2 * pt( abs(tstat), df = n-1, lower.tail=FALSE) # two-sided test p-value