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
#library("gglasso")
#library("plotrix")
#library("gridExtra")
#library("randomForest")
# setwd("~/Desktop/research/rrr_lasso_401k")
library(boot)
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
spec=0 #spec in (1-3) # 3 is making mlogit slow, 0 replicates AE98
data<-get_data(df,spec,quintile) #trimming like Farrell; different than Chernozhukov et al. 
# SS: need to directly edit the script 'specifications_kappa.R' to supplement grid; first column of Y is participation

Y=data[[1]]
T=data[[2]] #eligibility
X=data[[3]] #no intercept
p <- dim(X)[2]
# perc=data[[4]] #grid pts


##################
# Replicate AF13 (via analytical two-step estimator asymptotics)
mylogit <- glm(T ~ X, data = data, family = "binomial")
n <- length(T)
tau <- mylogit$fitted.values
D <- Y[,1]
kappa <- 1 - D*(1-T)/(1-tau) - (1-D)*T/tau
theta <- t(X[,1])%*%kappa/sum(kappa)

Mtheta <- 2*mean(kappa)
# Calculate standard error
dg <- -2*(X[,1]-theta)

#sigmoid function, inverse of logit
sigmoid <- function(z){1/(1+exp(-z))}
Fx <- sigmoid(cbind(1,X)%*%mylogit$coefficients)
dtau <- cbind(1,X)*((((1-Fx)*Fx)%*%rep(1,5)))
dkappa <- (as.matrix(- D*(1-T)/(1-tau)^2 + (1-D)*T/tau^2)%*%rep(1,5))*dtau
Mgamma <- t(dkappa)%*%dg/n
#gradient function
grad <- function(theta, X, y){
  
  h <- sigmoid(X%*%theta)
  grad <- (X*((y - h)%*%rep(1,dim(X)[2]))) 
  grad
}

hess <- function(theta, X, y){
  m <- length(y)
  h <- sigmoid(X%*%theta)
  X.h <- X*((((1-h)*h)%*%rep(1,dim(X)[2])))
  hess <- - t(X)%*%X.h/m
  hess
}
g <- grad(mylogit$coefficients,cbind(1,X),T)
# confirm they match with standard error from coef(summary(mylogit))
H <- hess(mylogit$coefficients,cbind(1,X),T)
fisher_info  <- -solve(H)
diag(fisher_info%*%(t(g)%*%g/n)%*%fisher_info)^0.5/sqrt(n) # doesn't replicate
diag(fisher_info)^0.5/sqrt(n) # glm uses information equality


#diag(ginv(t(g)%*%g/n))^0.5
psi <- g%*%solve(H) # influence function for the logit coefficient estimates

meat <- kappa*dg + psi%*%Mgamma
V <- mean(meat^2)/(Mtheta^2)
sqrt(V/n)


##################
# Replicate AF13 via bootstrap
# bootstrap using parallelization
# make a data frame
df_subset <- as.data.frame(cbind(Y,df$multi2nd,df$samesex,X))
colnames(df_subset) <- c("takeup","outcome","T1","T2","X1","X2","X3","X4","X5","X6")
boot.kappa<- function(data,indices){
  data.b<- data[indices,];  
  D=data.b$takeup #eligibility
  X=data.b[,-1:-4] #no intercept
  mylogit1 <- glm(T1 ~ X1+X2+X3+X4+X5+X6, data = data.b, family = "binomial")
  tau <- mylogit1$fitted.values
  T=data.b$T1 #eligibility
  kappa <- 1 - D*(1-T)/(1-tau) - (1-D)*T/tau
  theta <- t(X[,1])%*%kappa/sum(kappa)
  
  mylogit2 <- glm(T2 ~ X1+X2+X3+X4+X5+X6, data = data.b, family = "binomial")
  tau <- mylogit2$fitted.values
  T=data.b$T2 #eligibility
  kappa <- 1 - D*(1-T)/(1-tau) - (1-D)*T/tau
  theta <- cbind(theta,t(X[,1])%*%kappa/sum(kappa))
  return(theta/4) # if ageq2nd then divide by 4
}

result.boot.kappa<- boot(data=df_subset, statistic=boot.kappa, parallel="multicore", ncpus = 16, R=200)

educm_complier <- result.boot.kappa
ageq2nd_complier <- result.boot.kappa

save(educm_complier, ageq2nd_complier, file = "complier_mean_spec0_X6_full_kappa_weight_bs.RData")
 