#stage 1 function
l=0.1

RMD_dantzig <- function(M, G, D, lambda=0, sparse = TRUE) {
  
  p <- ncol(G)
  zp <- rep(0, p)
  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,T,X,TX)
  
  A <- solve(diag(D),G)
  R <- rbind(A, -A)
  
  a <- solve(diag(D),M)
  r <- c(a - lambda*L, -a - lambda*L)
  
  if(sparse) {
    Ip <- as(p, "matrix.diag.csr")
    R <- as.matrix.csr(R)
    f <- rq.fit.sfnc(Ip, zp, R = R, r = r)
  } else {
    Ip <- diag(p)
    f <- rq.fit.fnc(Ip, zp, R = R, r = r)
  }
  
  return(f)
}

RMD_lasso <- function(M, G, D, lambda=0, control = list(maxIter = 1000, optTol = 10^(-5), 
                                                        zeroThreshold = 10^(-6)), beta.start = NULL) {
  
  p <- ncol(G)
  Gt <- G
  Mt <- M

  L <-c(l,rep(1,p-1)) #dictionary is ordered (constant,T,X,TX)
  lambda_vec=lambda*L*D
  
  if (is.null(beta.start)) {
    beta <- rep(0,p) #vs low-dimensional initialization
  }
  else {
    beta <- beta.start
  }
  wp <- beta
  mm <- 1
  while (mm < control$maxIter) {
    beta_old <- beta
    for (j in 1:p) {
      rho=Mt[j]-Gt[j,]%*%beta+Gt[j,j]*beta[j]
      z=Gt[j,j]
      
      if (sum(is.na(rho)) >= 1) {
        beta[j] <- 0
        next
      }
      if (rho < -1 * lambda_vec[j]) 
        beta[j] <- (rho+lambda_vec[j])/z
      if (abs(rho) <= lambda_vec[j]) 
        beta[j] <- 0
      if (rho > lambda_vec[j]) 
        beta[j] <- (rho-lambda_vec[j])/z
    }
    wp <- cbind(wp, beta)
    if (sum(abs(beta - beta_old), na.rm = TRUE) < control$optTol) {
      break
    }
    mm <- mm + 1
  }
  w <- beta
  w[abs(w) < control$zeroThreshold] <- 0
  return(list(coefficients = w, coef.list = wp, num.it = mm))
}


get_D <- function(Y,T,X,m,rho_hat,b){
  n=nrow(X)
  p=length(b(T[1],X[1,]))
  
  df=matrix(0,p,n)
  for (i in 1:n){
    df[,i]=b(T[i],X[i,])*as.vector(rho_hat %*% b(T[i],X[i,]))-m(Y[i],T[i],X[i,],b)
  }
  df=df^2
  D2=rowMeans(df)
  
  D=sqrt(D2)
  return(D) #pass around D as vector
}

c=0.5
alpha=0.1
tol=1e-6


RMD_stable<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,is_alpha,is_lasso){
  
  k=1
  
  p=length(b(T[1],X[1,]))
  n=length(T)
  d = ncol(Y) #SS: Y is a matrix [D,etc]
  N_hat0_list <- vector("list", d)  
  # low-dimensional moments
  X0=X[,1:p0]
  for (dd in 1:d) {
    MNG0<-get_MNG(Y[,dd],T,X0,b)
    M_hat0=MNG0[[1]]
    N_hat0=MNG0[[2]]
    G_hat0=MNG0[[3]]
    N_hat0_list[[dd]] <- N_hat0
  }
  
  # initial estimate
  rho_hat0=solve(G_hat0,M_hat0)
  rho_hat=c(rho_hat0,rep(0,p-ncol(G_hat0)))
  beta_hat0 = matrix(0,ncol(G_hat0),d)
  beta_hat <- matrix(0,p,d)
  for (dd in 1:d) {
    beta_hat0[, dd] = solve(G_hat0, N_hat0_list[[dd]])
    beta_hat[, dd] = c(beta_hat0[, dd],rep(0,p-ncol(G_hat0)))
  }
  # moments
  N_hat_list <- vector("list", d)
  for (dd in 1:d) {
    MNG<-get_MNG(Y[,dd],T,X,b)
    M_hat=MNG[[1]]
    N_hat=MNG[[2]]
    G_hat=MNG[[3]]
    N_hat_list[[dd]] <- N_hat
  }
  
  # penalty
  lambda=c*qnorm(1-alpha/(2*p))/sqrt(n) # snippet
  
  if(is_alpha){ 
    ###########
    # alpha_hat
    ###########
    diff_rho=1
    while(diff_rho>tol & k<=max_iter){
      
      # previous values
      rho_hat_old=rho_hat+0
      
      # normalization
      D_hat_rho=get_D(Y[,1],T,X,m,rho_hat_old,b) #SS: pick Y[,1] for the code to go thru but Y does not matter for m
      D_hat_rho=pmax(D_LB,D_hat_rho)
      D_hat_rho=D_hat_rho+D_add
      
      # RMD estimate
      if(is_lasso){
        rho_hat=RMD_lasso(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }else{
        rho_hat=RMD_dantzig(M_hat, G_hat, D_hat_rho, lambda)$coefficients
      }
      
      # difference
      diff_rho=two.norm(rho_hat-rho_hat_old)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(rho_hat)
    
  } else { 
    ###########
    # gamma_hat
    ###########
    diff_beta_max=1
    while(diff_beta_max>tol & k<max_iter){
      
      # previous values
      beta_hat_old=beta_hat+0
      
      # normalization
      D_hat_beta_list <- vector("list", d)
      for (dd in 1:d) {
        D_hat_beta = get_D(Y[, dd],T,X,m2,beta_hat_old[, dd],b)
        D_hat_beta = pmax(D_LB,D_hat_beta)
        D_hat_beta = D_hat_beta + D_add
        D_hat_beta_list[[dd]] <-  D_hat_beta
      }      
      # RMD estimate
      if(is_lasso){
        for (dd in 1:d) {
          beta_hat[, dd] <- RMD_lasso(N_hat, G_hat, D_hat_beta_list[[dd]], lambda)$coefficients
        }
      }else{
        for (dd in 1:d) {
          beta_hat[, dd] <- RMD_dantzig(N_hat, G_hat, D_hat_beta_list[[dd]], lambda)$coefficients
        }
      }
      
      # difference
      diff_beta <- rep(0,d)
      for (dd in 1:d) {
        diff_beta[dd]=two.norm(beta_hat[,dd]-beta_hat_old[,dd])
      }
      diff_beta_max=max(diff_beta)
      k=k+1
      
    }
    
    print(paste0('k: '))
    print(paste0(k))
    return(beta_hat)
    
  }
}

arg_Forest<- list(clas_nodesize=1, reg_nodesize=5, ntree=1000, na.action=na.omit, replace=TRUE)
arg_Nnet<- list(size=8,  maxit=1000, decay=0.01, MaxNWts=10000,  trace=FALSE)

get_stage1<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,alpha_estimator,gamma_estimator){
  
  p=length(b(T[1],X[1,]))
  n=length(T)
  
  ###########
  # alpha hat
  ###########
  if(alpha_estimator==0){ # dantzig
    
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,0)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  } else if(alpha_estimator==1){ # lasso
    
    rho_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,1,1)
    alpha_hat<-function(d,z){
      return(b(d,z)%*%rho_hat)
    }
    
  }
  
  ###########
  # gamma hat
  ###########  
  MNG<-get_MNG(Y[,1],T,X,b) #SS: pick Y[,1] for the code to go thru but Y does not matter for B
  B=MNG[[4]]
  d = ncol(Y) #SS: Y is a matrix [D,etc]
  
  if(gamma_estimator==0){ # dantzig
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,0)
    gamma_hat_list <- vector()
    for (dd in 1:d) {
      gamma_hat<-function(d,z){
        return(b(d,z)%*%beta_hat[, dd])
      }
      gamma_hat_list <- c(gamma_hat_list, gamma_hat)
    }
    
  } else if(gamma_estimator==1){ # lasso
    
    beta_hat=RMD_stable(Y,T,X,p0,D_LB,D_add,max_iter,b,0,1)
    gamma_hat_list <- vector()
    for (dd in 1:d) {
      gamma_hat<-function(d,z){
        return(b(d,z)%*%beta_hat[, dd])
      }
      gamma_hat_list <- c(gamma_hat_list, gamma_hat)
    }
    
  } else if(gamma_estimator==2){ # random forest
    gamma_hat_list <- vector()
    for (dd in 1:d) {
      forest<- do.call(randomForest, append(list(x=B,y=Y[,dd]), arg_Forest))
      gamma_hat<-function(d,z){
        return(predict(forest,newdata=b(d,z), type="response"))
      }
      gamma_hat_list <- c(gamma_hat_list, gamma_hat)
    }
    
  } else if(gamma_estimator==3){ # neural net
    
    # scale down, de-mean, run NN, scale up, remean so that NN works well
    maxs_B <- apply(B, 2, max)
    mins_B <- apply(B, 2, min)
    
    maxs_Y<-apply(Y, 2, max)
    mins_Y<-apply(Y, 2, min)
    
    # hack to ensure that constant covariates do not become NA in the scaling
    const=maxs_B==mins_B
    keep=(1-const)*1:length(const)
    
    NN_B<-B
    NN_B[,keep]<-scale(NN_B[,keep], center = mins_B[keep], scale = maxs_B[keep] - mins_B[keep])

    NN_Y<-Y
    for (dd in 1:d) {
      NN_Y[,dd]<-scale(NN_Y[,dd], center = mins_Y[dd], scale = maxs_Y[dd] - mins_Y[dd])
    }
    gamma_hat_list <- vector()
    
    for (dd in 1:d) {
      nn<- do.call(nnet, append(list(x=NN_B,y=NN_Y[,dd]), arg_Nnet)) #why is it degenerate with fitted.values=1?
      gamma_hat<-function(d,z){
        
        test<-t(as.vector(b2(d,z)))
        NN_b<-test
        NN_b[,keep]<-scale(t(NN_b[,keep]), 
                           center = mins_B[keep], 
                           scale = maxs_B[keep] - mins_B[keep])
        
        NN_Y_hat<-predict(nn,newdata=NN_b)
        Y_hat=NN_Y_hat*(maxs_Y[dd]-mins_Y[dd])+mins_Y[dd]
        
        return(Y_hat)
      }
      gamma_hat_list <- c(gamma_hat_list, gamma_hat)

    }
    
  }
  
  return(list(alpha_hat,gamma_hat_list))
}


