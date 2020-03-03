L=5

rrr<-function(Y,T,X,p0,D_LB,D_add,max_iter,b,alpha_estimator,gamma_estimator){
  
  n=nrow(X)
  folds <- split(sample(n, n,replace=FALSE), as.factor(1:L))

  # theta=rep(0,L)
  # v_sq=rep(0,L)
  
  d = ncol(Y) #SS: Y is a matrix [D,etc]
  theta = matrix(0,L,d) #SS: might not need this after all...
  moments_full = matrix(0,0,d) #SS: store the Psi's for each obs
  for (l in 1:L){
    
    Y.l=Y[folds[[l]],]
    Y.nl=Y[-folds[[l]],]
    
    T.l=T[folds[[l]]]
    T.nl=T[-folds[[l]]]
    
    X.l=X[folds[[l]],]
    X.nl=X[-folds[[l]],]
    
    n.l=length(T.l)
    n.nl=length(T.nl)
    moments = matrix(0,n.l,d) #SS: store the Psi's for each fold
    
    # get stage 1 (on nl)
    stage1_estimators<-get_stage1(Y.nl,T.nl,X.nl,p0,D_LB,D_add,max_iter,b,alpha_estimator,gamma_estimator)

    
    # debugging
    print(paste0('fold: ',l))
    #print(paste0('beta_hat: '))
    #print(paste0(round(beta_hat,2)))
    #print(paste0('rho_hat: '))
    #print(paste0(round(rho_hat,2)))
    
    # parameter to linear function
    alpha_hat=stage1_estimators[[1]]
    
    
    for (dd in 1:d) {
      gamma_hat=stage1_estimators[[2]][[dd]]
      
      #get stage 2 (on l)
      Psi=rep(0,n.l)
      for (i in 1:n.l){
        Psi[i]=psi(Y.l[i,dd],T.l[i],X.l[i,],m,alpha_hat,gamma_hat)
      }
      moments[, dd] = Psi
      theta[l,dd]=mean(Psi)
    }
    moments_full = rbind(moments_full, moments)
  }
  #SS: denominator is the first stage, corr to first col of Y
  #SS: numerator is the reduced form, corr to the rest col of Y
  j_hat <- mean(moments_full[,1])
  theta_hat <- colMeans(moments_full[,2:d])/j_hat
  psi_hat <- moments_full[,2:d] - moments_full[,1] %*% t(theta_hat)
  omega_hat <- t(psi_hat) %*% psi_hat  / n
  
  return(list(theta_hat,j_hat,omega_hat))
  
}

printer<-function(spec1){
  print(paste(" treated: ",spec1[1], " untreated: ", spec1[2], "   ATE:    ",round(spec1[3],2), "   SE:   ", round(spec1[4],2), sep=""))
}

for_tex<-function(spec1){
  print(paste(" & ",spec1[1], " & ", spec1[2], "   &    ",round(spec1[3],2), "   &   ", round(spec1[4],2), sep=""))
}