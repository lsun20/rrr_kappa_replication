naive_kappa <- function(Y, T, X) {
  n = nrow(X)
  d = ncol(Y) #SS: Y is a matrix [D,etc]
  
  theta = matrix(0, d) #SS: might not need this after all...
  moments_full = matrix(0, n, d) #SS: store the Psi's for each obs
  
  #kappa reweighting via logistic regression
  pi_hat <- predict(glm(T ~ X, family = binomial), type = "response")
  
  # censor pi_hat
  # pi_hat[pi_hat > 1 - 10^-12] <- 1 - 10^-12;
  # pi_hat[pi_hat < 10^-12] <-  10^-12;
  
  alpha_hat = (T - pi_hat) / (pi_hat * (1 - pi_hat))
  
  for (dd in 1:d) {
    moments_full[, dd] <- Y[, dd] * alpha_hat
  }
  
  #create an index for trimming
  # keep = (pi_hat <= 0.9 & pi_hat >= 0.1 );
  keep = (pi_hat <= 1 - 10^-12 & pi_hat >= 10^-12 );
  moments_full <- moments_full[keep,]
  
  j_hat <- mean(moments_full[, 1])
  theta_hat <- colMeans(moments_full[, 2:d]) / j_hat
  return(list(theta_hat, j_hat))
}
