fy <- function(z,x) {
  return(0 + 2*x^2*z)
}
fd <- function(z,x) {
  return(0 + 1*x*z)
}
px <- function(x) {
  return(0.05*(x<0.5) + 0.95*(x>0.5))
}
b_x <- function(x) {
  return(cbind(x,x^2,x^3,x^4)) # later will add constant term
}
get_sim_data<-function(N,s){
  set.seed(s) # for replication
  x <- runif(N)
  z <- runif(N) < px(x)
  d <- runif(N) < fd(z,x)
  y <- rnorm(N) + fy(z,x)
  
  outcome <- y
  takeup <- d
  Y <- cbind(outcome, takeup) #SS: for LATE just need Y and participation
  #SS: create a grid

  # perc <- c(-1,0,1,2,3,4) #SS: grid1
  perc <- c(-2,-1,0,1,2,3) #SS: grid0
  
  Y <- takeup 
  for (i in 1:length(perc)) {
    #SS: create n x d matrix
    # Y_entry <- ( outcome <= perc[i] )*takeup #grid for treated ctnr outcome grid (D, D*1{Y<=y})
    Y_entry <- ( outcome <= perc[i] )*(takeup-1) #grid for untreated ctnr outcome grid (D, (D-1)*1{Y<=y})
    Y <- cbind(Y,Y_entry)
  }
  
  T <- z
  
  ## low dim specification
  X <- b_x(x)
  ## high dim specification. NOTE: original paper is this spec squared (pairwise interactions thereof?)
  
  return(list(Y,T,X,perc))
  
}