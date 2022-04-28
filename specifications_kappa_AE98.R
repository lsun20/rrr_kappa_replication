get_data<-function(df,spec,quintile){
  #outcome <- df[,"educm"]
  outcome <- df[,"ageq2nd"]/4
  n <- length(outcome)
  takeup <- df[,"morekids"]
  # Y <- cbind(takeup,outcome) #SS: for LATE just need Y and participation
  Y <- cbind(takeup,outcome*takeup) #SS: for complier covairate just need D*X and participation
  #SS: create a grid
  # tausAll <- seq(5,95,5)
  # pt <- quantile(outcome, tausAll/100)
  # perc <- unique(pt)
  # unik <- !duplicated(pt)  ## logical vector of unique values
  # taus <- tausAll[unik]/100 ## get unique percentile
  # # taus <- 0 # if not using a grid, set to sth so it goes through
  # Y <- takeup 
  # for (i in 1:length(perc)) {
  #   #SS: create n x d matrix
  #   Y_entry <- ( outcome <= perc[i] )*takeup #grid for covariate grid (D, D*1{X<=x})
  #   # Y_entry <- ( outcome <= perc[i] )*takeup #grid for treated ctnr outcome grid (D, D*1{Y<=y})
  #   # Y_entry <- ( outcome <= perc[i] )*(takeup-1) #grid for untreated ctnr outcome grid (D, (D-1)*1{Y<=y})
  #   Y <- cbind(Y,Y_entry)
  # }
  
  #T <- df[,"multi2nd"]
  T <- df[,"samesex"]
  
  ## low dim specification
  #X.L <- cbind(df[,"educm"], df[,"agem1"], df[,"agefstm"], df[,"boy1st"], df[,"boy2nd"])
  #X.L <- cbind(poly(df[,"educm"], 6, raw=TRUE))
  X.L <- cbind(poly(df[,"ageq2nd"], 6, raw=TRUE))
  ## high dim specification. NOTE: original paper is this spec squared (pairwise interactions thereof?)
  X.H <- cbind(poly(df[,"agem1"], 6, raw=TRUE),
               poly(df[,"agefstm"], 8, raw=TRUE),
               df[,"boy1st"], df[,"boy2nd"]) 
  # what does original line of code mean?
  # copied from EJ: "(poly(age, 6, raw=TRUE) + poly(inc, 8, raw=TRUE) + poly(educ, 4, raw=TRUE) + poly(fsize, 2, raw=TRUE) + marr + twoearn + db + pira + hown)^2"
  X.vH=model.matrix(~(poly(df[,"agem1"], 6, raw=TRUE) + 
                        poly(df[,"agefstm"], 8, raw=TRUE) + 
                        poly(df[,"educm"], 4, raw=TRUE) + 
                        poly(df[,"ageq2nd"], 6, raw=TRUE) +
                        df[,"boy1st"] + 
                        df[,"boy2nd"])^2)
  X.vH=X.vH[,-1]
  
  if (spec==1){
    X=X.L
  } else if (spec==2) {
    X=X.H
  } else {
    X=X.vH
  }
  
  X <- scale(X,center=FALSE,scale=TRUE)
  
  #impose common support
  p.1 <- multinom(T~X-1, trace=FALSE)$fitted.values
  #indexes.to.drop <- which(p.1 < min(p.1[T==1]) | max(p.1[T==1]) < p.1)
  indexes.to.drop <- which(p.1 < 0 | 1 < p.1) # trim no obs
  if (length(indexes.to.drop)==0) {indexes.to.drop <- n+1}	#R throws a wobbly if [-indexes.to.drop] is negating an empty set. 
  n.per.treatment <- as.vector(table(T[-indexes.to.drop]))
  n.trim <- n.per.treatment[1]+n.per.treatment[2]
  
  Y.trimmed=Y[-indexes.to.drop,]
  T.trimmed=T[-indexes.to.drop]
  X.trimmed=X[-indexes.to.drop,]
  
  #keep this group: 
  #inc_low - bottom quintile of income
  #inc_mid - middle quintile of income
  #inc_high - top quintile of income
  
  # if (spec==1){
  #   inc=X.trimmed[,2]
  # } else if (spec==2) {
  #   inc=X.trimmed[,7]
  # } else {
  #   inc=X.trimmed[,7]
  # }
  # 
  if (quintile>0){
    q <- ntile(inc, 5)
    Y.q=Y.trimmed[q==quintile,]
    T.q=T.trimmed[q==quintile]
    X.q=X.trimmed[q==quintile,]
  } else {
    Y.q=Y.trimmed
    T.q=T.trimmed
    X.q=X.trimmed
  }
  
  return(list(Y.q,T.q,X.q,perc))
  
}


