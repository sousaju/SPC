### Selecting logratios for the construction of selective pivot coordinates (SPCs)

require("robustbase")

####################################################################################################################################################

selectLR = function(X, N1, xi = 0.1){
  # INPUT: X - object of class data.frame or matrix of size (N x D); positive values only; 
  #            ordered so that the first N1 rows correspond to the one group (controls) and the remaining N-N1 to the second group (patients);
  #        N1 - the number of the observations in the first group;
  #        xi - adjustable value in the setting of the selective interval, default is 0.1;
  # OUTPUT: I - matrix with D rows and D columns; 
  #             each row contains indices of columns which indicates how the columns should be reordered for the construction of the respective SPCs;
  #         pos - vector of length D indicating the position of the SPCs of interest;
  #         Ts - matrix with D rows and D columns containing the Welch's t-statistics for the pairwise logratios 
  #              (row corresponds to the part in the numerator, column to the part in the denominator);
  #         W - matrix with D rows and D columns indicating if the respective logratio is included into the aggregation (1) or not (0);
  
  D = ncol(X)
  N = nrow(X)
  q = qt(1-0.05/2, N-2)
  I = matrix(0, D, D)
  pos = rep(0, D)
  Ts = matrix(0, D, D)
  W = matrix(1, D, D)
  diag(W) = NA
  for(i in 1:D){
    lr = log(X/X[,i])
    tstat = c()
    for(j in 1:D){
      if(j!=i){
        g1 = lr[1:N1, j]
        g2 = lr[(N1+1):N, j]
        tstat = c(tstat, t.test(g1,g2)$statistic)
      }
    }
    md = median(tstat)
    qn = Qn(tstat)
    theta1 = md-2*qn
    theta2 = md+2*qn
    if (quantile(tstat, xi) > q) {
      theta2 = max(tstat) 
    } else if (quantile(tstat, 1-xi) < -q) {
      theta1 = min(tstat)
    }
    tstat0 = c(tstat[1:(i-1)], NA, tstat[(i):(D-1)])
    if (i==1) {
      tstat0 = c(NA, tstat)
    } else if (i==D) {
      tstat0 = c(tstat, NA)
    }
    Ts[i, ] = tstat0 
    i0 = c(which(tstat0<theta1), which(tstat0>theta2))
    W[i, i0] = 0
    pos[i] = length(i0)+1
    i1 = c(1:D)[-c(i, i0)]
    I[i, ] = c(i0, i, i1)
  }
  return(list("I" = I, "pos" = pos, "Ts" = Ts, "W" = W))
}