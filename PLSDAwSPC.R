### PLS-DA on selective pivot coordinates (SPCs)

require("compositions")
require("doParallel")
require("pls")

####################################################################################################################################################

PLSDAonSPC = function(X, N1, xi = 0.1){
  # INPUT: X - object of class data.frame or matrix of size (N x D); positive values only; 
  #            ordered so that the first N1 rows correspond to the one group (controls) and the remaining N-N1 to the second group (patients)
  #        N1 - the number of the observations in the first group;
  #        xi - adjustable value in the setting of the selective interval, default is 0.1;
  # OUTPUT: BETA - matrix of size (K x D) with bootstrap PLS model coefficients (K being the number of bootstrap replicates);
  #         signif - vector of length D indicating the significance of the variables 
  #                  (-1 for significantly decreased in patients, 0 for non-significant, 1 for significantly increased in patients);
  #         SL - output from "selectLR" function

  D = ncol(X)
  N = nrow(X)
  cn = colnames(X)
  y = c(rep(0, N1), rep(1, N-N1)) # vector indicating the groups (0 for controls, 1 for patients)
  yc = as.vector(scale(y, scale = F)) 
  Xc = as.matrix(as.data.frame(acomp(X)-mean(acomp(X))))
  codes = matrix(rep(0, (D-1)*D), ncol=D)
  for(ico in 1:(D-1)){
    codes[ico,] = c(rep(0, ico-1), 1, rep(-1, D-ico))
  }   
  V = gsi.buildilrBase(t(codes)) # matrix of logcontrast coefficients
  Z = log(Xc)%*%V # ordinary pivot coordinates
  misclass = mvr(yc~Z, ncomp = 10, validation="CV")
  comp = selectNcomp(misclass) # optimal number of PLS components based on the randomization test approach
  SL = selectLR(X, N1, xi) # selecting logratios for the construction of SPCs
  I = SL$I
  pos = SL$pos
  no_cores = detectCores() # the number of cores
  no_cores = no_cores-1 # the number of cores used
  cl = makeCluster(no_cores)
  registerDoParallel(cl)
  # bootstrap
  K = 500 # the number of bootstrap replicates
  BETA = foreach(k=1:K,
                 .packages = c("compositions", "pls"),
                 .combine=rbind) %dopar% {
                   ind1 = sample(1:N1, N1, replace = T)  
                   ind2 = sample((N1+1):N, N-N1, replace = T)
                   ind = c(ind1, ind2)
                   testX = X[ind, ]
                   testXc = as.matrix(as.data.frame(acomp(testX)-mean(acomp(testX))))
                   beta = rep(NA, D)
                   for(i in 1:D){
                     I1 = I[i, ]
                     pos1 = pos[i]
                     testXc1 = testXc[, I1]
                     testZ = log(testXc1)%*%V # i-th system of SPCs
                     tplsda = mvr(yc~testZ, ncomp = comp)
                     beta1 = as.vector(coef(tplsda))
                     beta[i] = beta1[pos1] # regression coeff. corresponding to the SPC of interest
                   }
                  return(beta)
                }
  stopCluster(cl)
  betaM = apply(BETA, 2, mean)
  betaS = apply(BETA, 2, sd)
  stdC = betaM/betaS # standardized PLS model coefficients
  names(stdC) = cn
  pval = 2*pnorm(abs(stdC), lower.tail = F)
  b_pval = p.adjust (pval, method = "BH")
  alpha0 = 0.05
  signif = rep(0, D)
  for(i in 1:D){
    if((b_pval[i] < alpha0) & (stdC[i] >= 0)){signif[i] = 1 # potential biomarker increased in patients
    } else if ((b_pval[i] < alpha0) & (stdC[i] < 0)){signif[i] = -1} # potential biomarker decreased in patients
  }
  names(signif) = cn
  return(list("BETA" = BETA, "signif" = signif, "SL" = SL))
}