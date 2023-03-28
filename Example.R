### Example: PLS-DA on selective pivot coordinates (SPCs)

# Generating the data
gendata = function(N1 = 20, D = 500, R = 20,
                   umin = 1, umax = 100, vmin = 1, vmax = 10,
                   sig_k = 0.3, sig_h = 0.05, 
                   sig_f = 0.2, sig_g = 0.007, A = 1){ 
  N = 2*N1
  k = rnorm(N, 0, sig_k)%*%t(rep(1, D))
  u = runif(D, min = umin, max = umax)
  v = runif(D, min = vmin, max = vmax)
  conc = u/v
  v = rep(1, N)%*%t(v)
  C = rep(1, N)%*%t(conc)
  C[1:N1, 1:R] = C[1:N1, 1:R]+matrix(1, N1, R)*A
  f = matrix(rnorm(N*D, 0, sig_f), N, D)
  g = matrix(rnorm(N*D, 0, sig_g), N, D)
  h = matrix(rnorm(N*D, 0, sig_h), N, D)
  X = (1-k)*(C+f)*v*exp(g)+h
  return(X)
}

set.seed(1)
X = gendata(D = 200) # simulated data 
                     # the first 20 rows correspond to patients, the remaining 20 to controls
                     # the first 20 columns correspond to biomarkers increased in patients, the remaining 180 to non-biomarkers
colnames(X) = paste("V", 1:200, sep = "")
X = X[c(21:40, 1:20), ] # reordering the rows so that controls are first, then patients



# Performing PLS-DA on SPCs
source("selectLR.R")
source("PLSDAwSPC.R")

plsda = PLSDAonSPC(X = X, N1 = 20)

sig = plsda$signif

colnames(X)[sig == 1] # identified as potential biomarkers increased in patients
colnames(X)[sig == -1] # identified as potential biomarkers decreased in patients



# Visualizing the selected logratios
library(gplots)
library(RColorBrewer)

SL = plsda$SL
Ts = SL$Ts
colnames(Ts) = rownames(Ts) = colnames(X)
W = SL$W
WW = as.data.frame(W)
WW[WW==1] = NA
WW[WW==0] = "*"
RowMedians = rowMedians(Ts, na.rm = T)
ord = order(-RowMedians)
TsO = Ts[ord, ord]
WWO = WW[ord, ord]

pdf("tstat.pdf", width = 14, height = 14.5)
heatmap.2(TsO,
          col = brewer.pal(11,"BrBG"),
          trace = "none",
          dendrogram = "none",
          Rowv = F, Colv = F, 
          cexRow = 0.5, cexCol = 0.5, 
          sepcolor = "black",
          key.xlab = "Welch's t-statistics", 
          key.title = NA,
          cellnote = WWO,
          notecol = "black",
          lmat = matrix(c(0,0,2,4,3,1), 3, 2),
          lwid=c(0.01,0.4), lhei=c(0.8,0.01,3))
dev.off()
