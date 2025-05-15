# This function implements cross-validation using "np" packages.

KSE_CV <- function(X,Y,D,Z,bwt=1){
  n = length(Z)
  
  X1 = X[Z==1]
  X0 = X[Z==0]
  Y1 = Y[Z==1]
  Y0 = Y[Z==0]
  D1 = D[Z==1]
  D0 = D[Z==0]
  
  # f(Z|X)
  fhat_bw = npregbw(Z~X)$bw
  pi.hat = fitted(npreg(tydat=Z ,txdat=X, exdat=X, bws=fhat_bw*bwt))
  f.hat = Z*pi.hat + (1-Z)*(1-pi.hat)
  
  # E(D|Z,X)
  pD1.hat_bw = npregbw(D1~X1)$bw
  pD1.hat = fitted(npreg(tydat=D1 ,txdat=X1, exdat=X, bws=pD1.hat_bw*bwt))
  
  pD0.hat_bw = npregbw(D0~X0)$bw
  pD0.hat = fitted(npreg(tydat=D0 ,txdat=X0, exdat=X, bws=pD0.hat_bw*bwt))
  
  # E(Y|Z,X)
  pY1.hat_bw = npregbw(Y1~X1)$bw
  pY1.hat = fitted(npreg(tydat=Y1 ,txdat=X1, exdat=X, bws=pY1.hat_bw*bwt))
  
  pY0.hat_bw = npregbw(Y0~X0)$bw
  pY0.hat = fitted(npreg(tydat=Y0 ,txdat=X0, exdat=X, bws=pY0.hat_bw*bwt))
  
  deltaD.hat = pD1.hat - pD0.hat
  deltaY.hat = pY1.hat - pY0.hat
  
  # estimator
  IPW = (2*Z-1)/f.hat/deltaD.hat*Y
  IPW = mean(IPW)
  
  REG = deltaY.hat/deltaD.hat
  REG = mean(REG)
  
  MR = (2*Z-1)/f.hat/deltaD.hat*
    (Y- pY0.hat - (D-pD0.hat)*deltaY.hat/deltaD.hat ) +  deltaY.hat/deltaD.hat
  MR = mean(MR)
  
  res = list(IPW=IPW,REG=REG,MR=MR)
  
  return(res)
  
  
  
  
  
  
}