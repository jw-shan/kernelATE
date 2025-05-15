source("s1.func.R")
source("s1.truevalue.R")



DataGen = function(n,SEED){
  
  set.seed(SEED)
  
  # x: baseline covariates. 
  x1           = rep(1,n)  #intercept term
  x2.ind       = rbinom(n,1,0.5)
  x2           = ((-1)^x2.ind) * runif(n,0.5,1)
  x            = cbind(x1,x2)
  u            = rbinom(n, 1, 0.5)
  
  ### True values depend on data
  pix.true     = expit(x %*% gamma.true)
  delta.true   = tanh(x %*% alpha.true)
  delta.d.true = tanh(x %*% beta.true)
  delta.y.true = delta.true * delta.d.true
  logOP.y.true = x %*% zeta.true 
  logOP.d.true = x %*% eta.true 
  
  p0p1.y.true = mapply(getProbScalarRDiff,atanh(delta.y.true),logOP.y.true)
  p0p1.d.true = mapply(getProbScalarRDiff,atanh(delta.d.true),logOP.d.true)
  
  z = rbinom(n,1,pix.true)
  
  p.d.true = p0p1.d.true[1,]
  p.d.true[z==1] = p0p1.d.true[2,z==1]
  p.d.true = p.d.true + 0.1*(2*u-1)
  d = rbinom(n,1,p.d.true)
  
  p.y.true = p0p1.y.true[1,]
  p.y.true[z==1] = p0p1.y.true[2,z==1]
  p.y.true = p.y.true + 0.1*(2*u-1)
  y = rbinom(n,1,p.y.true)
  p.d.true = p0p1.d.true[1,]
  p.d.true[z==1] = p0p1.d.true[2,z==1]
  p.d.true = p.d.true + 0.1*(2*u-1)
  d = rbinom(n,1,p.d.true)
  
  p.y.true = p0p1.y.true[1,]
  p.y.true[z==1] = p0p1.y.true[2,z==1]
  p.y.true = p.y.true + 0.1*(2*u-1)
  y = rbinom(n,1,p.y.true)
  
  return(list(x=x[,2],z=z,d=d,y=y,delta.true=delta.true))
  
}
