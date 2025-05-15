logit <- function(prob){ log(prob) - log(1-prob)}
expit <- function(logodds){ 1/(1+exp(-logodds))}

DataGen = function(n,SEED=1){
  
  set.seed(SEED)
  
  x = runif(n, -1, 1)
  ps = pnorm(x)
  z = rbinom(n, 1, ps)
  u = rbinom(n, 1, 0.5)
  
  y= z/ps - (1-z)/(1-ps) + 0.3*(2*u-1)
  lg = pnorm(1.5*x)*0.5 + 0.4
  d.ps=lg *z + (1-z)*0.05  + 0.05*(2*u-1)
  d = rbinom(n,1,d.ps)
  delta2 = (z/ps-(1-z)/(1-ps))*y/(lg-0.05)
  
  return(list(x=x,ps=ps,z=z,y=y,d=d,d.ps=d.ps,delta2=delta2))
}

