### Table 1
### in "Kernel estimation of average treatment effects in models with unmeasured confounders"
### edited by Jiawei Shan (jwshan@ruc.edu.cn)
### last revised in Nov. 11, 2024

rm(list=ls())
library(parallel)
library(reshape2)
source("s1.datagen.R")
source("s1.estimators.R")


## Monto Carlo times and Sample Size###
seed = 11
J <- 500
N <- 500
truevalue<- 0.087

# parallel setting
# configure parallel environment
ncores = detectCores()
if (ncores<40) {
  cl = makeCluster(ncores)
}else{
  cl = makeCluster(40)
}
clusterExport(cl,ls())

## Estimation function 
estimation <- function(count) {

  Data<-DataGen(N,seed*count)
  hopt <- 1.06*sd(Data$x)* N^{-1/5}
  h <- hopt * N^{1/5} * N^{-2/7}
  X<-Data$x
  Z<-Data$z
  D<-Data$d
  Y<-Data$y
  
  Naive <- mean(Y[D==1])-mean(Y[D==0])
  KIPW <- KSE_1(X,Y,D,Z,h)
  KREG <- KSE_3(X,Y,D,Z,h)
  KMR  <- KSE_t(X,Y,D,Z,hopt)
  est <- cbind(Naive,KIPW,KREG,KMR)
  
  return(est)
}

est  <- parSapply(cl,1:J,estimation)
est  <- t(est)
colnames(est)<-c("Naive","KIPW","KREG","KMR")

result <- matrix(nrow = 4, ncol = 4)
colnames(result)<-c("bias","stdev","RMSE","CR")
rownames(result)<-c("Naive","KIPW","KREG","KMR")
for (i in 1:4) {
  Delta <- mean(est[,i])
  bias  <- Delta - truevalue
  mse   <- 1/J*(sum((est[,i]-truevalue)^2))
  stdev <- sqrt(1/J*(sum((est[,i]-Delta)^2)))
  rmse<- sqrt(mse)
  
  count<-0
  for(j in 1:J){
    if(est[j,i]> 0.087-1.96*stdev & est[j,i]< 0.087+1.96*stdev)
      count<- count+1
  }
  coverage_rate <<- count/J
  CR <- coverage_rate
  
  result[i,] <- cbind(bias,stdev,rmse,CR)
}

print(result)

# save data
save.image("res_s1.rdata")

stopCluster(cl)




