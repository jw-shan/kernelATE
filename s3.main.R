### Figure 2 and Table 3
### in "Kernel estimation of average treatment effects in models with unmeasured confounders"
### edited by ***
### last revised in Nov. 11, 2024
### !!!!! It takes a little longer to run this code !!!!!


rm(list=ls())
library(parallel)
library(reshape2)
library(ggplot2)
library(np)
library(dplyr)
library(tidyr)
source("s3.datagen.R")

# Current time
current_time <- Sys.time()
print(current_time)
formatted_time <- format(current_time, "%Y_%m_%d_%H_%M_%S")


## Monto Carlo setting
seed = 123
J = 1000
numt = c(1000,5000,20000)


# print information
print(paste0("seed=",seed))
print(DataGen)


# Compute truevalue
N=30000000
Data.com<-DataGen(N)
truevalue<-mean(Data.com$delta2)
rm(Data.com)


# parallel setting
# configure parallel environment
ncores = detectCores()
coremax = 120
cl = makeCluster(min(ncores,coremax))
clusterExport(cl,ls())
clusterEvalQ(cl,library(np))->tmp


## Estimation function 
estimation <- function(count) {
  res = matrix(NA, nrow = length(numt), ncol = 3)
  for (i in seq_along(numt)) {
    Data<-DataGen(numt[i],seed+count)
    x<-Data$x
    y<-Data$y
    z<-Data$z
    d<-Data$d

    zhat = fitted(npreg( npregbw(z~x) ))
    y1hat = fitted(npreg( bws=npregbw(y[z==1]~x[z==1]),exdat=x ))
    y0hat = fitted(npreg( bws=npregbw(y[z==0]~x[z==0]),exdat=x ))
    d1hat = fitted(npreg( bws=npregbw(d[z==1]~x[z==1]),exdat=x ))
    d0hat = fitted(npreg( bws=npregbw(d[z==0]~x[z==0]),exdat=x ))

    res[i,1] = mean( (z*y/zhat - (1-z)*y/(1-zhat))/(d1hat-d0hat) ) 
    res[i,2] = mean( (y1hat - y0hat)/(d1hat-d0hat) )
    res[i,3] = mean( (2*z-1)/(z*zhat+(1-z)*(1-zhat))/(d1hat-d0hat)*(y-y0hat-(d-d0hat)*(y1hat-y0hat)/(d1hat-d0hat)) + (y1hat - y0hat)/(d1hat-d0hat) )
  }
  return(res)
}
est  <- parSapply(cl,1:J,estimation)
est  <- t(est)

stopCluster(cl)

est.KIPW.df <- data.frame((est[,seq_along(numt)]-truevalue)%*%diag(sqrt(numt)))
est.KREG.df <- data.frame((est[,-seq_along(numt)][,seq_along(numt)]-truevalue)%*%diag(sqrt(numt)))
est.KMR.df <- data.frame((est[,-seq_along(numt)][,-seq_along(numt)]-truevalue)%*%diag(sqrt(numt)))
est.df.rootn <- cbind(est.KIPW.df,est.KREG.df,est.KMR.df)

combined_df <- bind_rows(
  est.KIPW.df %>% mutate(source = "KIPW"),
  est.KREG.df %>% mutate(source = "KREG"),
  est.KMR.df %>% mutate(source = "KMR")
)
colnames(combined_df)[seq_along(numt)] <- as.character(numt)

df_long <- gather(combined_df, key = "Size", value = "Bias", -source) %>% 
  mutate(source=factor(source, levels = c("KIPW", "KREG", "KMR"))) %>% 
  mutate(Size=factor(Size,levels=as.character(numt)))

# save data
save.image("res_s3.rdata")
# save.image(paste0("res_s3_",formatted_time,".RData"))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  
# ------------plot-----------
# 
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(dplyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(dplyr)
# library(hrbrthemes)
library(viridis)

vioplot_KIPW <- df_long %>% filter(source=="KIPW") %>% 
  ggplot(aes(x = Size, y = Bias, fill = Size)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.2) +
  scale_fill_manual(values = gray.colors(3, start = 1, end = 0.8))+
  theme_bw() + 
  ylim(-70,70)+
  theme(panel.border = element_rect(color = "black", fill = NA),  
        panel.grid = element_blank())+  
  geom_hline(yintercept = 0, color = "red",size=1) + 
  theme(
    plot.title = element_text(size = 11),
    axis.title.y = element_text(hjust = 0.5), 
    axis.title.x = element_text(hjust = 0.5), 
    # axis.text.x = element_blank(),            
    axis.ticks.x = element_blank() ,           
    legend.position = "none"                   
  ) +
  labs(x="Sample size",y="Bias",title="(a) KIPW")+
  scale_x_discrete(labels = numt) 

vioplot_KREG <- df_long %>% filter(source=="KREG") %>% 
  ggplot(aes(x = Size, y = Bias, fill = Size)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.2) +
  scale_fill_manual(values = gray.colors(3, start = 1, end = 0.8))+
  theme_bw() + 
  ylim(-70,70)+
  theme(panel.border = element_rect(color = "black", fill = NA),  
        panel.grid = element_blank())+  
  geom_hline(yintercept = 0, color = "red",size=1) + 
  theme(
    plot.title = element_text(size = 11),
    axis.title.y = element_text(hjust = 0.5), 
    axis.title.x = element_text(hjust = 0.5), 
    # axis.text.x = element_blank(),            
    axis.ticks.x = element_blank() ,           
    legend.position = "none"                   
  ) +
  labs(x="Sample size",y="Bias",title="(b) KREG")+
  scale_x_discrete(labels = numt) 

vioplot_KEIF <- df_long %>% filter(source=="KMR") %>% 
  ggplot(aes(x = Size, y = Bias, fill = Size)) +
  geom_violin(width = 1) +
  geom_boxplot(width = 0.7, color = "black", alpha = 0.2) +
  scale_fill_manual(values = gray.colors(3, start = 1, end = 0.8))+
  theme_bw() + 
  ylim(-70,70)+
  theme(panel.border = element_rect(color = "black", fill = NA),  
        panel.grid = element_blank())+  
  geom_hline(yintercept = 0, color = "red",size=1) + 
  theme(
    plot.title = element_text(size = 11),
    axis.title.y = element_text(hjust = 0.5), 
    axis.title.x = element_text(hjust = 0.5), 
    # axis.text.x = element_blank(),            
    axis.ticks.x = element_blank() ,           
    legend.position = "none"                   
  ) +
  labs(x="Sample size",y="Bias",title="(c) KEIF")+
  scale_x_discrete(labels = numt) 

plt_combined <- grid.arrange( vioplot_KIPW, vioplot_KREG, vioplot_KEIF, nrow = 1, ncol = 3)

ggsave(
  filename = paste0("vioplot_combined.pdf"),
  # filename = paste0("res_s3_",formatted_time,".pdf"),
  plot = plt_combined,
  width = 8.5,             
  height = 3.5,          
  units = "in",
  dpi = 300 
)





# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  
# --------- Summary table -------
#  
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_num = length(numt)

# initaialize 
final.df <- data.frame(samplesize = numt,
                       bias_ipw=NA, stdev_ipw=NA,RMSE_ipw=NA,CR_ipw=NA,
                       bias_reg=NA, stdev_reg=NA,RMSE_reg=NA,CR_reg=NA,
                       bias_mr=NA, stdev_mr=NA,RMSE_mr=NA,CR_mr=NA)


compute_bias <- function(col){
  # col <- col[abs(col)<=150]
  J = length(col)
  Delta <- median(col)
  bias  <- Delta 
  mse   <- 1/J*(sum((col)^2))
  stdev <- sqrt(1/J*(sum((col-Delta)^2)))
  rmse<- sqrt(mse)
  count <- sum(col > - 1.96 * stdev & col <  1.96 * stdev, na.rm = TRUE)
  coverage_rate <- count / J
  CR <- coverage_rate
  res <- cbind(bias,stdev,rmse,CR)
  res <- round(res,3)
  res <- format(res,nsmall = 3)
  return(res)
}

est.bias.df <- apply(est.df.rootn, 2, compute_bias)
final.df[,2:5] <- t(est.bias.df[,1:n_num])
final.df[,6:9] <- t(est.bias.df[,(n_num+1):(2*n_num)])
final.df[,10:13] <- t(est.bias.df[,(2*n_num+1):(3*n_num)])
final.df <- final.df[,-c(4,8,12)]
print(final.df)

## transform to latex
library(knitr)
kable(final.df,format = "latex")

print(Sys.time())
print(Sys.time()-current_time)


# source("s3.main.R")
