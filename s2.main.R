### Figure 1 and Table 2
### in "Kernel estimation of average treatment effects in models with unmeasured confounders"
### edited by ***
### last revised in Nov. 11, 2024


rm(list=ls())
library(parallel)
library(reshape2)
library(ggplot2)
library(np)
library(dplyr)
library(tidyr)



source("s2.datagen.R")
source("s2.cv.R")


# Current time
current_time <- Sys.time()
print(current_time)
formatted_time <- format(current_time, "%Y_%m_%d_%H_%M_%S")


## Monto Carlo times and Sample Size###
seed = 17
J <- 500
N <- 1000
truevalue<- 0.087
bwt = c(0.2,0.5,0.8,1,1.5,2,3,4)

# bw name
bw_name = vector()
for (i in seq_along(bwt)) {
  if (bwt[i] != 1) {
    bw_name[i] <- paste0(bwt[i], "h_cv")  # names
  }else{
    bw_name[i] = "h_cv"
  }
}
bw_name = factor(bw_name,levels = bw_name)


# parallel setting
# configure parallel environment
ncores = detectCores()
coremax = 32
if (ncores<coremax) {
  cl = makeCluster(ncores)
}else{
  cl = makeCluster(coremax)
}
clusterExport(cl,ls())
clusterEvalQ(cl,library(np))


## Estimation function 
estimation <- function(count) {
  
  Data<-DataGen(N,seed*count)
  X<-Data$x
  Z<-Data$z
  D<-Data$d
  Y<-Data$y
  
  # est.list  <- KSE_CV(X,Y,D,Z)
  
  res = matrix(NA, nrow = length(bwt), ncol = 3)
  for (i in seq_along(bwt)) {
    est.list  <- KSE_CV(X,Y,D,Z,bwt[i])
    KIPW = est.list$IPW
    KREG = est.list$REG
    KMR  = est.list$MR
    res[i,] = c(KIPW,KREG,KMR)
  }
  
  return(res)
}

est  <- parSapply(cl,1:J,estimation)
est  <- t(est)


stopCluster(cl)



# 

est.df <- data.frame(est)
est.KIPW.df <- data.frame(est[,seq_along(bwt)])
est.KREG.df <- data.frame(est[,-seq_along(bwt)][,seq_along(bwt)])
est.KMR.df <- data.frame(est[,-seq_along(bwt)][,-seq_along(bwt)])



# combine data
combined_df <- bind_rows(
  est.KIPW.df %>% mutate(source = "KIPW"),
  est.KREG.df %>% mutate(source = "KREG"),
  est.KMR.df %>% mutate(source = "KMR")
)

df_long <- gather(combined_df, key = "Bandwidth", value = "Estimates", -source)
df_long$source <- factor(df_long$source, levels = c("KIPW", "KREG", "KMR"))

# save data
save.image("res_s2.rdata")



# # -------------plot----------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(grid)
library(gridExtra)
library(latex2exp)

bw_name = bwt

## plot three boxplots
boxplot_KIPW <- ggplot(subset(df_long,source=="KIPW"), aes(x = Bandwidth, y = Estimates)) +
  geom_boxplot(outlier.alpha = 0,fill="lightgray") +
  # geom_crossbar(width = 0.3, fatten = 0.5) +
  labs(x = TeX("Bandwidths ($\\times h_{cv}$)"), y = "Estimates",title = "(a) KIPW") +
  ylim(-0.8,1.2)+
  geom_hline(yintercept = truevalue, color = "red", linetype=1)+
  scale_x_discrete(labels = bw_name)+
  theme_bw() + # 使用白色背景
  theme(panel.border = element_rect(color = "black", fill = NA),  # 设置四周边框为黑色
        panel.grid = element_blank())  # 隐藏网格线
boxplot_KIPW

boxplot_KREG <- ggplot(subset(df_long,source=="KREG"), aes(x = Bandwidth, y = Estimates)) +
  geom_boxplot(outlier.alpha = 0,fill="lightgray") +
  # geom_crossbar(width = 0.3, fatten = 0.5) +
  labs(x = TeX("Bandwidths ($\\times h_{cv}$)"), y = "Estimates",title = "(b) KREG") +
  ylim(-0.8,1.2)+
  geom_hline(yintercept = truevalue, color = "red", linetype=1)+
  scale_x_discrete(labels = bw_name)+
  theme_bw() + # 使用白色背景
  theme(panel.border = element_rect(color = "black", fill = NA),  # 设置四周边框为黑色
        panel.grid = element_blank())  # 隐藏网格线
boxplot_KREG

boxplot_KMR <- ggplot(subset(df_long,source=="KMR"), aes(x = Bandwidth, y = Estimates)) +
  geom_boxplot(outlier.alpha = 0,fill="lightgray") +
  # geom_crossbar(width = 0.3, fatten = 0.5) +
  labs(x = TeX("Bandwidths ($\\times h_{cv}$)"), y = "Estimates",title = "(c) KEIF") +
  ylim(-0.8,1.2)+
  geom_hline(yintercept = truevalue, color = "red", linetype=1)+
  scale_x_discrete(labels = bw_name)+
  theme_bw() + # 使用白色背景
  theme(panel.border = element_rect(color = "black", fill = NA),  # 设置四周边框为黑色
        panel.grid = element_blank())  # 隐藏网格线
boxplot_KMR


# combine
plt_combined <- grid.arrange( boxplot_KIPW, boxplot_KREG, boxplot_KMR, nrow = 1, ncol = 3)

ggsave(
  filename = paste0("plt_combined.pdf"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  plot = plt_combined,
  width = 8.5,             # 宽
  height = 2.5,            # 高
  units = "in",          # 单位
  dpi = 300              # 分辨率DPI
)


# -------- tables ----------- 
n_bw = length(bwt)

# initaialize res
final.df <- data.frame(bandwidth = bw_name,
                       bias_ipw=NA, stdev_ipw=NA,RMSE_ipw=NA,CR_ipw=NA,
                       bias_reg=NA, stdev_reg=NA,RMSE_reg=NA,CR_reg=NA,
                       bias_mr=NA, stdev_mr=NA,RMSE_mr=NA,CR_mr=NA)


result <- matrix(nrow = 3, ncol = 4)
colnames(result)<-c("bias","stdev","RMSE","CR")
rownames(result)<-c("IPW","REG","MR")

compute_bias <- function(col){
  threshold <- 2
  col = col[abs(col)<=threshold] 
  J = length(col)
  Delta <- mean(col)
  bias  <- Delta - truevalue
  mse   <- 1/J*(sum((col-truevalue)^2))
  stdev <- sqrt(1/J*(sum((col-Delta)^2)))
  rmse<- sqrt(mse)
  count <- sum(col > Delta - 1.96 * stdev & col < Delta + 1.96 * stdev, na.rm = TRUE)  
  coverage_rate <- count / J
  CR <- coverage_rate
  res <- cbind(bias,stdev,rmse,CR)
  res <- round(res,3)
  res <- format(res,nsmall = 3)
  return(res)
}

est.bias.df <- apply(est.df, 2, compute_bias)
final.df[,2:5] <- t(est.bias.df[,1:n_bw])
final.df[,6:9] <- t(est.bias.df[,(n_bw+1):(2*n_bw)])
final.df[,10:13] <- t(est.bias.df[,(2*n_bw+1):(3*n_bw)])

print(final.df)
