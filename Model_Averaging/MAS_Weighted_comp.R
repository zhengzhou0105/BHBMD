# Settings-------
# Settings-------
 mywd <- "/N/scratch/zhezhou/MAS_Expo5/"
# mywd <- paste0(getwd(),"/")
fileDir <- "/N/slate/zhezhou/"
# seed <- 47405      # use IUB postal code as seed for random simulation
DSName <- "Expo5"
batch <- 100
iter <- 5
Napproach <- 7
ApprList <- c("AIC","BIC_BMA","WAIC","LPML","Stacking","PBMA","PBMABB")

samplingiter <- 1E4*5   # number of iterations in stan sampling
warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 3       # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling
# setwd(mywd)
options(bitmapType='cairo')


source(file = paste0(fileDir,"BHBMD_model_utilities.R"),
       local = T)
batchsize <- iter
iterTotal <- batch * batchsize

# Load batch files------

weighted_batch <- vector("list",iterTotal)

for(b in 1:batch){
  for(i in 1:iter){
    index <- batchsize * (b-1)+i
    weighted_temp <- read.csv(file = paste0("Weighted_batch_",DSName,"_batch_",b,"_",i,".csv"))
    weighted_temp <- weighted_temp %>% select(-1)
    Nobs <- ncol(weighted_temp) / Napproach
    Xnames <- colnames(weighted_temp)[1:Nobs]
    Xvalue <- as.numeric(stringr::str_replace(Xnames,pattern = "X",""))
    ## transform into wide format
    weighted_wide_temp <- apply(weighted_temp,MARGIN = 1,function(x){
      data.frame(
        X = rep(Xvalue,times=Napproach),
        Appr = rep(ApprList,each = Nobs),
        Wtd = as.numeric(x)
      )
    })
    weighted_wide <- do.call(rbind,weighted_wide_temp) %>% mutate(
      iter = rep(1:iterEff,each = Nobs * Napproach),
    )

    weighted_batch[[index]] <- weighted_wide
    rm(weighted_temp)
    rm(weighted_wide)
    rm(weighted_wide_temp)
  }
}

Yraw <- vector("list",iterTotal)
# DataList_batch <- vector("list",iterTotal)
for(b in 1:batch){
  for(i in 1:iter){
    index <- batchsize * (b-1)+i
    Data_temp <- read.csv(file = paste0("DataList_",DSName,"_batch_",b,"_",i,".csv"))
    # DataList_batch[[index]] <- Data_temp
    Yraw_temp <- Data_temp$ymean
    Yraw[[index]] <- Yraw_temp
    rm(Yraw_temp)
    rm(Data_temp)
  }
}

# Distributional summary------


# # Distributional summary over iterations (of simulation)
# 
# # median
# median_wide_batch <- lapply(weighted_batch,function(x){
#   x %>% group_by(Appr,X) %>% summarise(median(Wtd))
# })
# Q5_wide_batch <- lapply(weighted_batch,function(x){
#   x %>% group_by(Appr,X) %>% summarise(quantile(Wtd,0.05))
# })
# Q95_wide_batch <- lapply(weighted_batch,function(x){
#   x %>% group_by(Appr,X) %>% summarise(quantile(Wtd,0.95))
# })
# 
# median_batch <- lapply(weighted_batch,function(x){
#   apply(x,MARGIN = 2,median)
# })
# 
# summary_median <- lapply(
#   median_batch,function(x){
#     temp <- as.data.frame(x)
#     Nobs <- nrow(temp) / Napproach        # number of data points
#     Xnames <- rownames(temp)[1:Nobs]      # X values
#     # reorganize into wide format
#     temp <- temp %>% mutate(
#       X = rep(Xnames,times = Napproach),
#       Appr = rep(1:Napproach,each = Nobs)
#     ) %>% rename(Median = x)
#     wide <- dcast(data = temp,
#                          X ~ Appr, value.var = "Median")
#     colnames(wide)[2:(Napproach+1)] <- ApprList
#     
#     return(wide)
#   }
# )
# 
# # 5th percentile
# Q5_batch <- lapply(weighted_batch,function(x){
#   apply(x,MARGIN = 2,quantile,0.05)
# })
# summary_Q5 <- lapply(
#   Q5_batch,function(x){
#     temp <- as.data.frame(x)
#     Nobs <- nrow(temp) / Napproach        # number of data points
#     Xnames <- rownames(temp)[1:Nobs]      # X values
#     # reorganize into wide format
#     temp <- temp %>% mutate(
#       X = rep(Xnames,times = Napproach),
#       Appr = rep(1:Napproach,each = Nobs)
#     ) %>% rename(Q5 = x)
#     wide <- dcast(data = temp,
#                   X ~ Appr, value.var = "Q5")
#     colnames(wide)[2:(Napproach+1)] <- ApprList
#     
#     return(wide)
#   }
# )
# 
# # 95th percentile
# Q95_batch <- lapply(weighted_batch,function(x){
#   apply(x,MARGIN = 2,quantile,0.95)
# })
# summary_Q95 <- lapply(
#   Q95_batch,function(x){
#     temp <- as.data.frame(x)
#     Nobs <- nrow(temp) / Napproach        # number of data points
#     Xnames <- rownames(temp)[1:Nobs]      # X values
#     # reorganize into wide format
#     temp <- temp %>% mutate(
#       X = rep(Xnames,times = Napproach),
#       Appr = rep(1:Napproach,each = Nobs)
#     ) %>% rename(Q95 = x)
#     wide <- dcast(data = temp,
#                   X ~ Appr, value.var = "Q95")
#     colnames(wide)[2:(Napproach+1)] <- ApprList
#     
#     return(wide)
#   }
# )
# 
# # # length of 90% CI
# # CI_batch <- lapply(weighted_batch,function(x){
# #   apply(x,MARGIN = 2,function(x){
# #     (quantile(x,0.95)-quantile(x,0.05))/median(x)
# #   })
# # })
# # summary_CI <- lapply(
# #   CI_batch,function(x){
# #     temp <- as.data.frame(x)
# #     Nobs <- nrow(temp) / Napproach        # number of data points
# #     Xnames <- rownames(temp)[1:Nobs]      # X values
# #     # reorganize into wide format
# #     temp <- temp %>% mutate(
# #       X = rep(Xnames,times = Napproach),
# #       Appr = rep(1:Napproach,each = Nobs)
# #     ) %>% rename(CI = x)
# #     wide <- dcast(data = temp,
# #                   X ~ Appr, value.var = "CI")
# #     colnames(wide)[2:(Napproach+1)] <- ApprList
# #     
# #     return(wide)
# #   }
# # )
# # 
# # Demonstration of WAIC---------
# 
# 
# # temp <- summary_CI[[4]]
# # p_t <- vector("numeric",length = Nobs)
# # for(n in 1:Nobs){
# #   temp_vec <- unlist(temp[n,-1])
# #   temp_vec <- temp_vec- temp_vec[3]
# #   # temp_vec <- temp_vec[-3]
# #   ptest <- t.test(temp_vec)
# #   p_t[n] <- ptest$p.value
# # }
# 
# 
# Rel_median_batch <- vector("list",iterTotal)
# for(i in 1:iterTotal){
#   Rel_median_temp <- matrix(NA,nrow = Nobs,ncol = Napproach)
#   for(m in 1:Napproach){
#     Rel_median_temp[,m] <- summary_median[[i]][,1+m] / Yraw[[i]]
#   }
#   Rel_median_median <- apply(Rel_median_temp,2,median)
#   Rel_median_batch[[i]] <- Rel_median_median
# }
# Rel_median_mtx <- matrix(unlist(Rel_median_batch),ncol = Napproach,byrow = T)
# 
# Rel_Q5_batch <- vector("list",iterTotal)
# for(i in 1:iterTotal){
#   Rel_Q5_temp <- matrix(NA,nrow = Nobs,ncol = Napproach)
#   for(m in 1:Napproach){
#     Rel_Q5_temp[,m] <- summary_Q5[[i]][,1+m] / Yraw[[i]]
#   }
#   Rel_Q5_median <- apply(Rel_Q5_temp,2,median)
#   Rel_Q5_batch[[i]] <- Rel_Q5_median
# }
# 
# Rel_Q95_batch <- vector("list",iterTotal)
# for(i in 1:iterTotal){
#   Rel_Q95_temp <- matrix(NA,nrow = Nobs,ncol = Napproach)
#   for(m in 1:Napproach){
#     Rel_Q95_temp[,m] <- summary_Q95[[i]][,1+m] / Yraw[[i]]
#   }
#   Rel_Q95_median <- apply(Rel_Q95_temp,2,median)
#   Rel_Q95_batch[[i]] <- Rel_Q95_median
# }


# Rel_Dev_summary <- data.frame(
#   median = unlist(Rel_median_batch),
#   Q5 = unlist(Rel_Q5_batch),
#   Q95 = unlist(Rel_Q95_batch),
#   Approach = rep(ApprList,times = iterTotal),
#   Iter = rep(1:iterTotal,each = Napproach)
# )

# Wide Tibble------
get_Rel_Dev_wide <- function(df_weighted_wide,raw){
  df_weighted_wide <- df_weighted_wide %>% group_by(iter) %>% mutate(
    Yraw = rep(raw,Napproach),
    Odds = Wtd / Yraw
  )
  Rel_Dev <- df_weighted_wide$Odds
  
  return(Rel_Dev)
}

Rel_Dev_wide <- mapply(get_Rel_Dev_wide,weighted_batch,Yraw)
Rel_Dev_vec <- as.vector(Rel_Dev_wide)

weighted_wide <- do.call(rbind,weighted_batch) %>% mutate(
  Rel_Dev = Rel_Dev_vec,
  batch = rep(1:iterTotal,each = Nobs*Napproach*iterEff)
  ) %>% rename(
  Approach = Appr
)

weighted_wide_median <- weighted_wide %>% group_by(Approach,iter,batch) %>% 
  summarize(Median = median(Wtd)) %>% arrange(batch,iter,Approach)


# Plotting-------
ApprOrder <- c("WAIC","PBMABB","PBMA","LPML","Stacking","AIC","BIC_BMA")

df_plot <- weighted_wide_median
  
plot <- ggplot(data = df_plot,
               # aes(x=factor(Approach,level = ApprOrder)))+
               aes(x = Approach))+
  geom_violin(aes(y = Median,group = batch))+
  # geom_errorbar(aes(ymin = Q5,ymax=Q95,group = Approach),col="grey")+
  # geom_point(aes(y = median))+
  scale_y_continuous(name = "Relative Deviance from True Value")+
  scale_x_discrete(name = "Averaging Approaches")+
  theme_classic(base_size = 25)

svg(filename = "Averaging_Robust.svg",
    width = 12,height = 12)
print(plot)
dev.off()
