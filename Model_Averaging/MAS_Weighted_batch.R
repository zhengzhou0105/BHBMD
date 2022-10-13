# Settings-------
# Settings-------
mywd <- "/N/scratch/zhezhou/MAS_Expo5/"
fileDir <- "/N/slate/zhezhou/"

# mywd <- paste0(getwd(),"/")
# fileDir <- mywd
# seed <- 47405      # use IUB postal code as seed for random simulation
DSName <- "Expo5"
batchsize <- 5

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# task_id <- 1
batchID <- task_id
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
batch <- 100
iter <- batchsize
iterTotal <- batch * batchsize

iterName <- paste0(batchID,"_",1:batchsize)

# Load batch files------

weighted_batch <- vector("list",batchsize)
Yraw <- vector("list",batchsize)

# load data
for(i in 1:batchsize){
  weighted_temp <- read.csv(file = paste0("Weighted_batch_",DSName,"_batch_",batchID,"_",i,".csv"))
  weighted_batch[[i]] <- weighted_temp
  rm(weighted_temp)
}
for(i in 1:batchsize){
  Data_temp <- read.csv(file = paste0("DataList_",DSName,"_batch_",batchID,"_",i,".csv"))
  Yraw_temp <- Data_temp$ymean
  Yraw[[i]] <- Yraw_temp
  rm(Yraw_temp)
  rm(Data_temp)
}

# transform to long format
weighted_batch <- lapply(weighted_batch,function(x) {
  weighted_temp <- x
  weighted_temp <- weighted_temp %>% select(-1)
  Nobs <- ncol(weighted_temp) / Napproach
  Xnames <- colnames(weighted_temp)[1:Nobs]
  Xvalue <- as.numeric(stringr::str_replace(Xnames,pattern = "X",""))
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
  return(weighted_wide)
  rm(weighted_temp)
  rm(weighted_wide)
  rm(weighted_wide_temp)
})

# Calculate relative deviance------
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

Nobs <- length(Yraw[[1]])
  
weighted_wide <- do.call(rbind,weighted_batch) 
weighted_wide <- weighted_wide %>% mutate(
  Rel_Dev = Rel_Dev_vec,
  batch = paste0(batchID,"_",rep(1:batchsize,each = Nobs*Napproach*iterEff))
) %>% rename(
  Approach = Appr
)

weighted_wide_summary <- weighted_wide %>% group_by(Approach,batch,iter) %>% 
  summarize(MedianDev = median(Rel_Dev)) %>% arrange(batch,Approach,iter)

assign(x = paste0("weighted_wide_summary_batch_",batchID),
weighted_wide_summary)

save(
  list = paste0("weighted_wide_summary_batch_",batchID),
  file = paste0("weighted_wide_summary_batch_",batchID,".rdata")
)
