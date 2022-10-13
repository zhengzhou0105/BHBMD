# Settings-------
mywd <- "/N/scratch/zhezhou/MAS_Hill/"
fileDir <- "/N/slate/zhezhou/"
iter <- 5
# seed <- 47405      # use IUB postal code as seed for random simulation
samplingiter <- 1E4*5   # number of iterations in stan sampling
DSName <- "Hill"


warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 3       # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling
setwd(mywd)


source(file = paste0(fileDir,"BHBMD_model_utilities.R"),
       local = T)
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# task_id <- 1
batchID <- task_id
batchsize <- iter

# Simulate Data------
data_iter <- vector("list",length = iter)
for(i in 1:iter){
  data_iter[[i]] <- fn_metaDRsimulation(S=10,G=3,
                                        ybound=c(0.1,20),
                                        sigmavalue =0.3)
}

# Fitting and Evaluation-----------
DataList_batch <- lapply(data_iter,function(x) {
  PrepareLnormData(fn_getsimuData(x,DSName),RR = F)
})
# Export datalist as csv files per iteration
mapply(write.csv,x = DataList_batch,
       file = paste0("DataList_",DSName,"_batch_",batchID,"_",1:batchsize, ".csv"))

# stan fitting
fitList_batch <- lapply(DataList_batch,fit_allmodels,DSName = DSName)
assign(
  paste0("fitList_batch_",batchID),
  fitList_batch
)
# export stanfit objects per iteration
save(list = paste0("fitList_batch_",batchID),file = paste0("fitList_",DSName,"_batch_",batchID,".rdata"))

weights_batch <- lapply(fitList_batch, getwts)
assign(
  paste0("weights_batch_",batchID),
  weights_batch
)
# Export weightings
save(list = paste0("weights_batch_",batchID),
     file = paste0("weights_",DSName,"_batch_",batchID,".rdata"))
mapply(write.csv,x = weights_batch,
       file = paste0("Weights_",DSName,"_batch_",batchID,"_",1:batchsize,".csv"))

# Performance evaluation
PerfMetrics_batch <- t(mapply(ModelEvalbyFd,      # a batchsize * 4 matrix
                              fitList = fitList_batch,
                              modelweights = weights_batch,
                              DataList = DataList_batch))

# Weighted_batch <- vector("list",batchsize)
# Weighted_mean_batch <- vector("list",batchsize)
# Metrics_batch  <- vector("list",batchsize)
# ChoiceList_batch <- vector("list",batchsize)
# 
# for(n in 1:batchsize){
#   Weighted_batch[[n]] <- PerfMetrics_batch[[n]]
#   Weighted_mean_batch[[n]] <- PerfMetrics_batch[[(n+batchsize)]]
#   Metrics_batch[[n]] <- PerfMetrics_batch[[(n+batchsize*2)]]
#   ChoiceList_batch[[n]] <- PerfMetrics_batch[[(n+batchsize*3)]]
# }
Weighted_batch <- PerfMetrics_batch[,1]
mapply(write.csv,x = Weighted_batch, 
       file = paste0("Weighted_batch_",DSName,"_batch_",batchID,"_",1:batchsize,".csv"))

Weighted_mean_batch <- PerfMetrics_batch[,2]
mapply(write.csv,x = Weighted_mean_batch,
       file = paste0("Weighted_mean_",DSName,"_batch_",batchID,"_",1:batchsize,".csv"),row.names =F)

Metrics_batch <- PerfMetrics_batch[,3]
mapply(write.csv,x = Metrics_batch,
       file = paste0("Metrics_",DSName,"_batch_",batchID,"_",1:batchsize,".csv"),row.names =T)

ChoiceList_batch <- PerfMetrics_batch[,4]
# SelectionCount_batch <- do.call("+",ChoiceList_batch)
mapply(
  write.csv,
  x = ChoiceList_batch,
  file = paste0("Choice_",DSName,"_batch_",batchID,"_",1:batchsize, ".csv"),
  row.names =T
)

# Generate figures
options(bitmapType='cairo')
mapply(ShowWeightedValues,
       DataList = DataList_batch,
       modelwts = weights_batch,
       weighted = Weighted_batch,
       dataName = paste0(DSName,"_batch_",batchID,"_",1:batchsize)
)


