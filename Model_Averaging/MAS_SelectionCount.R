# Settings-------
# Settings-------
mywd <- "C:/Users/bks01/Downloads/working/"
fileDir <- "C:/Users/bks01/Downloads/working/"
# seed <- 47405      # use IUB postal code as seed for random simulation
DSName <- "Linear"
batch <- 100
iter <- 5

samplingiter <- 1E4*5   # number of iterations in stan sampling
warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 4       # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling
setwd(mywd)


source(file = paste0(fileDir,"BHBMD_model_utilities.R"),
       local = T)
batchsize <- iter
iterTotal <- batch * batchsize

# Load batch files------

weights_batch <- vector("list",iterTotal)
Metrics_batch <- vector("list",iterTotal)
Choice_batch <- vector("list",iterTotal)

for(b in 1:batch){
  for(i in 1:iter){
    index <- batchsize * (b-1)+i
    weights_temp <- read.csv(file = paste0("Weights_",DSName,"_batch_",b,"_",i,".csv"))
    ModelList <- weights_temp[,1]
    MethodList <- colnames(weights_temp)
    weights_temp <- weights_temp %>% select(-1)
    rownames(weights_temp) <- ModelList
    weights_batch[[index]] <- weights_temp
  }
}

for(b in 1:batch){
  for(i in 1:iter){
    index <- batchsize * (b-1)+i
    metrics_temp <- read.csv(file = paste0("Metrics_",DSName,"_batch_",b,"_",i,".csv"))
    MetricList <- metrics_temp[,1]
    metrics_temp <- metrics_temp %>% select(-1)
    rownames(metrics_temp) <- MetricList
    Metrics_batch[[index]] <- metrics_temp
  }
}

for(b in 1:batch){
  for(i in 1:iter){
    index <- batchsize * (b-1)+i
    choice_temp <- read.csv(file = paste0("Choice_",DSName,"_batch_",b,"_",i,".csv"))
    MetricList <- choice_temp[,1]
    choice_temp <- choice_temp %>% select(-1)
    rownames(choice_temp) <- MetricList
    Choice_batch[[index]] <- choice_temp
  }
}

# Get iterative Summary----
weights_mean <- Reduce("+",weights_batch) / iterTotal
Metrics_mean <- Reduce("+",Metrics_batch) / iterTotal
Choice_count <- Reduce("+", Choice_batch)

write.csv(weights_mean,
          file = paste0("Mean_weights_on_",DSName,"_data.csv"))
write.csv(Metrics_mean,
          file = paste0("Mean_metrics_on_",DSName,"_data.csv"))
write.csv(Choice_count,
          file = paste0("Choice_count_on",DSName,"_data.csv"))
