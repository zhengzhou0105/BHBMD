## Estimate bladder cancer risks attributable to iAs exposure via rice intake among Chinese urban population
## Author: Zheng Zhou
## Date: Oct 25 2022

## Input-------
# mywd <- paste0(getwd(),"/")
# mywd <- "C:/Users/bks01/Downloads/Working/"
# fileDir <- mywd

mywd <- "/N/scratch/zhezhou/test/"
fileDir <- "/N/slate/zhezhou/"

DSName <- "simu_ADD"
batchsize <- 5
## Settings---------
# load utility functions
source(file = paste0(fileDir,"BHBMD_model_utilities.R"),
       local = T)
options(bitmapType='cairo')

iter <- 1e3                                  # iterations for  ADD simulation

task_ID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# task_ID <- 1

batch_ID <- task_ID

## Simulate target ADD----

## bioaccessibility from Zhou et al. 2021
## Exposure from Zhou et al. 2020
BC_iAs <- rbeta(iter,4.91,1.85)
Conc_iAs <- rlnorm(iter,-2.87,0.94)
IRBW_rice <- rlnorm(iter,0.14,1.42)
ADD_iAs_rice <- BC_iAs * Conc_iAs * IRBW_rice
# ADD_batch <- rep(list(ADD_iAs_rice),batchsize)
  
## Load data---

## model weights
Weights_batch <- vector("list",batchsize)
for(b in 1:batchsize){
  df_temp <- read.csv(
    paste0(mywd,"Weights_",DSName,"_batch_",batch_ID,"_",b,".csv")
  )
  df_temp <- df_temp %>% select(-"X")
  Weights_batch[[b]] <- df_temp
}

## model posterior parameters
load(
  paste0(mywd,"fitList_",DSName,"_batch_",batch_ID,".rdata")
)
fitList_batch <- eval(parse(text = paste0("fitList_batch_",batch_ID)))

## Estimate weighted relative risks------

# based on overarching parameters
for(b in 1:batchsize){
  fd_weighted <- get_Weighted(
    fitList =  fitList_batch[[b]],
    df_wts =  Weights_batch[[b]],input =  ADD_iAs_rice,
    level = "overarching",Index = NULL
  )
  save(list = "fd_weighted",
       file = paste0(mywd,"RR_weighted_overarching_batch_",batch_ID,"_",b,".rdata"))
  rm(fd_weighted)
}

# based on study 3
for(b in 1:batchsize){
  fd_weighted <- get_Weighted(
    fitList =  fitList_batch[[b]],
    df_wts =  Weights_batch[[b]],input =  ADD_iAs_rice,
    level = "specific",Index = 3
  )
  save(list = "fd_weighted",
       file = paste0(mywd,"RR_weighted_study3_batch_",batch_ID,"_",b,".rdata"))
  rm(fd_weighted)
}