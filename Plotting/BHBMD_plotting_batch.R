## Plotting iterative BHBMD results by batch
## Author: Zheng Zhou
## Date: Oct 24 2022

## Input-------
# mywd <- paste0(getwd(),"/")
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

task_ID <- Sys.getenv("SLURM_ARRAY_TASK_ID")
# task_ID <- 1

batch_ID <- task_ID
## Load data------
# raw data
df_raw <- read.csv(file = paste0(fileDir,"Metadata_bladder_iAs_censored.csv"))

# DataList with simulated ADD
DataList_batch <- vector("list",batchsize)
for(i in 1:batchsize){
  Data_temp <- read.csv(file = paste0(mywd, "DataList_",DSName,"_batch_",batch_ID,"_",i,".csv"))
  DataList_batch[[i]] <- Data_temp %>% mutate(
    RR.l = df_raw$RR.l,
    RR.u = df_raw$RR.u
  )
}

# stanfit posteriors
load(file = paste0(mywd, "fitList_",DSName,"_batch_",batch_ID,".rdata"))

## Plotting by DR model------

for(b in 1:batchsize){
  # assign(                ## create a fitList object with batch index
  #   paste0("fitList_batch_",batch_ID,"_",b),
  #   eval(parse(text = paste0("fitList_batch_",batch_ID)))[[b]]
  # )
  fitList <- eval(parse(text = paste0("fitList_batch_",batch_ID)))[[b]]  ## without batch index
  DataList <- DataList_batch[[b]]
  
  ## extract stanfit by dose response model
  ## Model list: Linear, Power, Hill, Expo5
  Fit_Linear <- fitList$Linear
  Fit_Power <- fitList$Power
  # Fit_MM <- fitList$MM
  Fit_Hill <- fitList$Hill
  Fit_Expo5 <- fitList$Expo5
  
  ## plot overarching and study-specific curves by DR model
  
  batchtag <- paste0("simu_ADD_batch_",batch_ID,"_",b)
  ## Linear
  plot_overarching_Linear(stanfit = Fit_Linear,DataList = DataList,
                          savetag = batchtag)
  plot_specific_Linear(Fit_Linear,DataList,savetag = batchtag)
  
  ## Power
  plot_overarching_Power(stanfit = Fit_Power,DataList = DataList,
                         savetag = batchtag)
  plot_specific_Power(Fit_Power,DataList,savetag = batchtag)
  
  # ## MM
  # plot_overarching_MM(stanfit = Fit_MM,DataList = DataList,
  #                     savetag = batchtag)
  # plot_specific_MM(Fit_MM,DataList,savetag = batchtag)

  ## Hill
  plot_overarching_Hill(stanfit = Fit_Hill,DataList = DataList,
                        savetag = batchtag)
  plot_specific_Hill(Fit_Hill,DataList,savetag = batchtag)
  
  ## Expo5
  plot_overarching_Expo5(stanfit = Fit_Expo5,DataList = DataList,
                         savetag = batchtag)
  plot_specific_Expo5(Fit_Expo5,DataList,savetag = batchtag)
  
  ## clear memory
  rm(list = c("fitList","DataList","batchtag",
                 "Fit_Linear","Fit_Power",
                 # "Fit_MM",
                 "Fit_Hill","Fit_Expo5"))
}

