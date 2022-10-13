# Settings-------
# Settings-------
mywd <- "/N/scratch/zhezhou/MAS_Expo5/"
fileDir <- "/N/slate/zhezhou/"

# mywd <- paste0(getwd(),"/")
# fileDir <- mywd
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
filenames <- as.list(dir(pattern = "weighted_wide_summary_batch_*"))
lapply(filenames,load,.GlobalEnv)

# put all into a list
weighted_summary_list <- vector("list",batch)
for(i in 1:batch){
  weighted_summary_list[[i]] <- eval(parse(text = paste0("weighted_wide_summary_batch_",i)))
}

weighted_summary <- do.call(rbind,weighted_summary_list)

# weighted_summary_wide <- reshape2::melt(
#   weighted_summary,
#   id.vars = c("Approach","batch")
# ) %>% rename(Stat = variable, Dev = value)

# Plotting-------
ApprOrder <- c("WAIC","PBMABB","PBMA","LPML","Stacking","AIC","BIC_BMA")

df_plot <- weighted_summary
  
plot <- ggplot(data = df_plot
               # aes(x=factor(Approach,level = ApprOrder)))+
               # aes(x = Approach)
               )+
  facet_wrap( ~ Approach,ncol = 3)+
  geom_boxplot(aes(x = MedianDev,y=batch),col = "grey")+
  scale_x_continuous(
    name = "Relative Deviance from the True Value",
    breaks = c(1.0)
    )+
  scale_y_discrete(
    name = "Iterations",breaks = c(0,500)
  )+
  theme_classic(base_size = 25)
  # geom_point(aes(x = Stat,y = Dev,group = batch))
  # geom_violin(aes(y = median))+
  # geom_violin(aes(y=Q5))+
  # geom_violin(aes(y=Q95))+
  # geom_errorbar(aes(ymin = Q5,ymax=Q95),col="grey")+
  # geom_point(aes(y = median))+
  # scale_y_continuous(name = "Relative Deviance from True Value")+
  # scale_x_discrete(name = "Averaging Approaches")

png( "Averaging_Robust.png",width = 1600,height = 1600)
print(plot)
dev.off()
