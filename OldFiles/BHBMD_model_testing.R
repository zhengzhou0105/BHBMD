# Settings----

# global values
# seed <- 47405      # use IUB postal code as seed for random simulation
samplingiter <- 1E3   # number of iterations in stan sampling
warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 4    # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling

# load utility functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
fileDir <- paste0(getwd(),"/")
source("BHBMD_model_utilities.R",local = T)

# load data--------
# load("Hill meta base data.Rdata")
load("BHBMD base data.Rdata")
# load("Simulated Data for Model Eval.rdata")

DataList <- PrepareLnormData(df_bladder_meta,RR=T)
dataName <- "bladder_meta"

# Linear------

modelname <- "Lognormal_Summary_Linear_meta_3"

# stan sampling and save samples in a Rdata file
assign(paste0("stanfit_",modelname,"_",dataName),
       get_metaStanfit(modelname = modelname,dataName = dataName,Data = DataList))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))

# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(eval(parse(text = paste0("stanfit_",modelname,"_",dataName))))
df_posteriors <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, mu_b, sigma_b)
)
summary(df_posteriors[,c("a","b")])

# display study specific curves
showcurve_specific_Linear(df_posteriors = df_posteriors,
                          df_input = df_bladder_meta,
                          aname = "astar", bname = "bstar")

# overarching median curve
showcurve_overarching_Linear(df_posteriors = df_posteriors,
                             df_input = df_bladder_meta,
                             aname = "a" , bname = "b",
                             xup = 10)
# Power--------
modelname <- "Lognormal_Summary_Power_meta_4"

# stan sampling and save samples in a Rdata file
assign(paste0("stanfit_",modelname,"_",dataName),
       get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))

# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(eval(parse(text = paste0("stanfit_",modelname,"_",dataName))))

# MM--------
modelname <- "Lognormal_Summary_MM_meta_4"

# stan sampling and save samples in a Rdata file
assign(paste0("stanfit_",modelname,"_",dataName),
       get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))

# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(eval(parse(text = paste0("stanfit_",modelname,"_",dataName))))


# Hill----
modelname <- "Lognormal_Summary_Hill_meta_5"

# stan sampling and save samples in a Rdata file
assign(paste0("stanfit_",modelname,"_",dataName),
       get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))


# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(eval(parse(text = paste0("stanfit_",modelname,"_",dataName))))
df_PostPars <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b),
  c = rlnorm(iterEff, meanlog = mu_c, sdlog = sigma_c),
  g = rlnorm(iterEff, meanlog = mu_g, sdlog = sigma_g)
) 
summary(df_PostPars[,c("a","b","c","g")])

# display study specific curves
showcurve_specific_Hill(df_posteriors = df_PostPars,df_input = df_bladder_meta,
                        aname = "a" , bname = "b" , cname = "c" ,gname = "g")

# overarching median curve
showcurve_overarching_Hill(df_posteriors = df_PostPars,df_input = df_bladder_meta,
                           aname = "a" , bname = "b" , cname = "c" ,gname = "g",
                           xup = 10, yup = NULL)


# Exponential 5-----
modelname <- "Lognormal_Summary_Expo5_meta_5"

# stan sampling and save samples in a Rdata file
assign(paste0("stanfit_",modelname,"_",dataName),
       get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))

# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(eval(parse(text = paste0("stanfit_",modelname,"_",dataName))))

# Model Averaging--------
# Compare Linear 2 and Hill 2
ModComp_1 <- getwts(list(
  stanfit_Lognormal_Summary_Linear_meta_2_meta,
  stanfit_Lognormal_Summary_Hill_meta_2_meta
));ModComp_1

# Examine hierarchical structure in Linear equations
# No hier, only on intercept, on both intercept and slope
ModComp_2 <- getwts(list(   
  stanfit_Lognormal_Summary_Linear_meta_0_meta,
  stanfit_Lognormal_Summary_Linear_meta_1_meta,
  stanfit_Lognormal_Summary_Linear_meta_2_meta
));ModComp_2
  
# Examine hierarchical structure in Hill equation
# No hier, and stepwise in a b c g

ModComp_3 <- getwts(list(
  stanfit_Lognormal_Summary_Hill_meta_0_meta,
  stanfit_Lognormal_Summary_Hill_meta_1_meta,
  stanfit_Lognormal_Summary_Hill_meta_2_meta,
  stanfit_Lognormal_Summary_Hill_meta_3_meta,
  stanfit_Lognormal_Summary_Hill_meta_4_meta
))

# Compare Linear 2 and Hill 4
ModComp_4 <- getwts(list(
  stanfit_Lognormal_Summary_Linear_meta_2_meta,
  stanfit_Lognormal_Summary_Hill_meta_4_meta
));ModComp_4

# Simulated Linear data 
fitList_simu_Linear <- list(
  Linear = stanfit_Lognormal_Summary_Linear_meta_2_simu_Linear,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_simu_Linear,
  Power = stanfit_Lognormal_Summary_Power_meta_3_simu_Linear,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_simu_Linear
)
weights_simu_Linear <- getwts(fitList_simu_Linear);

PerfMetrics_simu_Linear <- ModelEvalbyFd(fitList = fitList_simu_Linear,
              modelweights = weights_simu_Linear,
              DataList = DataList);
PerfMetrics_simu_Linear
weights_simu_Linear

Results_simu_Linear <- list(
  RawData = DataList,
  fitList_simu_Linear = fitList_simu_Linear,
  weights_simu_Linear = weights_simu_Linear,
  ModelMetrics =  PerfMetrics_simu_Linear$ModelMetrics,
  ModelChoice = PerfMetrics_simu_Linear$ModelChoice,
  WeightedMean = PerfMetrics_simu_Linear$WeightedMean,
  WeightedQuantities = PerfMetrics_simu_Linear$WeightedQuantities
)
save(Results_simu_Linear,file = "Simulated Linear Model Eval Results.rdata")


# Simulated Power data 
fitList_simu_Power <- list(
  Linear = stanfit_Lognormal_Summary_Expo5_meta_4_simu_Power,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_simu_Power,
  Power = stanfit_Lognormal_Summary_Power_meta_3_simu_Power,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_simu_Power
)
weights_simu_Power <- getwts(fitList_simu_Power);

PerfMetrics_simu_Power <- ModelEvalbyFd(fitList = fitList_simu_Power,
                                    modelweights = weights_simu_Power,
                                    DataList = DataList);
PerfMetrics_simu_Power
weights_simu_Power

Results_simu_Power <- list(
  RawData = DataList,
  fitList_simu_Power = fitList_simu_Power,
  weights_simu_Power = weights_simu_Power,
  ModelMetrics =  PerfMetrics_simu_Power$ModelMetrics,
  ModelChoice = PerfMetrics_simu_Power$ModelChoice,
  WeightedMean = PerfMetrics_simu_Power$WeightedMean,
  WeightedQuantities = PerfMetrics_simu_Power$WeightedQuantities
)
save(Results_simu_Power,file = "Simulated Power Model Eval Results.rdata")

# Simulated MM data 
fitList_simu_MM <- list(
  Linear = stanfit_Lognormal_Summary_Expo5_meta_4_simu_MM,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_simu_MM,
  Power = stanfit_Lognormal_Summary_Power_meta_3_simu_MM,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_simu_MM
)
weights_simu_MM <- getwts(fitList_simu_MM);

PerfMetrics_simu_MM <- ModelEvalbyFd(fitList = fitList_simu_MM,
                                 modelweights = weights_simu_MM,
                                 DataList = DataList)
PerfMetrics_simu_MM
weights_simu_MM

Results_simu_MM <- list(
  RawData = DataList,
  fitList_simu_MM = fitList_simu_MM,
  weights_simu_MM = weights_simu_MM,
  ModelMetrics =  PerfMetrics_simu_MM$ModelMetrics,
  ModelChoice = PerfMetrics_simu_MM$ModelChoice,
  WeightedMean = PerfMetrics_simu_MM$WeightedMean,
  WeightedQuantities = PerfMetrics_simu_MM$WeightedQuantities
)
save(Results_simu_MM,file = "Simulated MM Model Eval Results.rdata")

# Simulated Hill data
fitList_simu_Hill <- list(
  Linear = stanfit_Lognormal_Summary_Linear_meta_2_simu_Hill,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_simu_Hill,
  Power = stanfit_Lognormal_Summary_Power_meta_3_simu_Hill,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_simu_Hill
)
weights_simu_Hill <- getwts(fitList_simu_Hill);

PerfMetrics_simu_Hill <- ModelEvalbyFd(fitList = fitList_simu_Hill,
                                 modelweights = weights_simu_Hill,
                                 DataList = DataList);
PerfMetrics_simu_Hill
weights_simu_Hill

Results_simu_Hill <- list(
  RawData = DataList,
  fitList_simu_Hill = fitList_simu_Hill,
  weights_simu_Hill = weights_simu_Hill,
  ModelMetrics =  PerfMetrics_simu_Hill$ModelMetrics,
  ModelChoice = PerfMetrics_simu_Hill$ModelChoice,
  WeightedMean = PerfMetrics_simu_Hill$WeightedMean,
  WeightedQuantities = PerfMetrics_simu_Hill$WeightedQuantities
)
save(Results_simu_Hill,file = "Simulated Hill Model Eval Results.rdata")


# Simulated Expo5 data
fitList_simu_Expo5 <- list(
  Linear = stanfit_Lognormal_Summary_Linear_meta_2_simu_Expo5,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_simu_Expo5,
  Power = stanfit_Lognormal_Summary_Power_meta_3_simu_Expo5,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_simu_Expo5
)
weights_simu_Expo5 <- getwts(fitList_simu_Expo5);

PerfMetrics_simu_Expo5 <- ModelEvalbyFd(fitList = fitList_simu_Expo5,
                                  modelweights = weights_simu_Expo5,
                                  DataList = DataList);
PerfMetrics_simu_Expo5
weights_simu_Expo5

Results_simu_Expo5 <- list(
  RawData = DataList,
  fitList_simu_Expo5 = fitList_simu_Expo5,
  weights_simu_Expo5 = weights_simu_Expo5,
  ModelMetrics =  PerfMetrics_simu_Expo5$ModelMetrics,
  ModelChoice = PerfMetrics_simu_Expo5$ModelChoice,
  WeightedMean = PerfMetrics_simu_Expo5$WeightedMean,
  WeightedQuantities = PerfMetrics_simu_Expo5$WeightedQuantities
)
save(Results_simu_Expo5,file = "Simulated Expo5 Model Eval Results.rdata")
ShowWeightedValues(DataList = Results_simu_Expo5$RawData,
                   modelwts = Results_simu_Expo5$weights_simu_Expo5,
                   weighted = Results_simu_Expo5$WeightedQuantities,
                   dataName = "simu_Expo5")

# Raw bladder data
fitList_bladder <- list(
  Linear = stanfit_Lognormal_Summary_Linear_meta_2_bladder,
  Hill = stanfit_Lognormal_Summary_Hill_meta_4_bladder,
  Power = stanfit_Lognormal_Summary_Power_meta_3_bladder,
  Expo5 = stanfit_Lognormal_Summary_Expo5_meta_4_bladder
)
weights_bladder <- getwts(fitList_bladder);

PerfMetrics_bladder <- ModelEvalbyFd(fitList = fitList_bladder,
                                   modelweights = weights_bladder,
                                   DataList = DataList);
PerfMetrics_bladder
weights_bladder

Results_bladder <- list(
  RawData = DataList,
  fitList_bladder = fitList_bladder,
  weights_bladder = weights_bladder,
  ModelMetrics =  PerfMetrics_bladder$ModelMetrics,
  ModelChoice = PerfMetrics_bladder$ModelChoice,
  WeightedMean = PerfMetrics_bladder$WeightedMean,
  WeightedQuantities = PerfMetrics_bladder$WeightedQuantities
)
save(Results_bladder,file = "Raw bladder data Model Eval Results.rdata")

CalculateMetrics(baseline = Results_simu_Expo5$RawData$ymean, 
                 weighted = rstan::extract(Results_simu_Expo5$fitList_simu_Expo5$Hill)$fd)

ShowWeightedValues(DataList = Results_bladder$RawData,
                   modelwts = Results_bladder$weights_bladder,
                   weighted = Results_bladder$WeightedQuantities,
                   dataName = "bladder")

