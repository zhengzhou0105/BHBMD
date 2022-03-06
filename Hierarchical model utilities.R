# ----Settings----
#* Environment----
set.seed(seed = seed)

library("rstan")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options(scipen = 10)
library(dplyr)
library(ggplot2)
library(xlsx)

iterEff <- (iter*(1-warmup)*Nchain)

# ----Utility functions------

#* general-------
dnormalize_minmax <- function(d){     # normalize dose to 0 and 1
  dmin <- min(d)          # applicable to single study and meta regression
  dmax <- max(d)
  dm <- (d-dmin)/(dmax-dmin)
  return(dm)
}

dnormalize_meanSD <- function(d){      # normalize to standardized Z
  dmean <- mean(d)
  dsd <- sd(d)
  dm <- (d-dmean)/dsd
  return(dm)
}

get_ysdL <- function(ymean,ysd){      # estimated summarized ysd at log scale
  ysdL <- sqrt(log((ysd/ymean)^2+1))   # based on summarized ymean and ysd 
  return(ysdL)                         # at regular scale
}

get_ymeanL <- function(ymean,ysd){   # estimate summarized ymean at log scale
  ysdL <- sqrt(log((ysd/ymean)^2+1)) # based on summarized ymean and ysd
  ymeanL <- log(ymean) - 0.5*ysdL^2    # at regular scale
  return(ymeanL)
}


get_ysd <- function(ymean,ysdL){         # anti of get_ysdL
  ysd <- ymean * sqrt(exp(ysdL^2)-1)    # applicable to RR data
  return(ysd)                            # get ysd regular from ysd log
}

get_Hill <- function(x,a,b,c,g){     
  # calculate response given Hill parameters
  # can accommodate raw and normalized dose
  # but normalized dose recommended
  # parameters must be real numbers
  if(length(a) != 1) stop("length of a not 1") 
  if(length(b) != 1) stop("length of b not 1")
  if(length(c) != 1) stop("length of c not 1")                    
  if(length(g) != 1) stop("length of g not 1")                     
  x <- unlist(x)                # x can be a vector
  
  # returned vector has the same length of x
  y <- numeric(length(x))
  y <- a + (b * Re(as.complex(x)^g))/(c^g + Re(as.complex(x)^g))
  return(y)
}

get_Linear <- function(x,a,b){     
  # calculate response given Linear parameters
  # can accomendate raw and normalized dose
  # but normalized dose recommended
  # parameters must be vector of length 1
  if(length(a) != 1) stop("length of a not 1") 
  if(length(b) != 1) stop("length of b not 1")
  x <- unlist(x)                # x can be a vector
  
  # returned vector has the same length of x
  y <- numeric(length(x))
  y <- a + b * x
  return(y)
}

# *Data Treatment-------

# Standardize data into the df_input format
RRDataUnify <- function(df_raw,
                    Index,N,Study,
                    dose = NULL,dosed=NULL,
                    RR=NULL,RRlog=NULL,
                    RRlower = NULL,RRupper = NULL,
                    ysd=NULL,ysdL=NULL){
  # For single study data
  if(is.null(Index) & is.null(Study)){
    ndose <- nrow(df_raw)
    vec_Index <- rep(1,ndose)
    vec_Study <- rep("default",ndose)
  } else {
    vec_Index <- df_raw[,Index]
    vec_Study <- df_raw[,Study]
    
  }

  # Number of subjects
  vec_N <- df_raw[,N]
  
  # get mean log of response
  if(is.null(RR) & !is.null(RRlog)){
    vec_RRlog <- df_raw[,RRlog]
    vec_RR <- exp(vec_RRlog)
  } 
  if(is.null(RRlog) & !is.null(RR)){
    vec_RR <- df_raw[,RR]
    vec_RRlog <- log(vec_RR)
  } 
  # get SD log of response
  if(is.null(ysd) & !is.null(ysdL)){
    vec_sdL <- df_raw[,ysdL]
    vec_sd <- get_ysd(ymean = vec_RR,ysdL = vec_sdL)
  }
  if(is.null(ysdL) & !is.null(ysd)){
    vec_sd <- df_raw[,ysd]
    vec_sdL <- get_ysdL(ymean = vec_RR,ysd = vec_sd)
  } 
  # get standardized dose
  if(is.null(dose) & !is.null(dosed)){
    vec_dosed <- df_raw[,dosed]
    stop("Only standardized dose available! Input original dose.")
  } 
  if(is.null(dosed) & !is.null(dose)){
    vec_dose <- df_raw[,dose]
    min_dose <- min(vec_dose)
    max_dose <- max(vec_dose)
    vec_dosed <- (vec_dose - min_dose) / (max_dose - min_dose)
  } 
  # Lower and upper bound of RR
  if(is.null(RRlower) & is.null(RRupper)){
    ndose <- nrow(df_raw)
    vec_RRlow <- rep(NA,ndose)
    vec_RRup <- rep(NA,ndose)
  } else {
    vec_RRlow <- df_raw[,RRlower]
    vec_RRup <- df_raw[,RRupper]
  }
  
  # Output
  df_unified <- data.frame(
    Index = unlist(vec_Index),
    Study = unlist(vec_Study),
    N = unlist(vec_N),
    RR = unlist(vec_RR),
    RR.l = unlist(vec_RRlow),
    RR.u = unlist(vec_RRup),
    RRlog = unlist(vec_RRlog),
    ysdL = unlist(vec_sdL),
    dose_minmax = unlist(vec_dosed),
    dose = unlist(vec_dose)
  )
  
  return(df_unified)
}


# reorder data and output for testing
ReorderData <- function(df_input,OrderNew){
  # extract data by group into a list, which is ordered by the new order
  # then concatenate into a dataframe
  
  # number of groups
  S <- length(OrderNew)
  ls_all <- list()        # list for storage
  IndexNew <- integer()
  # extract data
  for(s in 1:S){
    ls_all[[s]] <- df_input %>% filter(Study == OrderNew[s])
    # store the new index
    IndexNew <- c(IndexNew,rep(s, nrow(ls_all[[s]])))
  }
  # concatenate
  df_output <- do.call("rbind",ls_all)
  # Re-indexing
  df_output <- df_output %>% mutate(
    Index = IndexNew
  )

  return(df_output)
}

ExtractPosteriorbyChain <- function(df_posteriors,chains){
  IterStart <- (chains - 1)*iter * (1-warmup)+1
  IterEnd <- chains * iter * (1-warmup)
  df_IterIndex <- as_tibble(cbind(IterStart,IterEnd))
  df_IterIndex <- df_IterIndex %>% rowwise() %>% mutate(
    IterExtract = list(seq(from = IterStart, to = IterEnd))
  )
  IterExtract <- unlist(df_IterIndex$IterExtract)
  df_extract <- df_posteriors[IterExtract,]
}

ExtractRawData <- function(df_input){
  # df_temp <- df_input %>% group_by(Study) %>% mutate(
  #   GID = row_number()
  # ) %>% filter(GID != 1)
  df_temp <- df_input
  
  NMeta <- nrow(df_temp)               # number of rows in long format
  NStudy <- length(unique(df_temp$Index))   # number of study
  Index <- df_temp$Index                   # study index
  Nsub <- df_temp$N                   # number of subject at each dose level
  dose <- df_temp$dose                 # dose
  dosemeta <- df_temp$dose_minmax         # dose standarized in minmax
  ymean <- df_temp$RR                   # response mean at regular scale
  ymeanL <- df_temp$RRlog               # response mean at log scale
  # ysd <- df_temp$ysd                        # response SD at regular scale
  ysdL <- df_temp$ysdL                  # response SD at log scale
  # nstart <- integer(length = NStudy)
  # nend <- nstart
  # for(s in 1:NStudy){
  #   if(s == 1){
  #     nstart[s] <- 1
  #     nend[s] <- table(SStudy)[1]
  #   } else {
  #     nstart[s] <- sum(table(SStudy)[1:(s-1)])+1
  #     nend[s] <- sum(table(SStudy)[1:s])
  #   }
  # }
  
  RawData <- list(
    N = NMeta,
    S = NStudy,
    Index = Index,
    Nsub = Nsub,
    dose = dose,
    dosed = dosemeta,
    ymeanL = ymeanL,
    ysdL = ysdL
    # nstart = as.integer(nstart),
    # nend = as.integer(nend),
  )
  return(RawData)
}

DFtoMatrix <- function(df_input){  # convert long data frame to wide matrices 
  # add dose group index
  df_temp <- df_input %>% group_by(Study) %>%
    mutate(GID = row_number())
  # df_temp <- filter(df_temp,GID != 1)
  
  G <- as.integer(table(df_temp$Study))    # # of dose groups each study
  Gmax <- max(G)                 # largest number of dose group
  S <- max(df_temp$Index)                 # number of study
  # Ginits <- integer(length = S)
  # for(s in 1:S){
  #   Ginits[s] = ifelse(s == 1, 1, sum(G[1:(s-1)])+1)
  # }
  # matrices to store dose, mean SD and N
  # one study per row
  # one dose group per column
  mtx_doseD <- matrix(NA,nrow = S,ncol = Gmax) 
  mtx_dose <- mtx_doseD
  mtx_ymeanL <- mtx_doseD
  mtx_ysdL <- mtx_doseD
  mtx_n <- mtx_doseD
  
  for(s in 1:S){
    # extract data by study
    df_study <- df_temp %>% filter(Index == s)
    g_temp <- nrow(df_study)              # number of dose groups in the study
    # extract dose, ymean ysd and N
    vec_doseD_raw <- df_study$dose_minmax     # vector of standarized dose
    vec_dose_raw <- df_study$dose           # dose
    vec_ymeanL_Raw <- df_study$RRlog        # vector of mean at log
    vec_ysdL_raw <- df_study$ysdL            # vector of SD at log
    vec_n_raw <- df_study$N                  # vector of N
    
    # supplement with NA if # dose group less than max # dose group 
    if(g_temp < Gmax){
      NACount <- Gmax - g_temp
      vec_doseD <- c(vec_doseD_raw,rep(2,NACount))
      vec_dose <- c(vec_dose_raw,rep(2,NACount))
      vec_ymeanL <- c(vec_ymeanL_Raw,rep(2,NACount))
      vec_ysdL <- c(vec_ysdL_raw,rep(2,NACount))
      vec_n <- c(vec_n_raw,rep(2,NACount))
    } else {
      vec_doseD <- vec_doseD_raw
      vec_dose <- vec_dose_raw
      vec_ymeanL <- vec_ymeanL_Raw
      vec_ysdL <- vec_ysdL_raw
      vec_n <- vec_n_raw
    }
    
    
    # fill in completed data into matrix
    mtx_doseD[s,] <- vec_doseD
    mtx_dose[s,] <- vec_dose
    mtx_ymeanL[s,] <- vec_ymeanL
    mtx_ysdL[s,] <- vec_ysdL
    mtx_n[s,] <- vec_n
  }
  
  # Output in a list
  ls_output <- list(mtx_n,mtx_doseD,mtx_dose,mtx_ymeanL,mtx_ysdL,G,S,
                    # Ginits,
                    Gmax)
  names(ls_output) <- c("mtx_N","mtx_doseD","mtx_dose","mtx_ymeanL","mtx_ysdL","G","S",
                        # "Ginits",
                        "Gmax")
  return(ls_output)
}

#* Hill specific------------
get_prior_a <- function(ymean,ysd){       # get hyperparameters for intercept a
  ymax <- max(ymean)                      # applicable to normalized response
  idmax <- which(ymean == max(ymean))
  ysd_ymax <- ysd[idmax]
  a_upper <- 2 * (ymax + 2 * ysd_ymax)
  prior_a <- c(0,a_upper)
  return(prior_a)
}

get_prior_b <- function(ymean,ysd,dose){   # get hyperparameter for slope b
  idmax <- which(dose == max(dose))        # applicable to normalized response
  ymean_maxd <- ymean[idmax]
  ysd_maxd <- ysd[idmax]
  
  idmin <- which(dose == min(dose))
  ymean_mind <- ymean[idmin]
  ysd_mind <- ysd[idmin]
  
  b_slope <- (ymean_maxd+2*ysd_maxd-ymean_mind-2*ysd_mind)/
    (dose[idmax]-dose[idmin])
  prior_b <- c(0,b_slope * 5)
  return(prior_b)
}

#* meta Hill-----------

# get hyperparameter intercept a in meta data
get_prior_meta_a <- function(index,ymean,ysd){   
  # S is the number of studies in the meta data
  S <- length(unique(index))
  # estimate hyperparameter a for each study
  meta_a <- vector(mode = "numeric",length = S)
  for(i in 1:S){
    # id_temp is a vector of logical indicator of selected study
    id_temp <- which(index == i)   
    meta_a[i] <- get_prior_a(ymean[id_temp],ysd[id_temp])[2]
  }
  return(meta_a)
}

# get hyperparameter slope b in meta data
get_prior_meta_b <- function(index,ymean,ysd,dose){  
  # estimate b for each study
  # S is the number of study in the meta data
  S <- length(unique(index))
  meta_b <- vector(mode = "numeric",length = S)
  for(i in 1:S){
    id_temp <- which(index == i)
    meta_b[i] <- get_prior_b(ymean[id_temp],ysd[id_temp],dose[id_temp])[2]
  }
  return(meta_b)
}

# Model fitting---------
get_metaStanfit <- function(modelname,df_input,dataName,inits = NULL){
  
  # load models
  source("Hill model scripts.R")
  
  # Convert long format dataframe to wide format matrix
  ls_widedata <- DFtoMatrix(df_input = df_input)
  
  # Raw data
  RawData <- ExtractRawData(df_input)
  
  # Merge Raw Data with converted matrices
  dataAll <- append(
    RawData,
    list(
      G = ls_widedata[["G"]],
      Gmax = ls_widedata[["Gmax"]],
      mtx_doseD = ls_widedata[["mtx_doseD"]],
      mtx_dose = ls_widedata[["mtx_dose"]],
      mtx_ymeanL = ls_widedata[["mtx_ymeanL"]],
      mtx_ysdL = ls_widedata[["mtx_ysdL"]],
      mtx_N = ls_widedata[["mtx_N"]]
    )
  )
  
  
  # Priors
  
  # Linear - Vectorized function declaration 
  prior_Linear_meta_fun_vec_rag <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2)
  )
  
  # Hill- Vectorized function declaration
  prior_Hill_meta_fun_vec_rag <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )
  
  
  # Do not transform dose, fd not stand alone
  prior_Hill_meta_fun_index <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )

  # Do not transform dose, fd stand alone
  prior_Hill_meta_fun_index_fd <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )  
  # array vector indexing fd within study before likelihood
  prior_Hill_meta_Shao_fd_old_ver2 <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
    
  )
  
  # dataframe style double indexing no fd
  prior_Hill_meta_Shao_old <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )
  
  # fd within study using double loop before likelihood
  prior_Hill_meta_Shao_old_ver2 <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )
  
  # fd within study using double loop within likelihood
  prior_Hill_meta_Shao_old_ver3 <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )
  
  # # dataframe style double indexing fd as parameter
  # prior_Hill_meta_Shao_fd_old <- list(
  # prior_sigma = c(0,5),
  # prior_mu_a = c(-2,2),
  # prior_sigma_a = c(0,3),
  # prior_mu_b = c(-3,3),
  # prior_sigma_b = c(0,3),
  # prior_c = c(0,10),
  # prior_g = c(0,50)
  # )
  
  # dataframe style double indexing fd generated
  prior_Hill_meta_Shao_fdgen_old <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50)
  )
  
  prior_Hill_meta_Shao_nonrag <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    # dose = ls_widedata[["mtx_dose"]],
    # ymeanL = ls_widedata[["mtx_ymeanL"]],
    # ysdL = ls_widedata[["mtx_ysdL"]],
    mtx_Nsub = ls_widedata[["mtx_N"]],
    Gmax = ls_widedata[["Gmax"]],
    g_constr = c(0)
  )

  prior_Hill_meta_Shao_rag <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    # prior_mu_c = c(-2,2),
    # prior_sigma_c = c(0,2),
    prior_g = c(0,50),
    # dose = ls_widedata[["mtx_dose"]],
    # ymeanL = ls_widedata[["mtx_ymeanL"]],
    # ysdL = ls_widedata[["mtx_ysdL"]],
    # mtx_Nsub = ls_widedata[["mtx_N"]],
    Gmax = ls_widedata[["Gmax"]],
    g_constr = c(0)
  )  
  # # dataframe style single indexing double loop fd as parameter
  # prior_Hill_meta_Shao_fd_old1 <- list(
  # prior_sigma = c(0,5),
  # prior_mu_a = c(-2,2),
  # prior_sigma_a = c(0,3),
  # prior_mu_b = c(-3,3),
  # prior_sigma_b = c(0,3),
  # prior_c = c(0,10),
  # prior_g = c(0,18)
  # 
  # )
  # 
  # # dataframe style single indexing double loop fd generated
  # prior_Hill_meta_Shao_fdgen_old1 <- list(
  # prior_sigma = c(0,5),
  # prior_mu_a = c(-2,2),
  # prior_sigma_a = c(0,3),
  # prior_mu_b = c(-3,3),
  # prior_sigma_b = c(0,3),
  # prior_c = c(0,10),
  # prior_g = c(0,18)
  # 
  # )
  
  # matrix style double loop no fd
  prior_Hill_meta_Shao_matrix <- list(
    prior_sigma = c(0,5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_c = c(0,10),
    prior_g = c(0,50),
    g_constr = c(0)
  )
  
  # # non centered gUnif priors
  # prior_Hill_meta_NC_gUnif <- list(
  #   prior_sigma = c(0,2.5),
  #   prior_k = c(0,15),
  #   prior_p = c(1,15),
  #   
  #   prior_g = c(0,median(get_prior_meta_a(index = SStudy,ymean = ymean, 
  #                                         ysd = ysd))*3),
  #   prior_v = c(0,median(get_prior_meta_b(index = SStudy,ymean = ymean,
  #                                         ysd = ysd,dose = dosemeta)))
  # )
  # 
  # # non centered gLnorm priors
  # prior_Hill_meta_NC_gLnorm <- list(
  #   prior_sigma = c(0,2.5),
  #   prior_k = c(0,15),
  #   prior_p = c(1,15),
  #   
  #   prior_g = c(0,10),
  #   prior_v = c(0,median(get_prior_meta_b(index = SStudy,ymean = ymean,
  #                                         ysd = ysd,dose = dosemeta)))
  # )
  
  # centered g follows lnorm
  prior_Hill_meta_CT_gLnorm <- list(
    prior_sigma = c(0,2.5),
    prior_c = c(0,15),
    prior_g = c(1,15),
    prior_mu_a = c(0,10),
    prior_sigma_a = c(0,2.5),
    prior_s_b = c(0,15),
    prior_r_b = c(0,15)
  )
  
  # centered g follows lnorm2
  prior_Hill_meta_CT_gLnorm2 <- list(
    prior_sigma = c(0,2.5),
    prior_c = c(0,15),
    prior_g = c(1,15),
    prior_sigma_a = c(0,2.5),
    prior_s_b = c(0,15),
    prior_r_b = c(0,15)
  )
  
  # CT gUnif
  prior_Hill_meta_CT_gUnif <- list(
    prior_sigma = c(0,2.5),
    prior_c = c(0,15),
    prior_g = c(1,15),
    prior_a_a = c(0,10),
    # prior_b_a = c(0,max(get_prior_meta_a(index = SStudy,
    #                                      ymean = ymean,
    #                                      ysd = ysd))*3),
    prior_b_a = c(0,30),
    prior_s_b = c(0,15),
    prior_r_b = c(0,15)
  )
  
  # Linear
  prior_Linear_meta <- list(
    prior_sigma = c(0,2.5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2)
  )
  
  # Hill_meta_full_CT
  prior_Hill_meta_full_CT <- list(
    prior_sigma = c(0,2.5),
    prior_mu_a = c(-1,1),
    prior_sigma_a = c(0,2),
    prior_mu_b = c(-3,3),
    prior_sigma_b = c(0,2),
    prior_mu_c = c(-1,1),
    prior_sigma_c = c(0,2),
    prior_mu_g = c(-4,4),
    prior_sigma_g = c(0,2)
  )
  
  DataList <- append(dataAll,
                     eval(parse(text = paste0("prior_",
                                              modelname)))
                     )
  
  inits <- ifelse(is.null(inits),"random",inits)
  
  # Compile stanmodels
  assign(paste0("model_",modelname), 
         # stan_model(file = paste0("stanmodel_",modelname,".stan"))
         stan_model(model_code = eval(parse(text = paste0("modelstring_",modelname))))
  )

  stanfit <- stan(
                  # object =  eval(parse(text = paste0("model_",
                  #                                       modelname))),
                  # model_name = paste0("model_",
                  #                     modelname),
                  model_code = eval(parse(text = paste0("modelstring_",
                                                        modelname))),
                  data = DataList,
                  iter = iter,
                  chains = Nchain,
                  warmup = warmup * iter,
                  thin = thin,
                  seed = seed,
                  init = inits,
                  # open_progress = T,
                  init_r = 2,
                  control = list(
                    # adapt_delta  = 0.999,
                    # max_treedepth = 20
                  )
                  )
  SaveFileName <- paste0("stanfit_",modelname,"_",dataName,".RData")
  save(stanfit,file = SaveFileName)
  return(stanfit)
}

get_fitsummary <- function(stanfit){
  print(stanfit)
  summary <- as.data.frame(summary(stanfit)[1])
  # keep only stats in the name
  colnames(summary) <- substring(colnames(summary),9,20)
  # write.csv(summary,file = "posterior summary.csv")
  # modelname <- substr(deparse(substitute(stanfit)),9,100)
  df_posteriors <- as.data.frame(stanfit)
  return(df_posteriors)
}

# relative change in central tendency
getBMD_rel_Hill <- function(BMR,refdose,a,b,c,g){
  # if all parameters are of the same length
  iter <- length(a)
  parslengthNO <- (length(b) != iter) | (length(c) != iter) | (length(g) != iter)
  if(parslengthNO) stop("Parameters are not of the same length")
  
  yrefdose <- numeric(iter)
  yBMD <- numeric(iter)
  BMD <- numeric(iter)
  for(i in 1:iter){
    # calculate response at refernce dose level
    # notice the refdose is the necessarily standarized
    yrefdose[i] <- get_Hill(refdose,
                            a[i],b[i],c[i],g[i])
    # hybrid BMR
    yBMD[i] <- (BMR + 1) * yrefdose[i]
    # the length of initial value of x is the same as
    # the length of fitted parameters
    # jac = whether to return Jacobian = NULL
    # BMD_d is estimated by solving the Hill equation 
    # at given parameter values from posterior
    BMD[i] <- nleqslv::nleqslv(x = 0.33,
                                 get_Hill,jac = NULL,
                                 a[i],b[i],c[i],g[i])[[1]]
  }
  # BMD is at the same scale as the reference dose
  # BMD is the same length of parameters= iter
  return(BMD)

}

getBMD_rel_Linear <- function(BMR,refdose,a,b){
  # if all parameters are of the same length
  iter <- length(a)
  parslengthNO <- (length(b) != iter)
  if(parslengthNO) stop("Parameters are not of the same length")
  
  yrefdose <- numeric(iter)
  yBMD <- numeric(iter)
  BMD <- numeric(iter)
  for(i in 1:iter){
    # calculate response at refernce dose level
    # notice the refdose is the necessarily standarized
    yrefdose[i] <- get_Linear(refdose,
                            a[i],b[i])
    # hybrid BMR
    yBMD[i] <- (BMR + 1) * yrefdose[i]
    # the length of initial value of x is the same as
    # the length of fitted parameters
    # jac = whether to return Jacobian = NULL
    # BMD_d is estimated by solving the Hill equation 
    # at given parameter values from posterior
    BMD[i] <- nleqslv::nleqslv(x = 0.33,
                               get_Linear,jac = NULL,
                               a[i],b[i])[[1]]
  }
  # BMD is at the same scale as the reference dose
  # BMD is the same length of parameters= iter
  return(BMD)
  
}

getBMD_meta_Hill <- function(BMR,refdose,df_posteriors,df_input){       
  iterEff <- nrow(df_posteriors)
  NStudy <- length(unique(df_input$Index))
  # normalize the reference dose
  mindose <- min(df_input$dose)
  maxdose <- max(df_input$dose)
  refdose_d <- (refdose - mindose)/(maxdose - mindose)  

  # Overarching
  # figure out column number of a and b
  Index_a_overarching <- which(colnames(df_posteriors) == "a")
  Index_b_overarching <- which(colnames(df_posteriors) == "b")
  a_overarching <- df_posteriors[,Index_a_overarching]
  b_overarching <- df_posteriors[,Index_b_overarching]
  c <- df_posteriors$c
  g <- df_posteriors$g
  BMD_d_overarching <- getBMD_rel_Hill(BMR,refdose_d,
                                  a_overarching,
                                  b_overarching,
                                  c,g)
  # transform BMD_d back to original scale
  BMD_overarching <- BMD_d_overarching * (maxdose - mindose) + mindose
  
  # Study_specific
  Index_a_specific <- numeric(iterEff)
  Index_b_specific <- numeric(iterEff)
  df_BMD <- matrix(NA,nrow = iterEff,ncol = NStudy + 1)
  for(s in 1:NStudy){
    Index_a_specific <- which(colnames(df_posteriors) == 
                                paste0("a[",s,"]"))
    Index_b_specific <- which(colnames(df_posteriors) == 
                                paste0("b[",s,"]"))
    a_specific <- df_posteriors[,Index_a_specific]
    b_specific <- df_posteriors[,Index_b_specific]
    BMD_d_specific <- getBMD_rel_Hill(BMR,refdose_d,
                                 a_specific,
                                 b_specific,
                                 c,g)
    # transform BMD_d back to original scale
    df_BMD[,s] <- BMD_d_specific * (maxdose - mindose) + mindose
  }
  
  # combine ovearching and specific
  df_BMD[,NStudy + 1] <- BMD_overarching
  df_BMD <- as.data.frame(df_BMD)
  colnames(df_BMD) <- c(
    paste0("Study",1:NStudy),
    "Overarching"
  )
  print(apply(df_BMD,2,summary))
  save(df_BMD,file = "df_BMD.Rdata")
  return(df_BMD)
}


# Visualization------

# show study-specific curves for Hill partial hierarchical model
showcurve_specific_Hill <- function(df_posteriors,df_input,
                                    aname,bname,cname,gname,
                               xup = NULL, yup = NULL){
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(df_input$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh, length.out = 1E3)
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(df_input$dose)
  maxdose <- max(df_input$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  # doseinput <- dosed
  doseinput <- xdose
  
  # get a blank plot base
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()+
    theme_classic(base_size = 30)
  
  NStudy <- length(unique(df_input$Index))          # number of studies
  c <- df_posteriors[,cname]
  g <- df_posteriors[,gname]
  iterEff <- nrow(df_posteriors)
  a_specific <- numeric(iterEff)
  b_specific <- a_specific
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(xdose))

  # number of dose group in each study
  G <- as.integer(table(df_input$Index))
  
  for(s in 1:NStudy){

    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    
    # Two ways of plotting
    
    # Method 1: dose ~ median of iterative RR by all parameters
    # each iter of pars fits a distribution of y then take descriptives
    # y_specific is the dataframe of y based on study-specific parameters with simulation
    for(i in 1:iterEff){
      # input dosed is standardized
      y_specific[i,] <- get_Hill(doseinput,a = a_specific[i], b = b_specific[i],
                        c = c[i], g = g[i])
    }

    # summary data
    df_plot <- data.frame(x = xdose,
                          # curve of median/median/quantiles
                          median = apply(y_specific,MARGIN =  2,
                                         FUN= function(x) median(x,
                                                                 na.rm = T)),
                          mean = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) mean(x,
                                                             na.rm = T)),
                          Q5 = apply(y_specific,MARGIN =  2,
                                       FUN= function(x) quantile(x,0.05,
                                                                 na.rm = T)),
                          Q95 = apply(y_specific,MARGIN =  2,
                                        FUN= function(x) quantile(x,0.95,
                                                                  na.rm = T)),
                          Q25 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.25,
                                                                 na.rm = T)),
                          Q75 = apply(y_specific,MARGIN = 2,
                                      FUN = function(x) quantile(x,0.75,
                                                                 na.rm = T))
                          )
    
    # # Method 2: dose ~ RR by median parameters
    # y_median <- get_Hill(doseinput, a = median(a_specific),b = median(b_specific),
    #                      c = median(c), g = median(g))
    # y_p5 <- get_Hill(doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05),
    #                  c = quantile(c,0.05), g = quantile(g,0.05))
    # y_p95 <- get_Hill(doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95),
    #                   c = quantile(c,0.95), g = quantile(g,0.95))
    # df_plot <- data.frame(
    #   x = xdose,
    #   # curve of median/5 and 95th percentiles
    #   median = y_median, Q5 = y_p5, Q95 = y_p95
    # )
    # df_temp <- df_input %>% group_by(Study) %>% mutate(
    #   GID = row_number()
    # ) %>% filter(GID != 1)
    df_temp <- df_input
    
    Dose = df_temp$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR = df_temp$RR[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l = df_temp$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u = df_temp$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    
    # data points of each study
    df_rawdata <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    df_rawdata$RR.l[is.na(df_rawdata$RR.l)] <- df_rawdata$RR
    df_rawdata$RR.u[is.na(df_rawdata$RR.u)] <- df_rawdata$RR
    
    Reference <- df_input$Study[(sum(G[1:s-1])+1)]
    xupper <- ifelse(is.null(xup),
                     ceiling(min(
                       max(Dose)*1.1, 
                        max(Dose)+1
                                )
                      ),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling (min(
                       max(RR.u,na.rm = T)*1.1,
                          max(RR.u,na.rm = T)+1
                        )
                     ),
                     yup)
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_rawdata,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_rawdata,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0("study #",s," ",Reference,".svg"),
        width = 10,height = 10)
    print(
      plot_specific_summary
    )
    dev.off()
  }
}

showcurve_overarching_Hill <- function(df_posteriors,df_input,
                                       aname,bname,cname,gname,
                                  xup = NULL, yup = NULL){
  doselow <- 0
  dosehigh <- ceiling(max(df_input$dose))
  xdose <- seq(from = doselow, to = dosehigh, length.out = 1E3)
  # convert xaxis coordinates
  mindose <- min(df_input$dose)
  maxdose <- max(df_input$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose 
  # doseinput <- dosed
  doseinput <- xdose
  
  NStudy <- length(unique(df_input$Index))          # number of studies
  a <- df_posteriors[,aname]
  b <- df_posteriors[,bname]
  c <- df_posteriors[,cname]
  g <- df_posteriors[,gname]
 
  # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(xdose))

  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Hill(doseinput,a[i],b[i],c[i],g[i])
  }

  df_plot <- data.frame(x = xdose,
                        median = apply(y_overarching,MARGIN =  2,
                                       FUN= function(x) median(x,
                                                               na.rm = T)),
                        mean = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) mean(x,
                                                           na.rm = T)),
                        Q5 = apply(y_overarching,MARGIN =  2,
                                     FUN= function(x) quantile(x,0.05,
                                                               na.rm = T)),
                        Q95 = apply(y_overarching,MARGIN =  2,
                                      FUN= function(x) quantile(x,0.95,
                                                                na.rm = T)),
                        Q25 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.25,
                                                              na.rm = T)),
                        Q75 = apply(y_overarching,MARGIN =  2,
                                    FUN= function(x) quantile(x,0.75,
                                                              na.rm = T))
                        )
  
  # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Hill(doseinput,a = median(a), b = median(b), c = median(c),
  #                      g = median(g))
  # y_p5 <- get_Hill(x = doseinput, a = quantile(a,0.05), b = quantile(b, 0.05),
  #                  c = quantile(c, 0.05), g = quantile(g, 0.05))
  # y_p95 <- get_Hill(x = doseinput, a = quantile(a,0.95), b = quantile(b, 0.95),
  #                   c = quantile(c, 0.95), g = quantile(g, 0.95))
  # df_plot <- data.frame(x = xdose,
  #                       median = y_median,Q5 = y_p5,Q95 = y_p95)
  
  df_input$RR.l[is.na(df_input$RR.l)] <- 1
  df_input$RR.u[is.na(df_input$RR.u)] <- 1
  
  xupper <- ifelse(is.null(xup),
                   ceiling (max(
                     median(df_input$dose)*1.1,
                     median(df_input$dose) + 1
                   )),
                   xup
  )

  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(df_input$RR.u,na.rm = T)*1.1,
                       max(df_input$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()+
    theme_classic(base_size = 30)
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "RR")+
    # median curve
    geom_line(data = df_plot,
      aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
    # # 5% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q5),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # # 95% curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = Q95),
    #           color = "grey", size = 0.8,linetype ="dotted")+
    # 25% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q25),
              color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = df_input,
               aes(x = dose,y = RR, group = Study,
                   color = Study, shape = Study),
               size = 5)+
    geom_segment(data = df_input,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = Study, color = Study),
                 size = 1)+
    # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))
  
    
    svg(filename = "Overarching curve.svg",
        width = 10,height = 10)
    suppressWarnings(
      print(plotcurves)
    )
    dev.off()
}

# show study-specific curves for Hill partial hierarchical model
showcurve_specific_Linear <- function(df_posteriors,df_input,
                                    aname,bname,
                                    xup = NULL, yup = NULL){
  # curves as xy plot
  doselow <- 0
  dosehigh <- ceiling(max(df_input$dose))
  # doselow and dosehigh are theoretical range of doses to be as wide as possible
  xdose <- seq(from = doselow, to = dosehigh, length.out = 1E3)
  # convert xaxis coordinates because parameters were fitted on standarized dose
  mindose <- min(df_input$dose)
  maxdose <- max(df_input$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose
  # doseinput <- dosed
  doseinput <- xdose
  
  # get a blank plot base
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()+
    theme_classic(base_size = 30)
  
  NStudy <- length(unique(df_input$Index))          # number of studies
  iterEff <- nrow(df_posteriors)
  a_specific <- numeric(iterEff)
  b_specific <- a_specific
  y_specific <- matrix(NA, nrow = iterEff,ncol = length(xdose))
  
  # number of dose group in each study
  G <- as.integer(table(df_input$Index))
  
  for(s in 1:NStudy){
    
    # calculate study specific y
    a_specific <- df_posteriors[,paste0(aname,"[",s,"]")]
    b_specific <- df_posteriors[,paste0(bname,"[",s,"]")]
    
    # Two ways of plotting
    
    # # Method 1: dose ~ median of iterative RR by all parameters
    # # each iter of pars fits a distribution of y then take descriptives
    # # y_specific is the dataframe of y based on study-specific parameters with simulation
    # for(i in 1:iterEff){
    #   # input dosed is standardized
    #   # y_specific[i,] <- get_Hill(doseinput,a = a_specific[i], b = b_specific[i],
    #   #                            c = c[i], g = g[i])
    #   y_specific[i,] <- get_Linear(doseinput, a = a_specific[i], b = b_specific[i])
    # }
    # 
    # # summary data
    # df_plot <- data.frame(x = xdose,
    #                       # curve of median/median/quantiles
    #                       median = apply(y_specific,MARGIN =  2,
    #                                      FUN= function(x) median(x,
    #                                                              na.rm = T)),
    #                       mean = apply(y_specific,MARGIN =  2,
    #                                    FUN= function(x) mean(x,
    #                                                          na.rm = T)),
    #                       Q5 = apply(y_specific,MARGIN =  2,
    #                                  FUN= function(x) quantile(x,0.05,
    #                                                            na.rm = T)),
    #                       Q95 = apply(y_specific,MARGIN =  2,
    #                                   FUN= function(x) quantile(x,0.95,
    #                                                             na.rm = T)),
    #                       Q25 = apply(y_specific,MARGIN = 2,
    #                                   FUN = function(x) quantile(x,0.25,
    #                                                              na.rm = T)),
    #                       Q75 = apply(y_specific,MARGIN = 2,
    #                                   FUN = function(x) quantile(x,0.75,
    #                                                              na.rm = T))
    # )
    
    # Method 2: dose ~ RR by median parameters
    y_median <- get_Linear(doseinput, a = median(a_specific),b = median(b_specific))
    y_p5 <- get_Linear(doseinput, a = quantile(a_specific,0.05),b = quantile(b_specific,0.05))
    y_p95 <- get_Linear(doseinput, a = quantile(a_specific,0.95),b = quantile(b_specific,0.95))
    df_plot <- data.frame(
      x = xdose,
      # curve of median/5 and 95th percentiles
      median = y_median, Q5 = y_p5, Q95 = y_p95
    )
    df_temp <- df_input %>% group_by(Study) %>% mutate(
      GID = row_number()
    ) %>% filter(GID != 1)
    df_temp <- df_input
    
    Dose = df_temp$dose[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR = df_temp$RR[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.l = df_temp$RR.l[(sum(G[1:s-1])+1):sum(G[1:s])]
    RR.u = df_temp$RR.u[(sum(G[1:s-1])+1):sum(G[1:s])]
    
    # data points of each study
    df_rawdata <- data.frame(
      Dose = Dose,
      RR.l = RR.l,
      RR.u = RR.u,
      RR = RR
    )
    df_rawdata$RR.l[is.na(df_rawdata$RR.l)] <- df_rawdata$RR
    df_rawdata$RR.u[is.na(df_rawdata$RR.u)] <- df_rawdata$RR
    
    Reference <- df_input$Study[(sum(G[1:s-1])+1)]
    xupper <- ifelse(is.null(xup),
                     ceiling(min(
                       max(Dose)*1.1, 
                       max(Dose)+1
                     )
                     ),
                     xup
    )
    
    yupper <- ifelse(is.null(yup),
                     ceiling (min(
                       max(RR.u,na.rm = T)*1.1,
                       max(RR.u,na.rm = T)+1
                     )
                     ),
                     yup)
    plot_specific_summary <- myplotbase+
      # title to display reference
      labs(title = paste0("study #",s," ",Reference),
           x = "ADD(ug/kg/D)",
           y = "RR")+
      # median curve
      geom_line(data = df_plot,
                aes(x = x , y = median),
                color = "blue", size = 1)+
      # # # mean curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = mean),
      #           color = "red", size = 2)+
      # 5% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q5),
                color = "green", size = 1,linetype ="dotted")+
      # 95% curve
      geom_line(data = df_plot,
                aes(x = x , y = Q95),
                color = "red", size = 1,linetype ="dotted")+
      # # 25% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q25),
      #           color = "green", size = 0.8,linetype ="dotted")+
      # # 75% curve
      # geom_line(data = df_plot,
      #   aes(x = x , y = Q75),
      #           color = "red", size = 0.8,linetype ="dotted")+
      # data points median
      geom_point(data = df_rawdata,
                 aes(x = Dose,y = RR),
                 size = 5, color = s)+
      # upper and lower bound of RR
      geom_segment(data = df_rawdata,
                   aes(x = Dose,xend = Dose,
                       y = RR.l,yend = RR.u),
                   size = 1, color= s)+
      # set axis coordinates bound
      coord_cartesian(xlim = c(0,xupper),
                      ylim = c(0,yupper))+
      scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
      scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))
    
    # create figs
    # dev.new("windows",noRStudioGD = T)
    # save as svg
    svg(filename = paste0("study #",s," ",Reference,".svg"),
        width = 10,height = 10)
    print(
      plot_specific_summary
    )
    dev.off()
  }
}

showcurve_overarching_Linear <- function(df_posteriors,df_input,
                                       aname,bname,
                                       xup = NULL, yup = NULL){
  doselow <- 0
  dosehigh <- ceiling(max(df_input$dose))
  xdose <- seq(from = doselow, to = dosehigh, length.out = 1E3)
  # convert xaxis coordinates
  mindose <- min(df_input$dose)
  maxdose <- max(df_input$dose)
  dosed <- (xdose - mindose) / (maxdose - mindose)
  # dosed[dosed < 0] <- 0
  
  # Use original or transformed dose 
  # doseinput <- dosed
  doseinput <- xdose
  
  NStudy <- length(unique(df_input$Index))          # number of studies
  a <- df_posteriors[,aname]
  b <- df_posteriors[,bname]

  # # plotting method 1: dose ~ median of iterative RR by all parameters
  # y_overarching <- matrix(NA, nrow = iterEff,ncol = length(xdose))
  # 
  # # each iter of pars fits a distribution of y then take descriptives
  # for(i in 1:iterEff){
  #   y_overarching[i,] <- get_Linear(doseinput,a[i],b[i])
  # }
  # 
  # df_plot <- data.frame(x = xdose,
  #                       median = apply(y_overarching,MARGIN =  2,
  #                                      FUN= function(x) median(x,
  #                                                              na.rm = T)),
  #                       mean = apply(y_overarching,MARGIN =  2,
  #                                    FUN= function(x) mean(x,
  #                                                          na.rm = T)),
  #                       Q5 = apply(y_overarching,MARGIN =  2,
  #                                  FUN= function(x) quantile(x,0.05,
  #                                                            na.rm = T)),
  #                       Q95 = apply(y_overarching,MARGIN =  2,
  #                                   FUN= function(x) quantile(x,0.95,
  #                                                             na.rm = T)),
  #                       Q25 = apply(y_overarching,MARGIN =  2,
  #                                   FUN= function(x) quantile(x,0.25,
  #                                                             na.rm = T)),
  #                       Q75 = apply(y_overarching,MARGIN =  2,
  #                                   FUN= function(x) quantile(x,0.75,
  #                                                             na.rm = T))
  # )
  
  # # Plotting method 2: dose ~ RR by median parameters
  y_median <- get_Linear(doseinput,a = median(a), b = median(b))
  y_p5 <- get_Linear(x = doseinput, a = quantile(a,0.05), b = quantile(b, 0.05))
  y_p95 <- get_Linear(x = doseinput, a = quantile(a,0.95), b = quantile(b, 0.95))
  y_p25 <- get_Linear(x = doseinput, a = quantile(a,0.25), b = quantile(b, 0.25))
  y_p75 <- get_Linear(x = doseinput, a = quantile(a,0.75), b = quantile(b, 0.75))
  df_plot <- data.frame(x = xdose,
                        median = y_median,Q5 = y_p5,Q95 = y_p95,
                        Q25 = y_p25,Q75 = y_p75)
  
  df_input$RR.l[is.na(df_input$RR.l)] <- 1
  df_input$RR.u[is.na(df_input$RR.u)] <- 1
  
  xupper <- ifelse(is.null(xup),
                   ceiling (max(
                     median(df_input$dose)*1.1,
                     median(df_input$dose) + 1
                   )),
                   xup
  )
  
  yupper <- ifelse(is.null(yup),
                   ceiling(
                     min(
                       max(df_input$RR.u,na.rm = T)*1.1,
                       max(df_input$RR.u,na.rm = T) +1
                     )
                   ),
                   yup)
  graphics.off()
  myplotbase <- ggplot(data.frame())+ geom_blank()+
    theme_classic(base_size = 30)
  
  plotcurves <- myplotbase+
    # title to display reference
    labs(title = "Overarching Curve",
         x = "ADD(ug/kg/D)",
         y = "RR")+
    # median curve
    geom_line(data = df_plot,
              aes(x = x , y = median),
              color = "blue", size = 2)+
    # # mean curve
    # geom_line(data = df_plot,
    #           aes(x = x , y = mean),
    #           color = "red", size = 2)+
  # # 5% curve
  # geom_line(data = df_plot,
  #           aes(x = x , y = Q5),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # # 95% curve
  # geom_line(data = df_plot,
  #           aes(x = x , y = Q95),
  #           color = "grey", size = 0.8,linetype ="dotted")+
  # 25% curve
  geom_line(data = df_plot,
            aes(x = x , y = Q25),
            color = "brown", size = 0.8,linetype ="dotted")+
    # 75% curve
    geom_line(data = df_plot,
              aes(x = x , y = Q75),
              color = "brown", size = 0.8,linetype ="dotted")+
    # data points
    geom_point(data = df_input,
               aes(x = dose,y = RR, group = Study,
                   color = Study, shape = Study),
               size = 5)+
    geom_segment(data = df_input,
                 aes(x = dose,xend = dose,
                     y = RR.l,yend = RR.u,
                     group = Study, color = Study),
                 size = 1)+
    # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))
  
  
  svg(filename = "Overarching curve.svg",
      width = 10,height = 10)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
}
