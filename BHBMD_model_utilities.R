# ----Settings----
#* Environment----
if(exists("seed")){
  stanseed = seed
  set.seed(seed = seed)
  } else {
  stanseed = round(rnorm(1,0,100),0)
}

library("rstan")
library("rstanarm")
library("loo")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library("shinystan")
library("later")
later:::ensureInitialized()

options(scipen = 10)
library("dplyr")
library("ggplot2")
library("readxl")
library("tidyr")
library("ggpubr")
library("reshape2")


iterEff <- (samplingiter*(1-warmup)*Nchain)

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
  CV <- ysd / ymean
  ysdL <- sqrt(log(CV^2 +1))   # based on summarized ymean and ysd 
  return(ysdL)                         # at regular scale
}

get_ymeanL <- function(ymean,ysd){   # estimate summarized ymean at log scale
  CV <- ysd / ymean
  ymeanL <- log(ymean/sqrt(1+CV^2))
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

get_Power <- function(x,a,b,g){
  # calculate response given Power parameters
  # can accomendate raw and normalized dose
  # but normalized dose recommended
  # parameters must be real numbers
  if(length(a) != 1) stop("length of a not 1") 
  if(length(b) != 1) stop("length of b not 1")
  if(length(g) != 1) stop("length of g not 1")
  x <- unlist(x)                # x can be a vector
  
  y <- numeric(length = length(x))
  y <- a + b * Re(as.complex(x)^g)
  return(y)
}

get_Expo5 <- function(x,a,b,c,d){
  # calculate response given Power parameters
  # can accomendate raw and normalized dose
  # but normalized dose recommended
  # parameters must be real numbers
  if(length(a) != 1) stop("length of a not 1") 
  if(length(b) != 1) stop("length of b not 1")
  if(length(c) != 1) stop("length of c not 1") 
  if(length(d) != 1) stop("length of d not 1")
  x <- unlist(x)                # x can be a vector
  
  y <- numeric(length = length(x))
  y <- a * (c - (c-1) * exp(-(b*x)^d))
  return(y)
}

get_MM <- function(x,a,b,c){
  # calculate response given Power parameters
  # can accomendate raw and normalized dose
  # but normalized dose recommended
  # parameters must be real numbers
  if(length(a) != 1) stop("length of a not 1") 
  if(length(b) != 1) stop("length of b not 1")
  if(length(c) != 1) stop("length of c not 1") 
  x <- unlist(x)                # x can be a vector
  
  y <- numeric(length = length(x))
  y <- a + b * x / (c + x)
  return(y)
}

getRRSE <- function(RR.l,RR.u){
  SE <- (log(RR.u)-log(RR.l))/2/1.96
  SE[is.nan(SE)] <- 0
  return(SE)
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

fit_Linear <- function(pars,input,y){
  yhat <- pars[1] + pars[2] * input
  SSE <- sum((yhat - y)^2)
  if(is.finite(SSE)) {return(SSE)} else {
    return(1000 + sum(pars^2))
  } 
}

fit_Hill <- function(pars,input,y){
  yhat <- pars[1] + (pars[2] * input^pars[4]) / (pars[3]^pars[4] + input^pars[4])
  SSE <- sum((yhat - y)^2)
  if(is.finite(SSE)) {return(SSE)} else {
    return(1000 + sum(pars^2))
  } 
}

fit_Power <- function(pars,input,y){
  yhat <- pars[1] + pars[2] * (input^pars[3])
  SSE <- sum((yhat - y)^2)
  if(is.finite(SSE)) {return(SSE)} else {
    return(1000 + sum(pars^2))
  } 
}

fit_MM <- function(pars,input,y){
  yhat <- pars[1] + pars[2] * input / (pars[3] + input)
  SSE <- sum((yhat - y)^2)
  if(is.finite(SSE)) {return(SSE)} else {
    return(1000 + sum(pars^2))
  } 
}

fit_Expo5 <- function(pars,input,y){
  yhat <- pars[1] * (pars[3] - (pars[3]-1) * exp(-1 * (pars[2] * input)^pars[4]))
  SSE <- sum((yhat - y)^2)
  if(is.finite(SSE)) {return(SSE)} else {
    return(1000 + sum(pars^2))
  } 
}
  
# *Data Treatment-------

# Simulate meta DR data based on four DR models
fn_metaDRsimulation <- function(S = 10 , G = 3,Sub = c(250,210,180),
                                Dose = c(10,50,90),
                                ybound = c(exp(-10),20),
                                sigmavalue = 0.2){
  N <- S*G
  Index <- rep(1:S,each = G)
  Study <- LETTERS[Index]
  # alldose <- replicate(S,seq(from = 0, to = max(dose), lengt = 1E3))
  Nsub <- ceiling(replicate(S,Sub+runif(G,0,60)))
  dose <- replicate(S,Dose + runif(G,-1,1)*10)
  x_bound <- c(min(dose),max(dose))
  y_bound <- sort(runif(2,min = ybound[1],max = ybound[2]))
  noisesize <- 0.1
  
  # Linear
  popfit_linear <- optim(c(1,1),fit_Linear, input = x_bound, y = y_bound,
                         method = c("L-BFGS-B"),lower = c(0,0))
  pars_linear <- popfit_linear$par
  # get_Linear(x_bound,pars_linear[1],pars_linear[2])
  a_Linear <- pars_linear[1] + rnorm(S,0,noisesize * pars_linear[1])
  b_Linear <- pars_linear[2] + rnorm(S,0,noisesize * pars_linear[2])
  
  # Hill
  popfit_Hill <- optim(c(1,1,50,1),fit_Hill, input = x_bound, y = y_bound,
                       method = c("L-BFGS-B"),lower = c(0,0,0,0),upper = c(Inf,Inf,Inf,18))
  pars_Hill <- popfit_Hill$par
  # get_Hill(x_bound,pars_Hill[1],pars_Hill[2],pars_Hill[3],pars_Hill[4])
  a_Hill <- pars_Hill[1] + rnorm(S,0,noisesize * pars_Hill[1])
  b_Hill <- pars_Hill[2] + rnorm(S,0,noisesize * pars_Hill[2])
  c_Hill <- pars_Hill[3] + rnorm(S,0,noisesize * pars_Hill[3])
  g_Hill <- pars_Hill[4] + rnorm(S,0,noisesize * pars_Hill[4])
  
  # Power
  popfit_Pow <- optim(c(1,1,1),fit_Power, input = x_bound, y = y_bound,
                      method = "L-BFGS-B",lower = c(0,0,0), upper = c(Inf,Inf,18))
  pars_Pow<- popfit_Pow$par
  # get_Power(x_bound,pars_Pow[1],pars_Pow[2],pars_Pow[3])
  a_Pow <- pars_Pow[1] + rnorm(S,0,noisesize * pars_Pow[1])
  b_Pow <- pars_Pow[2] + rnorm(S,0,noisesize * pars_Pow[2])
  g_Pow <- pars_Pow[3] + rnorm(S,0,noisesize * pars_Pow[3])
  
  # MM
  popfit_MM <- optim(c(0.5,0.5,0.5),fit_MM, input = x_bound, y = ybound,
                     method = "L-BFGS-B",lower = c(0,0,0)
  )
  pars_MM <- popfit_MM$par
  a_MM <- pars_MM[1] + rnorm(S,0,noisesize * pars_MM[1])
  b_MM <- pars_MM[2] + rnorm(S,0,noisesize * pars_MM[2])
  c_MM <- pars_MM[3] + rnorm(S,0,noisesize * pars_MM[3])
  
  # Exponential 5
  popfit_Expo5 <- optim(c(0.5,0.5,1,1),fit_Expo5, input = x_bound, y = y_bound,
                        method = "L-BFGS-B",lower = c(0,0,1,1), upper = c(Inf,Inf,Inf,18))
  pars_Expo5 <- popfit_Expo5$par
  # get_Expo5(x_bound,pars_Expo5[1],pars_Expo5[2],pars_Expo5[3],pars_Expo5[4])
  a_Expo5 <- pars_Expo5[1] + rnorm(S,0,noisesize * pars_Expo5[1])
  b_Expo5 <- pars_Expo5[2] + rnorm(S,0,noisesize * pars_Expo5[2])
  c_Expo5 <- pars_Expo5[3] + rnorm(S,0,noisesize * pars_Expo5[3])
  d_Expo5 <- pars_Expo5[4] + rnorm(S,0,noisesize * pars_Expo5[4])
  
  # Plotting and simulation fd
  fd_Linear <- vector("list",S)
  fd_Hill <- vector("list",S)
  fd_MM <- vector("list",S)
  fd_Power <- vector("list",S)
  fd_Expo5 <- vector("list",S)
  
  xdose <- dose           # use dose for simulation/ alldose for plotting
  
  # plot(0,xlim = c(0,max(dose)), ylim = c(0,y_bound[2]))
  # for(s in 1:S){
  #   fd_Linear[[s]] <- get_Linear(xdose[,s],a = a_Linear[s], b = b_Linear[s])
  #   lines(xdose[,s],fd_Linear[[s]],col = "blue")
  #   fd_Hill[[s]] <- get_Hill(xdose[,s],a = a_Hill[s], b = b_Hill[s], c = c_Hill[s], g = g_Hill[s])
  #   lines(xdose[,s],fd_Hill[[s]], col = "red")
  #   fd_Power[[s]] <- get_Power(xdose[,s], a = a_Pow[s], b = b_Pow[s], g = g_Pow[s])
  #   lines(xdose[,s],fd_Power[[s]], col = "green")
  #   fd_Expo5[[s]] <- get_Expo5(xdose[,s], a = a_Expo5[s], b = b_Expo5[s], c = c_Expo5[s], d = d_Expo5[s])
  #   lines(xdose[,s],fd_Expo5[[s]],col = "yellow")
  # }
  for(s in 1:S){
    fd_Linear[[s]] <- get_Linear(xdose[,s],a = a_Linear[s], b = b_Linear[s])
    fd_Hill[[s]] <- get_Hill(xdose[,s],a = a_Hill[s], b = b_Hill[s], c = c_Hill[s], g = g_Hill[s])
    fd_Power[[s]] <- get_Power(xdose[,s], a = a_Pow[s], b = b_Pow[s], g = g_Pow[s])
    fd_MM[[s]] <- get_MM(xdose[,s], a = a_MM[s], b = b_MM[s], c = c_MM[s])
    fd_Expo5[[s]] <- get_Expo5(xdose[,s], a = a_Expo5[s], b = b_Expo5[s], c = c_Expo5[s], d = d_Expo5[s])
  }
  
  y_Linear <- vector("list",N)
  y_Hill <- vector("list",N)
  y_MM <- vector("list",N)
  y_Power <- vector("list",N)
  y_Expo5 <- vector("list",N)
  for(g in 1:N){
    y_Linear[[g]] <- rlnorm(Nsub[g],meanlog = log(fd_Linear[[Index[g]]][g-3*(Index[g]-1)]),sdlog = sigmavalue)
    y_Hill[[g]] <- rlnorm(Nsub[g],meanlog = log(fd_Hill[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
    y_Power[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_Power[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
    y_Expo5[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_Expo5[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
    y_MM[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_MM[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
  }
  y_Linear_summary <- lapply(y_Linear,function(x){
    cbind(mean(x),sd(x))
  })
  y_Hill_summary <- lapply(y_Hill,function(x){
    cbind(mean(x),sd(x))
  })
  y_Power_Summary <- lapply(y_Power,function(x){
    cbind(mean(x),sd(x))
  })
  y_MM_Summary <- lapply(y_MM,function(x){
    cbind(mean(x),sd(x))
  })
  y_Expo5_Summary <- lapply(y_Expo5, function(x){
    cbind(mean(x),sd(x))
  })
  
  df_simu_meta_Linear <- data.frame(
    Index = Index,
    Study = Study,
    dose = as.vector(dose),
    Nsub = as.vector(Nsub),
    ymean = unlist(lapply(y_Linear_summary,function(x) x[,1])),
    ysd = unlist(lapply(y_Linear_summary,function(x) x[,2]))
  )
  
  df_simu_meta_Hill <- data.frame(
    Index = Index,
    Study = Study,
    dose = as.vector(dose),
    Nsub = as.vector(Nsub),
    ymean = unlist(lapply(y_Hill_summary,function(x) x[,1])),
    ysd = unlist(lapply(y_Hill_summary,function(x) x[,2]))
  )
  
  df_simu_meta_Power <- data.frame(
    Index = Index,
    Study = Study,
    dose = as.vector(dose),
    Nsub = as.vector(Nsub),
    ymean = unlist(lapply(y_Power_Summary,function(x) x[,1])),
    ysd = unlist(lapply(y_Power_Summary,function(x) x[,2]))
  )
  df_simu_meta_MM <- data.frame(
    Index = Index,
    Study = Study,
    dose = as.vector(dose),
    Nsub = as.vector(Nsub),
    ymean = unlist(lapply(y_MM_Summary,function(x) x[,1])),
    ysd = unlist(lapply(y_MM_Summary,function(x) x[,2]))
  )  
  df_simu_meta_Expo5 <- data.frame(
    Index = Index,
    Study = Study,
    dose = as.vector(dose),
    Nsub = as.vector(Nsub),
    ymean = unlist(lapply(y_Expo5_Summary,function(x) x[,1])),
    ysd = unlist(lapply(y_Expo5_Summary,function(x) x[,2]))
  )
  
  baseinfo <- list(
    x_bound,y_bound,
    a_Linear,b_Linear,
    a_Pow,b_Pow,g_Pow,
    a_MM,b_MM,c_MM,
    a_Hill,b_Hill,c_Hill,g_Hill,
    a_Expo5,b_Expo5,c_Expo5,d_Expo5
  )
  ls_data <- list(
    Linear = df_simu_meta_Linear,
    Power = df_simu_meta_Power,
    MM = df_simu_meta_MM,
    Hill = df_simu_meta_Hill,
    Expo5 = df_simu_meta_Expo5,
    baseinfo = baseinfo
  )
  return(ls_data)
}

# Extract Simulated Data by Model
fn_getsimuData <- function(ls_data,DSName){
  df_input <- as.data.frame(ls_data[DSName])
  colnames(df_input) <- stringr::str_replace(colnames(df_input),pattern = paste0(DSName,"."),replacement = "")
  return(df_input)
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

PrepareRRData <- function(df_input){
  Index <- df_input$Index
  Study <- df_input$Study
  dose <- df_input$dose
  Nsub <- df_input$Nsub
  ymean <- df_input$RR
  
  ymeanL <- log(ymean)
  ysdL <- getRRSE(RR.l = df_input$RR.l,RR.u = df_input$RR.u)
  ysd <- get_ysd(ymean,ysdL)
  doseZ <- dnormalize_minmax(dose)
  
  N <- nrow(df_input)
  S <- length(unique(df_input$Index))
  Gcount <- df_input %>% group_by(Index) %>% count()
  G <- as.data.frame(ungroup(Gcount))[,2]
  Data <- list(N = N,S = S,G = G,
               Index = Index,Study = Study,Nsub = Nsub,
               ymean = ymean,ymeanL = ymeanL,
               ysd = ysd, ysdL = ysdL,
               doseZ = doseZ,dose = dose)
  return(Data)
}

PrepareNRRLnormData <- function(df_input){
  Index <- df_input$Index
  Study <- df_input$Study
  dose <- df_input$dose
  Nsub <- df_input$Nsub
  ymean <- df_input$ymean
  ysd <- df_input$ysd
  
  ymeanL <- get_ymeanL(ymean,ysd)
  ysdL <- get_ysdL(ymean,ysd)
  doseZ <- dnormalize_minmax(dose)
  
  N <- nrow(df_input)
  S <- length(unique(df_input$Index))
  Gcount <- df_input %>% group_by(Index) %>% count()
  G <- as.data.frame(ungroup(Gcount))[,2]
  Data <- list(N = N,S = S,G = G,
               Index = Index,Study = Study,Nsub = Nsub,
               ymean = ymean,ymeanL = ymeanL,
               ysd = ysd,  ysdL = ysdL,
               doseZ = doseZ,dose = dose
               )
  return(Data)
}

PrepareLnormData <- function(df_input, RR=F){
  if(RR){
    PreparedData <- PrepareRRData(df_input)
  } else {
    PreparedData <- PrepareNRRLnormData(df_input)
  }
  return(PreparedData)
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
get_metaStanfit <- function(modelname,dataName,DataList,
                            inits = NULL){
  # Priors
  # Calculate data-based control values for hyperparameters
  ybound <- c(
    max(DataList$ymean),
    sd(DataList$ymean)
  )
  a_bound <- log(ybound[1] + 3 * ybound[2])*2
  Xrange <- sort(c(max(DataList$doseZ),min(DataList$doseZ)))
  xbound <- c(
    Xrange[1],
    Xrange[2]
  )
  
  source(paste0(fileDir,"/BHBMD_priors.R"),local = T)

  inits <- ifelse(is.null(inits),"random",inits)
  
  # Compile stanmodels
  # assign(paste0("model_",modelname), 
  #        # stan_model(file = paste0("stanmodel_",modelname,".stan"))
  #        stan_model(model_code = eval(parse(text = paste0("modelstring_",modelname))))
  # )
  
  assign(paste0("model_",modelname),
         stan_model(file = paste0(fileDir, modelname,".stan")))
  
  SaveFileName <- paste0("stanfit_",modelname,"_",dataName)
  
  assign(SaveFileName,
         rstan::sampling(
           object =  eval(parse(text = paste0("model_",
                                              modelname))),
           data = append(DataList,
                         eval(parse(text = paste0("prior_",
                                                  modelname)))
           ),
           iter = samplingiter,
           chains = Nchain,
           warmup = warmup * samplingiter,
           thin = thin,
           seed = stanseed,
           init = inits,
           # init_r = 2,
           control = list(
             # adapt_delta  = 0.99,
             # max_treedepth = 12
           ),
           cores = parallel::detectCores()
         ))
  # save(list = SaveFileName,file = paste0(SaveFileName,".rdata"))
  saveRDS(
    object = eval(parse(text = paste0(SaveFileName))),
    file = paste0(SaveFileName,".rds")
  )
  return(eval(parse(text = SaveFileName)))
}

# wrapper to fit all four dose response models given a DataList
fit_allmodels <- function(DSName,DataList){
  dataName <- paste0("simu_",DSName)
  
  modelname <- "Lognormal_Summary_Linear_meta_3"
  assign(paste0("stanfit_",modelname,"_",dataName),
         get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
  modelname <- "Lognormal_Summary_Power_meta_4"
  assign(paste0("stanfit_",modelname,"_",dataName),
         get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
  modelname <- "Lognormal_Summary_Hill_meta_5"
  assign(paste0("stanfit_",modelname,"_",dataName),
         get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
  modelname <- "Lognormal_Summary_Expo5_meta_5"
  assign(paste0("stanfit_",modelname,"_",dataName),
         get_metaStanfit(modelname = modelname,dataName = dataName,DataList = DataList))
  
  assign(
    paste0("fitList_",dataName),
    list(
      Linear = eval(parse(text = paste0("stanfit_Lognormal_Summary_Linear_meta_3_",dataName))),
      Hill = eval(parse(text = paste0("stanfit_Lognormal_Summary_Hill_meta_5_",dataName))),
      Power = eval(parse(text = paste0("stanfit_Lognormal_Summary_Power_meta_4_",dataName))),
      Expo5 = eval(parse(text = paste0("stanfit_Lognormal_Summary_Expo5_meta_5_",dataName)))
    )
  )
  return(eval(parse(text = paste0("fitList_",dataName))))
}

# wrapper to fit four dose response models to a list of simulated data
fit_simuData <- function(simuData,DSName){
  df_input <- fn_getsimuData(simuData,DSName)
  DataList <- PrepareLnormData(df_input,RR=F)
  dataName <- paste0("simu_",DSName)
  
  fitList <- fit_allmodels(DSName = DSName,DataList = DataList)
  return(fitList)
}

get_fitsummary <- function(stanfit){
  print(stanfit)
  summary <- as.data.frame(summary(stanfit)[1])
  # keep only stats in the name
  colnames(summary) <- substring(colnames(summary),9,20)
  # write.csv(summary,file = "posterior summary.csv")
  # modelname <- substr(deparse(substitute(stanfit)),9,100)
  df_posteriors <- as.data.frame(stanfit)
  write.csv(df_posteriors,
            row.names = F,
            file = paste0("Posteriors_",modelname,"_",dataName,".csv"))
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

# Model Evaluation-------
# getNpars <- function(stanfit){
#   TryNpars <- try(get_num_upars(stanfit),silent = T)
#   if("try-error" %in% class(TryNpars)){
#     fit <- rstan::extract(stanfit)
#     Nfit <- dim(as.data.frame(fit))[2]
#     Nobs <- dim(fit$log_lik)[2]
#     Npars <- Nfit - Nobs * 2 - 1
#   } else {
#     Npars <- get_num_upars(stanfit)
#   }
#   Npars
# }

getAIC <- function(log_lik,pars){
  if(is.array(log_lik)){
    # log_lik_array <- log_lik
    log_lik_mtx <- apply(log_lik,3,c)
  } else {
    stop("Log likelihood must be an array")
  }
  # dev <- apply(log_lik_array,c(1,2),sum)
  dev <- apply(log_lik_mtx,1,sum)
  AIC <- dev - pars
  # AIC_mean <- mean(AIC)
  # return(AIC_mean)
  return(AIC)
}

getBIC <- function(log_lik,pars){
  if(is.array(log_lik)){
    log_lik_mtx <- apply(log_lik,3,c)
  } else {
    stop("Log likelihood must be an array")
  }
  # dev <- apply(log_lik_array,c(1,2),sum)
  dev <- apply(log_lik_mtx,1,sum)
  BIC <- dev - 0.5 * pars * log(dim(log_lik)[3]) 
  # BIC_mean <- mean(BIC)
  # return(BIC_mean)\
  return(BIC)
}

getLPML <- function(log_lik){
  if(is.array(log_lik)){
    Niter <- dim(log_lik)[1] * dim(log_lik)[2]
    Nobs <- dim(log_lik)[3]
    LL_mtx <- apply(log_lik,3,c)
  } else {
    stop("Log likelihood must be an array")
  }
  LCPO <- rep(NA,Nobs)
  for(n in 1:Nobs){
    LL_max <- max(LL_mtx[,n])
    LL_Z <- LL_mtx[,n] - LL_max
    LCPO[n] <- log(Niter)+LL_max-log(sum(exp(-LL_Z)))
  }
  LPML <- sum(LCPO)
  return(LPML)
}
getwts <- function(fitList){
  Nmod <- length(fitList)
  Npars <- as.numeric(unlist(lapply(fitList,get_num_upars)))
  log_lik_array <- lapply(fitList,FUN = function(x) extract_log_lik(x,merge_chains = F))
  log_lik_list <- lapply(fitList,extract_log_lik)
  Nobs <- dim(log_lik_array[[1]])[3]
  Niter <- nrow(log_lik_list[[1]])
  r_eff_list <- lapply(fitList,function(x) {
    ll_array <- extract_log_lik(x,merge_chains = F)
    relative_eff(exp(ll_array))
  })
  LPML_list <- unlist(lapply(log_lik_array,getLPML))
  
  AIC_list <- matrix(NA,nrow = Niter,ncol = Nmod)
  BIC_list <- matrix(NA,nrow = Niter,ncol = Nmod)
  loo_list <- vector("list",Nmod)
  lpd_point <- matrix(NA,nrow = Nobs,ncol = Nmod)
  WAIC_list <- vector("list",Nmod)
  elpd_waic <- matrix(NA,Nobs,Nmod)
  # elpd_loo <- rep(NA,Nmod)
  
  for(i in 1:Nmod){
    AIC_list[,i] <- getAIC(log_lik = log_lik_array[[i]],pars = Npars[[i]])
    BIC_list[,i] <- getBIC(log_lik_array[[i]],Npars[[i]])
    loo_list[[i]] <- loo(log_lik_list[[i]],r_eff = r_eff_list[[i]])
    lpd_point[,i] <- loo_list[[i]]$pointwise[, "elpd_loo"] 
    # elpd_loo[i] <- loo_list[[i]]$estimates["elpd_loo", "Estimate"]
    WAIC_list[[i]] <- waic(log_lik_list[[i]])
    elpd_waic[,i] <- WAIC_list[[i]]$pointwise[,1]
  }
  
  lpd_point_Z <- lpd_point - apply(lpd_point,1,max)
  
  
  # AIC weights
  AIC_Z <- AIC_list - apply(AIC_list,1,max)
  AIC_wts <- apply(exp(AIC_Z) / apply(exp(AIC_Z),1,sum),2,mean)

  # BIC weights
  BIC_Z <- BIC_list - apply(BIC_list,1,max)
  BIC_wts <- apply(exp(BIC_Z) / apply(exp(BIC_Z),1,sum),2,mean)
  
  # WAIC weights
  # WAIC_list <- lapply(log_lik_list,waic)
  # WAICValue_list <- lapply(WAIC_list,function(x){
  #   x$estimates["elpd_waic",1]
  # })
  # WAIC_Z <- unlist(WAICValue_list) - max(unlist(WAICValue_list))
  # WAIC_wts <- exp(WAIC_Z)/sum(exp(WAIC_Z))
  elpd_waic_Z <- elpd_waic - apply(elpd_waic,1,max)
  WAIC_wts <- apply(exp(elpd_waic_Z) / apply(exp(elpd_waic_Z),1,sum),2,mean)
  
  # LPML weights
  LPML_Z <- LPML_list - max(LPML_list)
  LPML_Wts <- exp(LPML_Z) / sum(exp(LPML_Z))
  
  # loo
  wts_stacking <- stacking_weights(lpd_point_Z)
  wts_pbma <- pseudobma_weights(lpd_point_Z,BB=F)
  wts_pbmabb <- pseudobma_weights(lpd_point_Z,BB=T)
  
  # Weights
  mtx_wts <- round(cbind(
    AIC = AIC_wts,
    BIC_BMA = BIC_wts,
    WAIC = WAIC_wts,
    LPML = LPML_Wts,
    Stacking = wts_stacking,
    PseudoBMA = wts_pbma,
    PseudoBMA_BB = wts_pbmabb
  ),7)  
  return(mtx_wts)
}

MatchTrueModel <- function(modelweights,dataName){
  require("stringr")

  dataLabel <- str_replace(dataName,pattern = "simu_",replacement = "")
  TrueMIndex <- which(rownames(modelweights) == dataLabel)
  
  
  # consider correct if the true model is given the max weight
  TrueByMaxWt <- apply(modelweights,MARGIN = 2,function(x){
    x[TrueMIndex] == max(x)
  })
  
  # consider correct if the true model weight >0.9
  TrueBy90Wt <- apply(modelweights,MARGIN = 2,function(x){
    x[TrueMIndex] >= 0.9
  })
}

CalculateMetrics <- function(baseline,weighted,
                             Metrics = c("RMSE","MAE")){
  
  ls_compare <- weighted
  if(!exists("Method_list")){
    Method_list <- names(ls_compare)
  } 
  
  y_error <- lapply(ls_compare,function(x){
    t(apply(x,MARGIN = 1,function(x){
      x - baseline
    }))
  })
  
  metric_list <- Metrics
  Nmetrics <- length(metric_list)
  
  RMSE_iter <- lapply(y_error,function(x) {
    sqrt(apply(x,MARGIN = 1,function(x) mean(x^2)))})
  RMSE <- lapply(RMSE_iter,mean)
  
  MAE_iter <- lapply(y_error,function(x) {
    apply(x,MARGIN = 1,function(x) mean(abs(x)))
  })
  MAE <- lapply(MAE_iter,mean)
  
  # MRB_iter <- lapply(y_weighted,function(x){
  #   apply(x,MARGIN = 1,function(x){abs(x - y_original) / y_original})
  # })
  # MRB <- lapply(MRB_iter,median)
    
  CompareMethods <- as.data.frame(rbind(
    unlist(RMSE),
    unlist(MAE)
    # unlist(MRB)
  ))
  colnames(CompareMethods) <- Method_list
  rownames(CompareMethods) <- metric_list
  
  return(CompareMethods)
}

ModelEvalbyFd <- function(fitList,modelweights,DataList){
  if(exists("modelweights")){
    wts_1 <- modelweights
  } else {
    wts_1 <- getwts(fitList)
  }
  input_quantity <- DataList$dose
  y_original <- DataList$ymean
  
  Method_list <- colnames(wts_1)       # methods of averaging
  modelnames <- rownames(modelweights) # models used to fit data
  Nstan <- length(fitList)            # number of fitted models
  Nmethods <- ncol(wts_1)             # number of methods
  
  # stanfit_Linear <- fitList$Linear
  # stanfit_Hill <- fitList$Hill
  # stanfit_Pow <- fitList$Power
  # stanfit_Expo5 <- fitList$Expo5
  
  fd_list <- vector("list",Nstan)
  fd_list <- lapply(fitList,function(x) rstan::extract(x)$fd)

  # fd_Linear <- extract(stanfit_Linear)$fd
  # fd_Hill <- extract(stanfit_Hill)$fd
  # fd_Pow <- extract(stanfit_Pow)$fd
  # fd_Expo5 <- extract(stanfit_Expo5)$fd
  
  Nobs <- dim(fd_list[[1]])[2]
  Niter <- dim(fd_list[[1]])[1]          # stan iterations
  
  y_weighted <- vector("list",length = Nmethods)
  df_mean_weighted <- matrix(NA,nrow = Nmethods,ncol = Nobs)

  
  for(m in 1:Nmethods){
    y_weighted[[m]] <- matrix(0,Niter,Nobs)
    colnames(y_weighted[[m]]) <- input_quantity
    for(s in 1:Nstan){
      y_weighted[[m]] <- y_weighted[[m]] + fd_list[[s]] * wts_1[s,m] 
    }
    df_mean_weighted[m,] <- apply(y_weighted[[m]],MARGIN = 2,mean)
  }
  df_mean_weighted <- rbind(df_mean_weighted,y_original)
  rownames(df_mean_weighted) <- c(Method_list,"Raw Values")
  colnames(df_mean_weighted) <- input_quantity
  
  Metric_list <- c("RMSE","MRB")
  Nmetrics <- 2
  CompareMethods <- CalculateMetrics(baseline = y_original,
                                     weighted = y_weighted,
                                     Metrics = Metric_list)
  colnames(CompareMethods) <- Method_list
  
  MethodChoice <- replicate(Nmethods,rep(0,Nmetrics))
  min_Index <- apply(CompareMethods,1,function(x){
    which(x == min(x))
  })
  for(m in 1:Nmetrics){
    MethodChoice[m,c(min_Index[[m]])] <- 1
  }
  colnames(MethodChoice) <- Method_list
  rownames(MethodChoice) <- Metric_list
  
  # min_Index <- apply(CompareMethods,1,function(x) {
  #   which(x == min(x))
  #   # paste0("Min= ",round(min(x),5)," in ", paste(Method_list[which(x == min(x))],collapse = " & ") )
  #   paste(Method_list[which(x == min(x))],collapse = " & ")
  # })
  # CompareMethods$Choice <- min_Index
  # 
  
  return(list(
      WeightedQuantities = y_weighted,
      WeightedMean = df_mean_weighted,
      ModelMetrics = CompareMethods,
      ModelChoice = MethodChoice
  ))
} 


# wrapper for model averaging evaluation based on simulated data
MAS_iter <- function(DSName,data_iter){
  # Prepare data into datalist for rstan
  DataList_iter <- lapply(data_iter,function(x) 
    PrepareLnormData(fn_getsimuData(x,DSName = DSName),RR=F))
  # Export datalist as csv files per iteration
  mapply(write.csv,x = DataList_iter, file = paste0("DataList_",DSName,"_",1:iter,".csv"))
  # stan fitting
  fitList_iter <- lapply(data_iter,function(x) fit_datalist(x,DSName = DSName))
  # export stanfit objects per iteration
  mapply(save,list = "fitList_iter",file = paste0("fitList_",DSName,"_",1:iter,".rdata"))
  
  weights_iter <- lapply(fitList_iter, getwts)
  # Export weightings
  mapply(write.csv,x = weights_iter, file = paste0("Weights_",DSName,"_",1:iter,".csv"))
  # Performance evaluation
  PerfMetrics_iter <- t(mapply(ModelEvalbyFd,      # a iter * 4 matrix
                               fitList = fitList_iter,
                               modelweights = weights_iter,
                               DataList = DataList_iter))
  Weighted_iter <- PerfMetrics_iter[,1]
  mapply(write.csv,x = Weighted_iter, file = paste0("Weighted_iter_",DSName,"_",1:iter,".csv"))
  
  Weighted_mean <- PerfMetrics_iter[,2]
  mapply(write.csv,x = Weighted_mean, file = paste0("Weighted_mean_",DSName,"_",1:iter,".csv"),row.names =F)
  
  Metrics <- PerfMetrics_iter[,3]
  mapply(write.csv,x = Metrics, file = paste0("Metrics_",DSName,"_",1:iter,".csv"),row.names =F)
  
  ChoiceList <- PerfMetrics_iter[,4]
  SelectionCount <- do.call("+",ChoiceList)
  write.csv(SelectionCount,file = paste0("Selection_Count_",DSName,".csv"),row.names =F)
  
  options(bitmapType='cairo')
  mapply(ShowWeightedValues,
         DataList = DataList_iter,
         modelwts = weights_iter,
         weighted = Weighted_iter,
         dataName = paste0("simu_",DSName,"_",1:iter)
  )
  
}

# Presentation------

# show study-specific curves for Hill partial hierarchical model
showcurve_specific_Hill <- function(df_posteriors,DataList,
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
  
  # Use original or standarized dose 
  doseinput <- dosed
  # doseinput <- xdose
  
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
  
  
  
  # Plotting
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
    # set axis coordinates bound
    coord_cartesian(xlim = c(0,xupper),
                    ylim = c(0,yupper))+
    scale_x_continuous(breaks = seq(from = 0, to = xupper,by = xupper/10))+
    scale_y_continuous(breaks = seq(from = 0, to = yupper,by = yupper/10))+
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
                 size = 1)

  
    
    svg(filename = "Overarching curve_Hill.svg",
        width = 10,height = 10)
    # suppressWarnings(
    #   print(plotcurves)
    # )
    print(plotcurves)
    dev.off()
}

# show study-specific curves for Hill partial hierarchical model
showcurve_specific_Linear <- function(df_posteriors,DataList,
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
  doseinput <- dosed
  # doseinput <- xdose
  
  NStudy <- length(unique(df_input$Index))          # number of studies
  a <- df_posteriors[,aname]
  b <- df_posteriors[,bname]
  
  # # plotting method 1: dose ~ median of iterative RR by all parameters
  y_overarching <- matrix(NA, nrow = iterEff,ncol = length(xdose))

  # each iter of pars fits a distribution of y then take descriptives
  for(i in 1:iterEff){
    y_overarching[i,] <- get_Linear(doseinput,a[i],b[i])
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
  
  # # # Plotting method 2: dose ~ RR by median parameters
  # y_median <- get_Linear(doseinput,a = median(a), b = median(b))
  # y_p5 <- get_Linear(x = doseinput, a = quantile(a,0.05), b = quantile(b, 0.05))
  # y_p95 <- get_Linear(x = doseinput, a = quantile(a,0.95), b = quantile(b, 0.95))
  # y_p25 <- get_Linear(x = doseinput, a = quantile(a,0.25), b = quantile(b, 0.25))
  # y_p75 <- get_Linear(x = doseinput, a = quantile(a,0.75), b = quantile(b, 0.75))
  # df_plot <- data.frame(x = xdose,
  #                       median = y_median,Q5 = y_p5,Q95 = y_p95,
  #                       Q25 = y_p25,Q75 = y_p75)
  
  
  # Plotting
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
  
  
  svg(filename = "Overarching curve_Linear.svg",
      width = 10,height = 10)
  suppressWarnings(
    print(plotcurves)
  )
  dev.off()
}

printEvalResults <- function(dataName,ls_results){
  require("xlsx")
  ls_results <- eval(parse(text = paste0("Results_",dataName)))
  
  df_data <- data.frame(
    Index = ls_results$RawData$Index,
    Study = ls_results$RawData$Study,
    Nsub = ls_results$RawData$Nsub,
    dose = ls_results$RawData$dose,
    ymean = ls_results$RawData$ymean,
    ysd = ls_results$RawData$ysd
  )
  
  filename <- paste0(dataName," model evaluation results.xlsx")
  write.xlsx(x = df_data,
             file = filename,
             sheetName = "data",
             row.names = F)
  write.xlsx(x = data.frame(summary(ls_results[[2]]$Linear)[1]),
             file = filename,
             sheetName = "stanfit_Linear",
             append = T)
  write.xlsx(x = data.frame(summary(ls_results[[2]]$Power)[1]),
             file = filename,
             sheetName = "stanfit_Power",
             append = T)
  write.xlsx(x = data.frame(summary(ls_results[[2]]$Hill)[1]),
             file = filename,
             sheetName = "stanfit_Hill",
             append = T)
  write.xlsx(x = data.frame(summary(ls_results[[2]]$Expo5)[1]),
             file = filename,
             sheetName = "stanfit_Expo5",
             append = T)
  
  write.xlsx(x = ls_results[[3]],
             file = filename,
             sheetName = "Weights by approaches",
             append = T)
  
  write.xlsx(x = ls_results$WeightedMean,
             file = filename,
             sheetName = "mean of weighted",
             append = T
  )
  
  write.xlsx(x = ls_results$ModelMetrics,
             file = filename,
             sheetName = "Weighting metrics",
             append = T)
  # addDataFrame(
  #   x = as.data.frame(ls_results$ModelChoice),
  #   sheet = "Weighting metrics",
  #   startRow = 5
  # )
  
}

ShowWeightedValues <- function(DataList,modelwts,weighted,
                               dataName = "data",
                               dev = "png",
                               saveunit = "in",
                               size = 15){
  # # Data used for testing and debugging
  # dose <- c(1,2,3,1,2,3,1,2,3)
  # Nstudy <- 3
  # doseZ <- dnormalize_minmax(dose)
  # Index <- rep(1:3,each = 3)
  # ymean <- c(1,1.1,1.2,1,1.2,1.3,1,1.5,1.6)
  # ysd <- c(0,0.2,0.3,0,0.1,0.15,0,0.2,0.1)
  # Nmethods <- 2
  # methodlist <- c("A","B")
  # Nmodel <- 4
  # modellist <- c("alpha","beta","gamma","eta")
  # modelwts <- data.frame(
  #   alpha = c(0.1,0.5),
  #   beta = c(0.2,0.1),
  #   gamma = c(0.3,0.1),
  #   eta = c(0.4,0.3)
  # )
  # rownames(modelwts) <- methodlist
  # Niter <- 10
  # Nobs <- 9
  # weighted <- list(
  #   data.frame(replicate(Nobs,runif(Niter,1,1.5))),
  #   data.frame(replicate(Nobs,runif(Niter,1,1.5)))
  # )
  # y_weighted <- unlist(lapply(weighted,function(x){
  #   gather(as.data.frame(t(x)),"y")$value
  # }))
  # df_wts <- cbind(
  #   gather(modelwts,"model","weight") ,
  #   Approach = rep(methodlist,times = Nmethods),
  #   Row = "1"
  # ) 
  
  require("ggplot2")
  require("dplyr")
  require("tidyr")
  require("ggpubr")
  
  dose <- DataList$dose
  doseZ <- DataList$doseZ
  Index <- DataList$Index
  ymean <- DataList$ymean
  ysd <- DataList$ysd
  Nstudy <- DataList$S
  y_weighted <- unlist(lapply(weighted,function(x){
    gather(as.data.frame(t(x)),"y")$value
  }))
  
  modellist <- rownames(modelwts)
  Nmodel <- nrow(modelwts)
  methodlist <- colnames(modelwts)
  Nmethods <- ncol(modelwts)
  Niter <- dim(weighted[[1]])[1]
  Nobs <- dim(weighted[[2]])[2]
  
  
  df_weighted <- data.frame(
    Approach = rep(methodlist,each = Niter * Nobs),
    iter = rep(rep(1:Niter,each = Nobs),times = Nmethods),
    dose = rep(dose,times = Nmethods * Niter),
    doseZ = rep(doseZ,times = Nmethods * Niter),
    Index = rep(Index,times = Nmethods * Niter),
    y = y_weighted,
    ymean = rep(ymean,times = Nmethods * Niter) ,
    ysd = rep(ysd, times = Nmethods * Niter)
  )
  # df_wts <- cbind(
  #   gather(as.data.frame(modelwts),"Approach","weight") ,
  #   model = rep(modellist,times = Nmodel),
  #   Row = "1"
  # ) 
  
  plot <- ggplot(data = df_weighted)+
    # ggplot(data = df_weighted, aes(x = dose))+
    geom_line(aes(x = dose,y = y,group = iter),col = "grey")+
    stat_summary(aes(x = dose,y = y),fun = "mean",geom="line",colour="red")+
    stat_summary(aes(x = dose,y = y), fun = function(x){quantile(x,c(0.05))},geom = "line",colour = "blue")+
    stat_summary(aes(x = dose,y = y), fun = function(x){quantile(x,c(0.95))},geom = "line",colour = "blue")+
    geom_point(aes(x = dose,y = ymean),col="black")+
    geom_errorbar(aes(x = dose,ymin = ymean - ysd, ymax = ymean + ysd))+
    facet_wrap(Index ~ Approach,scales = "free",
               nrow = max(Index))+    
    scale_x_continuous(name = "Dose")+
    scale_y_continuous(name = "Response")+
    theme_classic(base_size = 15)

  if(saveunit %in% c("cm","in")){
    savescale <- 1
  } else if(savescale == "mm"){
    savescale <- 10
  } else if(saveunit == "px"){
    savescale <- 100
  }
  
  filename <- paste0(dataName," Weighted Quantities across Averaging Approaches")
  # save main plot
  ggsave(filename = paste0(filename,".",dev),
         plot = plot,
         device = dev,
         width = size * savescale,
         height = size* savescale)
  # ggsave(filename = "Weights.png",
  #        plot = plot1,
  #        device = "png",
  #        width = 12,height = 12)
  # 
}
