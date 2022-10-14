# Generate example data
# Used in testing for comparison between BBMD and my script
# All study indexed as in meta format, i.e. including study index

library("dplyr")
library("xlsx")
iter <- 1E3
Nchain <- 4
warmup <- 0.5
source("BHBMD model utilities v2.r")

# Normal----

# hypothetical study info
# PFOA relative kidney weight change male rats
# dose: 0,10,50,250
# N   : 50,50,50,50

# Normal individual
df_norm_individual <- data.frame(
  N = c(50,50,50,50),

  Index = rep(1,sum(c(50,50,50,50))),
  dose = c(
    rep(0,50),
    rep(10,50),
    rep(50,50),
    rep(250,50)
  ),
  Y = c(
    rnorm(50,2.5,0.2),
    rnorm(50,2.6,0.3),
    rnorm(50,4.7,0.2),
    rnorm(50,10.1,0.4)
  )
)

# Normal Summary
df_norm_summary <- data.frame(
  N = c(50,50,50,50),
  Index = rep(1,4),
  Y = group_by(df_norm_individual,dose) %>% summarise(mean = mean(Y), sd = sd(Y))
) %>% rename(
  dose = Y.dose,
  ymean = Y.mean,
  ysd = Y.sd
)

# Lognormal----
# Hypothetical study information
# PFHS developmental rats
# dose: 0,10,50,250
# N   : 50,50,50,50

# Lognormal individual
df_lnorm_individual <- data.frame(
  N = c(50,50,50,50),
  
  Index = rep(1,sum(c(50,50,50,50))),
  dose = c(
    rep(0,50),
    rep(10,50),
    rep(50,50),
    rep(250,50)
  ),
  Y = c(
    rlnorm(50,1.1,0.12),
    rlnorm(50,1.2,0.14),
    rlnorm(50,1.5,0.09),
    rlnorm(50,2.2,0.21)
  )
)

# lognormal summary
df_lnorm_summary <- data.frame(
  N = c(50,50,50,50),
  Index = rep(1,4),
  Y = group_by(df_lnorm_individual,dose) %>% summarise(mean = mean(Y), sd = sd(Y))
) %>% rename(
  dose = Y.dose,
  ymean = Y.mean,
  ysd = Y.sd
)

write.xlsx(df_norm_individual,file = "Example Dose-response Data.xlsx",
           sheetName = "Norm Individual",row.names = F)
write.xlsx(df_norm_summary,file = "Example Dose-response Data.xlsx",
           sheetName = "Norm Summary",row.names = F,append = T)
write.xlsx(df_lnorm_individual,file = "Example Dose-response Data.xlsx",
           sheetName = "Lnorm Individual",row.names = F,append = T)
write.xlsx(df_lnorm_summary,file = "Example Dose-response Data.xlsx",
           sheetName = "Lnorm Summary",row.names = F,append = T)

save(
  df_norm_individual,df_norm_summary,
  df_lnorm_individual,df_lnorm_summary,
  file = "Example Data.Rdata"
)

# Model based simulation----
S <- 10
G <- 3
sigmavalue <- 0.21
Nsub <- ceiling(replicate(S,c(250,210,180)+runif(G,0,60)))
dose <- replicate(S,sort(runif(G,0,100)))
y_bound <- sort(
  runif(2,min = exp(-10),max = 20)
)

N <- S*G
Index <- rep(1:S,each = G)
Study <- LETTERS[Index]
alldose <- replicate(S,seq(from = 0, to = max(dose), lengt = 1E3))
x_bound <- c(min(dose),max(dose))
noisesize <- 0.1

# fdL_upper <- runif(1,9,10)               # fdL is not log
# fdL_lower <- runif(1,0.8,1.2)
# slope_upper <- (fdL_upper - fdL_lower) / max(dose)


# Linear
# b_Linear <- (1 - 0.2 * runif(S,0,1)) * slope_upper
# a_Linear <- fdL_lower * (1 + runif(S,0,1)*0.5) 
popfit_linear <- optim(c(1,1),fit_Linear, input = x_bound, y = y_bound)
pars_linear <- popfit_linear$par
get_Linear(x_bound,pars_linear[1],pars_linear[2])
a_Linear <- pars_linear[1] + rnorm(S,0,noisesize * pars_linear[1])
b_Linear <- pars_linear[2] + rnorm(S,0,noisesize * pars_linear[2])

# Hill
popfit_Hill <- optim(c(1,1,median(dose),1),fit_Hill, input = x_bound, y = y_bound,
                   method = c("L-BFGS-B"),lower = c(-Inf,-Inf,0,0),upper = c(Inf,Inf,Inf,18))
pars_Hill <- popfit_Hill$par
get_Hill(x_bound,pars_Hill[1],pars_Hill[2],pars_Hill[3],pars_Hill[4])
a_Hill <- pars_Hill[1] + rnorm(S,0,noisesize * pars_Hill[1])
b_Hill <- pars_Hill[2] + rnorm(S,0,noisesize * pars_Hill[2])
c_Hill <- pars_Hill[3] + rnorm(S,0,noisesize * pars_Hill[3])
g_Hill <- pars_Hill[4] + rnorm(S,0,noisesize * pars_Hill[4])

# a_Hill <- fdL_lower * (1 + runif(S,0,1)*0.5) 
# b_Hill <- fdL_upper * (1- runif(S,0,1)*0.1) - fdL_lower
# c_Hill <-  median(dose) * ( 1 + runif(S,-1,1) * 0.4)
# g_Hill <- log(b_Hill,base = 2) + runif(S,-1,1) * 0.1

# Power
# a_Pow <- fdL_lower * (1 + runif(S,0,1)*0.5) 
# g_Pow <- runif(S,1,3)
# b_Pow <- (fdL_upper * runif(S,0.75,1.25) - a_Pow) / max(dose)^g_Pow
popfit_Pow <- optim(c(1,1,1),fit_Power, input = x_bound, y = y_bound,
                    method = "L-BFGS-B",lower = c(-Inf,-Inf,0), upper = c(Inf,Inf,18))
pars_Pow<- popfit_Pow$par
get_Power(x_bound,pars_Pow[1],pars_Pow[2],pars_Pow[3])
a_Pow <- pars_Pow[1] + rnorm(S,0,noisesize * pars_Pow[1])
b_Pow <- pars_Pow[2] + rnorm(S,0,noisesize * pars_Pow[2])
g_Pow <- pars_Pow[3] + rnorm(S,0,noisesize * pars_Pow[3])

# Michaelis Menten 
popfit_MM <- optim(c(1,1,1),fit_MM, input = x_bound, y = y_bound,
                    method = "L-BFGS-B",lower = c(0,-Inf,0), upper = c(Inf,Inf,18))
pars_MM<- popfit_MM$par
get_MM(x_bound,pars_MM[1],pars_MM[2],pars_MM[3])
a_MM <- pars_MM[1] + rnorm(S,0,noisesize * pars_MM[1])
b_MM <- pars_MM[2] + rnorm(S,0,noisesize * pars_MM[2])
c_MM <- pars_MM[3] + rnorm(S,0,noisesize * pars_MM[3])

# Exponential 5
popfit_Expo5 <- optim(c(1,1,1,1),fit_Expo5, input = x_bound, y = y_bound,
                      method = "L-BFGS-B",lower = c(0,0,1,1), upper = c(Inf,Inf,Inf,18))
pars_Expo5 <- popfit_Expo5$par
get_Expo5(x_bound,pars_Expo5[1],pars_Expo5[2],pars_Expo5[3],pars_Expo5[4])
a_Expo5 <- pars_Expo5[1] + rnorm(S,0,noisesize * pars_Expo5[1])
b_Expo5 <- pars_Expo5[2] + rnorm(S,0,noisesize*0.5 * pars_Expo5[2])
c_Expo5 <- pars_Expo5[3] + rnorm(S,0,noisesize * pars_Expo5[3])
d_Expo5 <- pars_Expo5[4] + rnorm(S,0,noisesize * pars_Expo5[4])

# Plotting and simulation fd
fd_Linear <- vector("list",S)
fd_Hill <- vector("list",S)
fd_MM <- vector("list",S)
fd_Power <- vector("list",S)
fd_Expo5 <- vector("list",S)

xdose <- dose           # use dose for simulation/ alldose for plotting

plot(0,xlim = c(0,max(dose)), ylim = c(0,y_bound[2]))
for(s in 1:S){
  fd_Linear[[s]] <- get_Linear(xdose[,s],a = a_Linear[s], b = b_Linear[s])
  lines(xdose[,s],fd_Linear[[s]],col = "blue")
  fd_Hill[[s]] <- get_Hill(xdose[,s],a = a_Hill[s], b = b_Hill[s], c = c_Hill[s], g = g_Hill[s])
  lines(xdose[,s],fd_Hill[[s]], col = "red")
  fd_Power[[s]] <- get_Power(xdose[,s], a = a_Pow[s], b = b_Pow[s], g = g_Pow[s])
  lines(xdose[,s],fd_Power[[s]], col = "green")
  fd_MM[[s]] <- get_MM(xdose[,s], a = a_MM[s], b = b_MM[s], c = c_MM[s])
  lines(xdose[,s],fd_MM[[s]], col = "purple")
  fd_Expo5[[s]] <- get_Expo5(xdose[,s], a = a_Expo5[s], b = b_Expo5[s], c = c_Expo5[s], d = d_Expo5[s])
  lines(xdose[,s],fd_Expo5[[s]],col = "yellow")
}

y_Linear <- vector("list",N)
y_Hill <- vector("list",N)
y_Power <- vector("list",N)
y_MM <- vector("list",N)
y_Expo5 <- vector("list",N)
for(g in 1:N){
  y_Linear[[g]] <- rlnorm(Nsub[g],meanlog = log(fd_Linear[[Index[g]]][g-3*(Index[g]-1)]),sdlog = sigmavalue)
  y_Hill[[g]] <- rlnorm(Nsub[g],meanlog = log(fd_Hill[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
  y_MM[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_MM[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
  y_Power[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_Power[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
  y_Expo5[[g]] <- rlnorm(Nsub[g], meanlog = log(fd_Expo5[[Index[g]]][g-3*(Index[g]-1)]), sdlog = sigmavalue)
}
y_Linear_Summary <- lapply(y_Linear,function(x){
  cbind(mean(x),sd(x))
})
y_Hill_Summary <- lapply(y_Hill,function(x){
  cbind(mean(x),sd(x))
})
y_MM_Summary <- lapply(y_MM,function(x){
  cbind(mean(x),sd(x))
})
y_Power_Summary <- lapply(y_Power,function(x){
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
  ymean = unlist(lapply(y_Linear_Summary,function(x) x[,1])),
  ysd = unlist(lapply(y_Linear_Summary,function(x) x[,2]))
)

df_simu_meta_Hill <- data.frame(
  Index = Index,
  Study = Study,
  dose = as.vector(dose),
  Nsub = as.vector(Nsub),
  ymean = unlist(lapply(y_Hill_Summary,function(x) x[,1])),
  ysd = unlist(lapply(y_Hill_Summary,function(x) x[,2]))
)

df_simu_meta_MM <- data.frame(
  Index = Index,
  Study = Study,
  dose = as.vector(dose),
  Nsub = as.vector(Nsub),
  ymean = unlist(lapply(y_MM_Summary,function(x) x[,1])),
  ysd = unlist(lapply(y_MM_Summary,function(x) x[,2]))
)

df_simu_meta_Power <- data.frame(
  Index = Index,
  Study = Study,
  dose = as.vector(dose),
  Nsub = as.vector(Nsub),
  ymean = unlist(lapply(y_Power_Summary,function(x) x[,1])),
  ysd = unlist(lapply(y_Power_Summary,function(x) x[,2]))
)

df_simu_meta_Expo5 <- data.frame(
  Index = Index,
  Study = Study,
  dose = as.vector(dose),
  Nsub = as.vector(Nsub),
  ymean = unlist(lapply(y_Expo5_Summary,function(x) x[,1])),
  ysd = unlist(lapply(y_Expo5_Summary,function(x) x[,2]))
)

# Export data
baseinfo <- list(
  x_bound,y_bound,
  a_Linear,b_Linear,
  a_Pow,b_Pow,g_Pow,
  a_MM,b_MM,c_MM,
  a_Hill,b_Hill,c_Hill,g_Hill,
  a_Expo5,b_Expo5,c_Expo5,d_Expo5
)
simulated_data <- list(
  Linear = df_simu_meta_Linear,
  Power = df_simu_meta_Power,
  MM = df_simu_meta_MM,
  Hill = df_simu_meta_Hill,
  Expo5 = df_simu_meta_Expo5,
  baseinfo = baseinfo
)
save(
  simulated_data,
  file = "Simulated Data for Model Eval.rdata"
)
