# Settings-------
# mywd <- "/N/scratch/zhezhou/MAS_Expo5/"
# fileDir <- "/N/slate/zhezhou/"

mywd <- paste0(getwd(),"/")
fileDir <- mywd
DSName <- "bladder_meta"

# seed <- 47405      # use IUB postal code as seed for random simulation
iter_ppd <- 1e3 *2 
samplingiter <- 1E4*5   # number of iterations in stan sampling
warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 3       # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling
options(bitmapType='cairo')
source(file = paste0(fileDir,"BHBMD_model_utilities.R"),
       local = T)

# Local functions--------
library(EnvStats)
library(extraDistr)

PPD_scatter <- function(modelname,yarm){
  ypred <- eval(parse(text = paste0("y_",modelname))) 
  df_plot <- melt(ypred) %>% rename(iter = Var1, Ypred = value) %>% 
    arrange(iter) %>% mutate(
      Xvalue = rep(Xpred,times=iter_ppd)  )
  Xselect <- c(0,0.2,0.4,0.6,0.8,1)
  df_plot2 <- df_plot %>% filter(Xvalue %in% Xselect) %>% arrange(iter,Xvalue)
  
  # scatter plot
  plot <- ggplot(data = df_plot) +
    geom_point(aes(x = Xvalue,y = log(Ypred),group = iter),col = "light blue")+
    # add reference lines
    # geom_hline(yintercept = log(500),col="green")+
    # geom_hline(yintercept = log(0.1),col="green")+
    # annotate(geom = "text",
    #          label = "log(RR=0.1)",
    #          x= 1, y = -2,size = 8, col = "red")+
    # annotate(geom = "text",
    #          label = "log(RR=500)",
    #          x = 1, y = 7,size = 8, col = "red")
    labs(title = paste0("Simulations from Prior Predictive Distribution of ",modelname," DRmodel"))+
    # trim coordinates
    # coord_cartesian(ylim = c(-8,10))+
    scale_y_continuous(
      # breaks = c(
      #   seq(from = round(quantile(log(ypred),0.0001),0),
      #       to = round(quantile(log(ypred),0.9999),0),
      #       by = 1)),
      limits = c(-yarm,yarm),
      name = "Logarithm of Relative Risks"
      
    )+
    scale_x_continuous(name = "Standardized Dose [0,1]",
                       limits = c(0,1),
                       breaks = Xselect)+
    geom_vline(xintercept = Xselect)+
    # improve format
    theme_classic(base_size = 25)
  
  png(paste0("Scatterplot_",modelname,".png"),
      width =  1600, height =  1200)
  print(plot)
  dev.off()
  
}
PPD_kernel <- function(modelname,yarm){
  ypred <- eval(parse(text = paste0("y_",modelname))) 
  df_plot <- melt(ypred) %>% rename(iter = Var1, Ypred = value) %>% 
    arrange(iter) %>% mutate(
      Xvalue = rep(Xpred,times=iter_ppd)  )
  Xselect <- c(0,0.2,0.4,0.6,0.8,1)
  df_plot2 <- df_plot %>% filter(Xvalue %in% Xselect) %>% arrange(iter,Xvalue)
  
  # kernel density at sparse data points
  plot2 <- ggplot(df_plot2) + 
    # kernel density
    geom_density(mapping = aes(y = log(Ypred)),col="red",size = 3)+
    facet_wrap(~Xvalue,ncol = length(Xselect),scales = "free_x")+
    scale_y_continuous(limits = c(-yarm,yarm))+
    theme_classic(base_size = 25)
  
  png(paste0("Kernel_density_",modelname,".png"),
      width =  1600, height =  1200)
  print(plot2)
  dev.off()
  
}
# Load data----
df_DataList <- read.csv(file = paste0("DataList_",DSName,".csv"))

# Get prior range------
DataList <- PrepareLnormData(df_DataList,RR=T)
Xrange <- sort(c(max(DataList$doseZ),min(DataList$doseZ)))

Xpred <- seq(from = Xrange[1] , to= Xrange[2],by = 0.001)
# Xpred <- DataList$dose
Nobs <- length(Xpred)

ybound <- c(
  max(DataList$ymean),
  sd(DataList$ymean)
)
xbound <- c(
  Xrange[1],
  Xrange[2]
)
a_bound <- log(ybound[1] + 3 * ybound[2])*2

# Hyperparameters Linear-----

prior_Lognormal_Summary_Linear_meta_3 <- list(
  # prior_sigma = c(0,2),
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  # prior_sigma_a = c(0,2),
  prior_sigma_a = c(1),
  prior_scale_b =c(1)
  # prior_a_b = c(1,10),
  # prior_b_b = c(1,10)
)

# Hyperparameters Power-----
prior_Lognormal_Summary_Power_meta_4 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  # prior_a_b = c(1,10),
  # prior_b_b = c(1,10),
  prior_scale_b =c(1),
  prior_mu_g = c(-1,1),
  prior_sigma_g = c(1),
  g_constr = c(1,18)
)

# Hyperparameters MM-----
prior_Lognormal_Summary_MM_meta_4 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  # prior_mu_b = c(-2,2),
  # prior_sigma_b = c(0,2),
  # prior_a_b = c(1,10),
  # prior_b_b = c(1,10),
  prior_scale_b =c(1),
  prior_mu_c = c(-2,2),
  prior_sigma_c = c(1)
)

# Hyperparameters Hill-----
prior_Lognormal_Summary_Hill_meta_5 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  # prior_mu_b = c(-2,2),
  # prior_sigma_b = c(0,1.5),
  # prior_a_b = c(1,10),
  # prior_b_b = c(1,10),
  prior_scale_b = c(1),
  prior_mu_c = c(-2,2),
  prior_sigma_c = c(1),
  prior_mu_g = c(-1,1),
  prior_sigma_g = c(1),
  g_constr = c(0,18)
)

# Hyperparameters Expo5-----
prior_Lognormal_Summary_Expo5_meta_5 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  # prior_a_b = c(1,10),
  # prior_b_b = c(1,10),
  prior_scale_b =c(1),
  # prior_mu_c = c(-2,2),
  # prior_sigma_c = c(0,2),
  prior_scale_c =c(1),
  c_lower = c(1),
  # prior_mu_c = c(-2,2),
  # prior_sigma_c = c(0,1.5),
  prior_mu_d = c(-1,1),
  prior_sigma_d = c(1),
  d_constr = c(0,18)
)

# Test-----
x_test <- seq(0,5,by=0.01)
y_test <- get_Power(x_test,a=1,b=1,g = 3)  # intercept slope power
y_test <- get_MM(x_test,a=1,b=4,c=0.7)   # intercept upperbound half
y_test <- get_Hill(x_test,a=0,b=4,     # intercept upperbound
                   c=0.5,g=0.5)        # half power
y_test <- get_Expo5(x_test,
                    a=1,b=2,          # intercept slope
                    c=3,d=2)          # upperbound  power
plot(x_test,y_test,ylim=c(0,max(y_test)))

# PPD Linear-------
sigma <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Linear_meta_3$prior_sigma[1],
                 prior_Lognormal_Summary_Linear_meta_3$prior_sigma[1]
              )
mu_a <- runif(iter_ppd,
              prior_Lognormal_Summary_Linear_meta_3$prior_mu_a[1],
              prior_Lognormal_Summary_Linear_meta_3$prior_mu_a[2])
sigma_a <- rexp(iter_ppd,
                   # prior_Lognormal_Summary_Linear_meta_3$prior_sigma_a[1],
                   prior_Lognormal_Summary_Linear_meta_3$prior_sigma_a[1])
# a_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Linear_meta_3$prior_a_b[1],
#               prior_Lognormal_Summary_Linear_meta_3$prior_a_b[2])
# b_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Linear_meta_3$prior_b_b[1],
#               prior_Lognormal_Summary_Linear_meta_3$prior_b_b[2])
scale_b <- rexp(iter_ppd,
                prior_Lognormal_Summary_Linear_meta_3$prior_scale_b[1])

a <- rlnorm(iter_ppd,meanlog = mu_a,sdlog = sigma_a)
# b <- rlnorm(iter_ppd,meanlog = mu_b,sdlog = sigma_b)
# b <- rinvgamma(iter_ppd,a_b,b_b)
b <- rhcauchy(iter_ppd,scale_b)
fd_Linear <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
y_Linear <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
for(i in 1:iter_ppd){
  fd_Linear[i,] <- get_Linear(Xpred,a[i],b[i])
  y_Linear[i,] <- rlnorm(Nobs,log(fd_Linear[i,]),sigma[i])
}

# PPD Power-------
sigma <- rexp(iter_ppd,
               # prior_Lognormal_Summary_Power_meta_4$prior_sigma[1],
               prior_Lognormal_Summary_Power_meta_4$prior_sigma[1]
               )
mu_a <- runif(iter_ppd,
              prior_Lognormal_Summary_Power_meta_4$prior_mu_a[1],
              prior_Lognormal_Summary_Power_meta_4$prior_mu_a[2])
sigma_a <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Power_meta_4$prior_sigma_a[1],
                 prior_Lognormal_Summary_Power_meta_4$prior_sigma_a[1])
# mu_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Power_meta_4$prior_mu_b[1],
#               prior_Lognormal_Summary_Power_meta_4$prior_mu_b[2])
# sigma_b <- runif(iter_ppd,
#                  prior_Lognormal_Summary_Power_meta_4$prior_sigma_b[1],
#                  prior_Lognormal_Summary_Power_meta_4$prior_sigma_b[2])
# a_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Power_meta_4$prior_a_b[1],
#               prior_Lognormal_Summary_Power_meta_4$prior_a_b[2])
# b_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Power_meta_4$prior_b_b[1],
#               prior_Lognormal_Summary_Power_meta_4$prior_b_b[2])
scale_b <- rexp(iter_ppd,
                prior_Lognormal_Summary_Power_meta_4$prior_scale_b[1])
mu_g <- runif(iter_ppd,
              prior_Lognormal_Summary_Power_meta_4$prior_mu_g[1],
              prior_Lognormal_Summary_Power_meta_4$prior_mu_g[2])
sigma_g <- rexp(iter_ppd,
              # prior_Lognormal_Summary_Power_meta_4$prior_sigma_g[1],
              prior_Lognormal_Summary_Power_meta_4$prior_sigma_g[1])
# s_g <- runif(iter_ppd,
#              prior_Lognormal_Summary_Power_meta_4$prior_s_g[1],
#              prior_Lognormal_Summary_Power_meta_4$prior_s_g[2])
# r_g <- runif(iter_ppd,
#              prior_Lognormal_Summary_Power_meta_4$prior_r_g[1],
#              prior_Lognormal_Summary_Power_meta_4$prior_r_g[2])
 
a <- rlnorm(iter_ppd,mu_a,sigma_a)
# b <- rlnorm(iter_ppd,mu_b,sigma_b)
# b <- rinvgamma(iter_ppd,a_b,b_b)
b <- rhcauchy(iter_ppd,scale_b)
# g <- rtgamma(iter_ppd,shape = s_g,scale = 1/r_g,
#              min = prior_Lognormal_Summary_Power_meta_4$g_constr[1],
#              max = prior_Lognormal_Summary_Power_meta_4$g_constr[2])
# g <- EnvStats::rlnormTrunc(iter_ppd,
#                            meanlog = mu_g,sdlog=sigma_g,
#                            min = prior_Lognormal_Summary_Power_meta_4$g_constr[1],
#                            max = prior_Lognormal_Summary_Power_meta_4$g_constr[2])

g <- vector(length = iter_ppd)
fd_Power <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
y_Power <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
for(i in 1:iter_ppd){
  g[i] <- rlnorm(1,mu_g[i],sigma_g[i])
  # g[i] <- rgamma(1,s_g[i],r_g[i])
  # g[i] <- rhcauchy(1,sigma_g)
  while((g[i] < prior_Lognormal_Summary_Power_meta_4$g_constr[1]) |
        (g[i] > prior_Lognormal_Summary_Power_meta_4$g_constr[2])){
    g[i] <- 1
  }
}
for(i in 1:iter_ppd){
  fd_Power[i,] <- get_Power(Xpred,a[i],b[i],g[i])
  y_Power[i,] <- rlnorm(Nobs,log(fd_Power[i,]),sigma[i])
}

# PPD MM------
sigma <- rexp(iter_ppd,
               prior_Lognormal_Summary_MM_meta_4$prior_sigma[1]
               # prior_Lognormal_Summary_MM_meta_4$prior_sigma[2]
              )
mu_a <- runif(iter_ppd,
              prior_Lognormal_Summary_MM_meta_4$prior_mu_a[1],
              prior_Lognormal_Summary_MM_meta_4$prior_mu_a[2])
sigma_a <- rexp(iter_ppd,
                 prior_Lognormal_Summary_MM_meta_4$prior_sigma_a[1]
                 # prior_Lognormal_Summary_MM_meta_4$prior_sigma_a[2]
                )
# mu_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_MM_meta_4$prior_mu_b[1],
#               prior_Lognormal_Summary_MM_meta_4$prior_mu_b[2])
# sigma_b <- runif(iter_ppd,
#                  prior_Lognormal_Summary_MM_meta_4$prior_sigma_b[1],
#                  prior_Lognormal_Summary_MM_meta_4$prior_sigma_b[2])
# a_b <- runif(iter_ppd,
#              prior_Lognormal_Summary_MM_meta_4$prior_a_b[1],
#              prior_Lognormal_Summary_MM_meta_4$prior_a_b[2])
# b_b <- runif(iter_ppd,
#              prior_Lognormal_Summary_MM_meta_4$prior_b_b[1],
#              prior_Lognormal_Summary_MM_meta_4$prior_b_b[2])
scale_b <- rexp(iter_ppd,
                prior_Lognormal_Summary_MM_meta_4$prior_scale_b[1])
mu_c <- runif(iter_ppd,
              prior_Lognormal_Summary_MM_meta_4$prior_mu_c[1],
              prior_Lognormal_Summary_MM_meta_4$prior_mu_c[2])
sigma_c <- rexp(iter_ppd,
                 prior_Lognormal_Summary_MM_meta_4$prior_sigma_c[1]
                 # prior_Lognormal_Summary_MM_meta_4$prior_sigma_c[2]
                )

a <- rlnorm(iter_ppd,mu_a,sigma_a)
# b <- rlnorm(iter_ppd,mu_b,sigma_b)
# b <- rinvgamma(iter_ppd,a_b,b_b)
b <- rhcauchy(iter_ppd,scale_b)
c_par <- rlnorm(iter_ppd,mu_c,sigma_c)

fd_MM <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
y_MM <- matrix(NA,nrow = iter_ppd,ncol = Nobs)

for(i in 1:iter_ppd){
  fd_MM[i,] <- get_MM(Xpred,a[i],b[i],c_par[i])
  y_MM[i,] <- rlnorm(Nobs,log(fd_MM[i,]),sigma[i])
}

# PPD Hill----
sigma <- rexp(iter_ppd,
               # prior_Lognormal_Summary_Hill_meta_5$prior_sigma[1],
               prior_Lognormal_Summary_Hill_meta_5$prior_sigma[1])
mu_a <- runif(iter_ppd,
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_a[1],
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_a[2])
sigma_a <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Hill_meta_5$prior_sigma_a[1],
                 prior_Lognormal_Summary_Hill_meta_5$prior_sigma_a[1])
# a_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Hill_meta_5$prior_a_b[1],
#               prior_Lognormal_Summary_Hill_meta_5$prior_a_b[2])
# b_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Hill_meta_5$prior_b_b[1],
#               prior_Lognormal_Summary_Hill_meta_5$prior_b_b[2])
scale_b <- rexp(iter_ppd,
                prior_Lognormal_Summary_Hill_meta_5$prior_scale_b[1])
mu_c <- runif(iter_ppd,
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_c[1],
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_c[2])
sigma_c <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Hill_meta_5$prior_sigma_c[1],
                 prior_Lognormal_Summary_Hill_meta_5$prior_sigma_c[1])
mu_g <- runif(iter_ppd,
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_g[1],
              prior_Lognormal_Summary_Hill_meta_5$prior_mu_g[2])
sigma_g <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Hill_meta_5$prior_sigma_g[1],
                 prior_Lognormal_Summary_Hill_meta_5$prior_sigma_g[1])

a <- rlnorm(iter_ppd,mu_a,sigma_a)
# b <- rinvgamma(iter_ppd,a_b,b_b)
b <- rhcauchy(iter_ppd,scale_b)
c_par <- rlnorm(iter_ppd,mu_c,sigma_c)
g <- vector(length = iter_ppd)
for(i in 1:iter_ppd){
  g[i] <- rlnorm(1,mu_g[i],sigma_g[i])
  # g[i] <- rgamma(1,s_g[i],r_g[i])
  # g[i] <- rhcauchy(1,sigma_g)
  while((g[i] < prior_Lognormal_Summary_Hill_meta_5$g_constr[1]) |
        (g[i] > prior_Lognormal_Summary_Hill_meta_5$g_constr[2])){
    g[i] <- 1
  }
}
fd_Hill <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
y_Hill <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
for(i in 1:iter_ppd){
  fd_Hill[i,] <- get_Hill(Xpred,a[i],b[i],c_par[i],g[i])
  y_Hill[i,] <- rlnorm(Nobs,log(fd_Hill[i,]),sigma[i])
}

# PPD Expo5----
sigma <- rexp(iter_ppd,
               # prior_Lognormal_Summary_Expo5_meta_5$prior_sigma[1],
               prior_Lognormal_Summary_Expo5_meta_5$prior_sigma[1])
mu_a <- runif(iter_ppd,
              prior_Lognormal_Summary_Expo5_meta_5$prior_mu_a[1],
              prior_Lognormal_Summary_Expo5_meta_5$prior_mu_a[2])
sigma_a <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_a[1],
                 prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_a[1])
# a_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Expo5_meta_5$prior_a_b[1],
#               prior_Lognormal_Summary_Expo5_meta_5$prior_a_b[2])
# b_b <- runif(iter_ppd,
#               prior_Lognormal_Summary_Expo5_meta_5$prior_b_b[1],
#               prior_Lognormal_Summary_Expo5_meta_5$prior_b_b[2])
scale_b <- rexp(iter_ppd,
                prior_Lognormal_Summary_Expo5_meta_5$prior_scale_b[1])
# mu_c <- runif(iter_ppd,
#                      prior_Lognormal_Summary_Expo5_meta_5$prior_mu_c[1],
#                      prior_Lognormal_Summary_Expo5_meta_5$prior_mu_c[2])
# sigma_c <- rexp(iter_ppd,
#                      # prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_c[1],
#                      prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_c[1])
scale_c <- rexp(iter_ppd,
                prior_Lognormal_Summary_Expo5_meta_5$prior_scale_c[1])
mu_d <- runif(iter_ppd,
              prior_Lognormal_Summary_Expo5_meta_5$prior_mu_d[1],
              prior_Lognormal_Summary_Expo5_meta_5$prior_mu_d[2])
sigma_d <- rexp(iter_ppd,
                 # prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_d[1],
                 prior_Lognormal_Summary_Expo5_meta_5$prior_sigma_d[1])
a <- rlnorm(iter_ppd,mu_a,sigma_a)
b <- rhcauchy(iter_ppd,scale_b)
c_par <- vector(length = iter_ppd)
d <- vector(length = iter_ppd)
for(i in 1:iter_ppd){
  d[i] <- rlnorm(1,mu_d[i],sigma_d[i])
  while(d[i] > prior_Lognormal_Summary_Expo5_meta_5$d_constr[2]){
    d[i] <- d[i]/2
  }
  
  c_par[i] <- rhcauchy(1,scale_c)
  while(c_par[i] < prior_Lognormal_Summary_Expo5_meta_5$c_lower[1]){
    c_par[i] <- c_par[i]+0.5
  }
}

fd_Expo5 <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
y_Expo5 <- matrix(NA,nrow = iter_ppd,ncol = Nobs)
for(i in 1:iter_ppd){
  fd_Expo5[i,] <- get_Expo5(Xpred,a[i],b[i],c_par[i],d[i])
  y_Expo5[i,] <- rlnorm(Nobs,log(fd_Expo5[i,]),sigma[i])
}
# Visual check on PPD------
modelname <- "Hill"
PPD_kernel(modelname,yarm = 20)
PPD_scatter(modelname,yarm = 20)

