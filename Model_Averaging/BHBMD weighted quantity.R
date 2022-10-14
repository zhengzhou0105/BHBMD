# Settings-----
iter <- 1E4   # number of iterations in stan sampling
warmup <- 0.5      # proportion of steps used for warmup in stan sampling
Nchain <- 4       # number of chains in stan sampling
thin <- 1          # factor of thin in stan sampling
source("BHBMD model utilities v2.r")

# Input Data-----
fitList <- fitList_simu_Linear
wts_1 <- ModComp_Linear
DataList <- DataList

# Model averaging---------
input_quantity <- DataList$dose
y_original <- DataList$ymean

Method_list <- c("AIC","BIC","WAIC","LPML","Stacking","PseudoBMA","PseudoBMA_BB")

stanfit_Linear <- fitList$Linear
stanfit_Hill <- fitList$Hill
stanfit_Pow <- fitList$Power
stanfit_Expo5 <- fitList$Expo5

Nstan <- nrow(wts_1)
Nmethods <- ncol(wts_1)

fd_Linear <- extract(stanfit_Linear)$fd
fd_Hill <- extract(stanfit_Hill)$fd
fd_Pow <- extract(stanfit_Pow)$fd
fd_Expo5 <- extract(stanfit_Expo5)$fd


# df_posterior_Linear_0 <- as.data.frame(stanfit_Lognormal_Summary_Linear_meta_0_meta) %>% select(
#   a = a, b = b, sigma = sigma
# )
# 
# df_posterior_Linear_1 <- as.data.frame(stanfit_Lognormal_Summary_Linear_meta_1_meta) %>% mutate(
#   a = rlnorm(iterEff,meanlog = mu_a, sdlog = sigma_a)
# ) %>% select(  a = a, b = b, sigma = sigma)

# df_posterior_Linear <- as.data.frame(stanfit_Linear) %>% mutate(
#   a = rlnorm(iterEff,meanlog = mu_a, sdlog = sigma_a),
#   b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b)
# ) %>% select(a = a, b = b, sigma = sigma,fd = fd)
# 
# df_posterior_Hill <- as.data.frame(stanfit_Hill) %>% mutate(
#   a = rlnorm(iterEff,meanlog = mu_a, sdlog = sigma_a),
#   b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b),
#   c = rlnorm(iterEff,meanlog = mu_c, sdlog = sigma_c),
#   g = rlnorm(iterEff, meanlog = mu_g, sdlog = sigma_g)
# ) %>% select(a = a, b = b, sigma = sigma, c = c , g = g)
# 
# df_posterior_Pow <- as.data.frame(stanfit_Pow) %>% mutate(
#   a = rlnorm(iterEff,meanlog = mu_a, sdlog = sigma_a),
#   b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b),
#   g = rlnorm(iterEff, meanlog = mu_g, sdlog = sigma_g)
# ) %>% select(a = a, b = b, sigma = sigma, g = g) 
# 
# df_posterior_Expo5 <- as.data.frame(stanfit_Expo5) %>% mutate(
#   a = rlnorm(iterEff,meanlog = mu_a, sdlog = sigma_a),
#   b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b),
#   cminusone = rlnorm(iterEff, meanlog = mu_cminusone, sdlog = sigma_cminusone),
#   c = cminusone + 1,
#   d = rlnorm(iterEff, meanlog = mu_d, sdlog = sigma_d)
# ) %>% select(a = a, b = b, sigma = sigma, c = c , d = d)
  
# fd_Linear_0 <- matrix(NA,nrow = iterEff,ncol = length(input_quantity))
# fd_Linear_1 <- y_Linear_0
# fd_Linear <- matrix(NA,nrow = iterEff,ncol = length(input_quantity))
# fd_Hill <- matrix(NA,nrow = iterEff,ncol = length(input_quantity))
# fd_Pow <- matrix(NA,nrow = iterEff,ncol = length(input_quantity))
# fd_Expo5 <- matrix(NA,nrow = iterEff,ncol = length(input_quantity))

# for(i in 1:iterEff){
#   # fd_Linear_0[i,] <- get_Linear(x = input_quantity,a = df_posterior_Linear_0$a[i],
#   #                              b = df_posterior_Linear_0$b[i])
#   # fd_Linear_1[i,] <- get_Linear(x = input_quantity,a = df_posterior_Linear_1$a[i],
#   #                              b = df_posterior_Linear_1$b[i])
#   fd_Linear[i,] <- get_Linear(x = input_quantity,a = df_posterior_Linear$a[i],
#                                b = df_posterior_Linear$b[i])
#   fd_Hill[i,] <- get_Hill(x = input_quantity, a = df_posterior_Hill$a[i],
#                            b = df_posterior_Hill$b[i], c = df_posterior_Hill$c[i],
#                            g = df_posterior_Hill$g[i])
#   fd_Pow[i,] <- get_Power(x = input_quantity, a = df_posterior_Pow$a[i],
#                          b = df_posterior_Pow$b[i], g = df_posterior_Pow$g[i])
#   fd_Expo5[i,] <- get_Expo5(x = input_quantity, a = df_posterior_Expo5$a[i],
#                            b = df_posterior_Expo5$b[i], c = df_posterior_Expo5$c[i],
#                            d = df_posterior_Expo5$d[i])
# }

y_weighted <- vector("list",length = Nmethods)
for(j in 1:Nmethods){
  # y_weighted[[j]] <- y_Linear_0 * wts_1[1,j] + y_Linear_1 * wts_1[2,j] +
  #   y_Linear_2 * wts_1[3,j]
  # y_weighted[[j]] <- y_Hill_4 * wts_1[1,j] + y_Linear_2 * wts_1[2,j]
  y_weighted[[j]] <- fd_Linear * wts_1[1,j] + fd_Hill * wts_1[2,j]+
    fd_Pow * wts_1[3,j] + fd_Expo5 * wts_1[4,j]
}
y_error <- lapply(y_weighted,function(x) x - y_original)

# Evaluation metrics-------
RMSE_iter <- lapply(y_error,function(x) {sqrt(apply((x)^2,MARGIN = 1,mean))})
RMSE <- lapply(RMSE_iter,median)

MAE_iter <- lapply(y_error,function(x) {
  apply(abs(x),MARGIN = 1,mean)
})
MAE <- lapply(MAE_iter,median)

MRB_iter <- lapply(y_weighted,function(x){
  apply(abs(x - y_original) / y_original,MARGIN = 1,mean)
})
MRB <- lapply(MRB_iter,median)

CompareMethods <- as.data.frame(rbind(unlist(RMSE),unlist(MAE),unlist(MRB)))

min_Index <- apply(CompareMethods,1,function(x) {
  which(x == min(x))
})
colnames(CompareMethods) <- Method_list
CompareMethods$Choice <- Method_list[min_Index[min_Index != 8]]
rownames(CompareMethods) <- c("RMSE","MAE","MRB")
