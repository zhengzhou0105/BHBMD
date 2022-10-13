

# Linear----
prior_Linear_meta_fun_vec_rag <- list(
  prior_sigma = c(0,5),
  prior_mu_a = c(-1,1),
  prior_sigma_a = c(0,2),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2)
)

# Linear_meta_0 no hier
prior_Lognormal_Summary_Linear_meta_0 <- list(
  prior_sigma = c(0,2.5),
  prior_a = c(-1,1),
  prior_b = c(-2,2)
)

# Linear_meta_1 hier on a
prior_Lognormal_Summary_Linear_meta_1 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(0,10),
  prior_sigma_a = c(0,2.5),
  prior_b = c(1,2)
)

# Linear_meta_2 hier on a and b
prior_Lognormal_Summary_Linear_meta_2 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(-1*a_bound,a_bound),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2.5)
)

# Linear_meta_3 hier on a and b alternative priors
prior_Lognormal_Summary_Linear_meta_3 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  prior_scale_b =c(1)
)

# Power----

prior_Lognormal_Summary_Power_meta_3 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2.5),
  prior_mu_g = c(-3,3),
  prior_sigma_g = c(0,2.5),
  g_constr = c(1,18)
)

prior_Lognormal_Summary_Power_meta_4 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  prior_scale_b =c(1),
  prior_mu_g = c(-1,1),
  prior_sigma_g = c(1),
  g_constr = c(1,18)
)

# MM--------
prior_Lognormal_Summary_MM_meta_3 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2.5),
  prior_mu_c = c(-1,1),
  prior_sigma_c = c(0,2.5)
)

prior_Lognormal_Summary_MM_meta_4 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  prior_scale_b =c(1),
  prior_mu_c = c(-2,2),
  prior_sigma_c = c(1)
)

# Hill------
prior_Hill_meta_rag_loglik <- list(
  prior_sigma = c(0,5),
  prior_mu_a = c(-1,1),
  prior_sigma_a = c(0,2),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2),
  prior_c = c(0,10),
  prior_g = c(0,18),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_2 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(0,10),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(0,10),
  prior_sigma_b = c(0,2.5),
  prior_c = c(0,5),
  prior_g = c(0,18),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_1 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(0,10),
  prior_sigma_a = c(0,2.5),
  prior_b = c(-3,5),
  prior_c = c(0,5),
  prior_g = c(0,18),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_0 <- list(
  prior_sigma = c(0,2.5),
  prior_a = c(-2,3),
  prior_b = c(-3,5),
  prior_c = c(0,5),
  prior_g = c(0,18),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_3 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(0,10),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(0,10),
  prior_sigma_b = c(0,2.5),
  prior_mu_c = c(0,10),
  prior_sigma_c = c(0,2.5),
  prior_g = c(0,18),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_4 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2.5),
  prior_mu_c = c(-3,3),
  prior_sigma_c = c(0,2.5),
  prior_mu_g = c(-3,3),
  prior_sigma_g = c(0,2.5),
  g_constr = c(0,18)
)

prior_Lognormal_Summary_Hill_meta_5 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  prior_scale_b = c(1),
  prior_mu_c = c(-2,2),
  prior_sigma_c = c(1),
  prior_mu_g = c(-1,1),
  prior_sigma_g = c(1),
  g_constr = c(1,18)
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
# 
# # centered g follows lnorm
# prior_Hill_meta_CT_gLnorm <- list(
#   prior_sigma = c(0,2.5),
#   prior_c = c(0,15),
#   prior_g = c(1,15),
#   prior_mu_a = c(0,10),
#   prior_sigma_a = c(0,2.5),
#   prior_s_b = c(0,15),
#   prior_r_b = c(0,15)
# )
# 
# # centered g follows lnorm2
# prior_Hill_meta_CT_gLnorm2 <- list(
#   prior_sigma = c(0,2.5),
#   prior_c = c(0,15),
#   prior_g = c(1,15),
#   prior_sigma_a = c(0,2.5),
#   prior_s_b = c(0,15),
#   prior_r_b = c(0,15)
# )
# 
# # CT gUnif
# prior_Hill_meta_CT_gUnif <- list(
#   prior_sigma = c(0,2.5),
#   prior_c = c(0,15),
#   prior_g = c(1,15),
#   prior_a_a = c(0,10),
#   # prior_b_a = c(0,max(get_prior_meta_a(index = SStudy,
#   #                                      ymean = ymean,
#   #                                      ysd = ysd))*3),
#   prior_b_a = c(0,30),
#   prior_s_b = c(0,15),
#   prior_r_b = c(0,15)
# )

# Exponential------
prior_Lognormal_Summary_Expo5_meta_4 <- list(
  prior_sigma = c(0,2.5),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(0,2.5),
  prior_mu_b = c(-3,3),
  prior_sigma_b = c(0,2.5),
  prior_mu_cminusone = c(0,2.5),
  prior_sigma_cminusone = c(0,2.5),
  prior_mu_c = c(-3,3),
  prior_sigma_c = c(0,2.5),
  prior_mu_d = c(-3,3),
  prior_sigma_d = c(0,2.5),
  d_constr = c(0,18)
)

prior_Lognormal_Summary_Expo5_meta_5 <- list(
  prior_sigma = c(1),
  prior_mu_a = c(-a_bound,a_bound),
  prior_sigma_a = c(1),
  prior_scale_b =c(1),
  prior_scale_c =c(1),
  c_lower = c(1),
  prior_mu_d = c(-1,1),
  prior_sigma_d = c(1),
  d_constr = c(0,18)
)
