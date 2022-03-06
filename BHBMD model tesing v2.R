# Settings----

# global values
seed <- 47405
iter <- 1E4
warmup <- 0.5
Nchain <- 4
thin <- 1

# load utility functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Hierarchical model utilities.R")

# load data
load("Hill meta base data.Rdata")

# Validation of Single Study Model-----
#* use example data
dataList <- list(
  G = nrow(df_example),
  n = df_example$N,
  ymean = df_example$ymeanL,
  ysd = df_example$ysdL,
  dose = df_example$dose,
)

prior_example <- list(
  prior_sigma = c(0,2.5),
  prior_a = get_prior_a(df_example$ymeanL,df_example$ysdL),
  prior_b = get_prior_b(df_example$ymeanL,df_example$ysdL,df_example$dose),
  prior_c = c(0,15),
  prior_g = c(1,15)
)

stanfit_Hill_example <- stan(model_code = modelstring_Hill_one,
                          data = append(dataList,prior_example),
                          iter = iter,
                          chains = 3,
                          warmup = warmup * iter,
                          thin = thin)

print(stanfit_Hill_example)


#* testing single study

# mannual input of study index
test_study <- bladder_meta_all %>% filter(Index == 1)
dataList <- list(
  G = nrow(test_study),
  n = test_study$N,
  # ymeanL = get_ymeanL(bladder_chen$RR,bladder_chen$ysd),
  ymean =test_study$ymeanL,
  # ysd = get_ysdL(bladder_chen$RR,c(0.1,bladder_chen$ysd[c(2,3)])) ,
  ysd = test_study$ysdL,
  dose = test_study$dose_minmax_single)

prior_single_study <- list(
  prior_sigma = c(0,2.5),
  prior_a = get_prior_a(test_study$RR,test_study$ysd),
  prior_b = get_prior_b(test_study$RR,test_study$ysd,test_study$dose_minmax_single),
  prior_c = c(0,15),
  prior_g = c(1,15)
)


stanfit_Hill_test_study <- stan(model_code = modelstring_Hill_one,
                          data = append(dataList,prior_single_study),
                          iter = iter,
                          chains = Nchain,
                          warmup = warmup * iter,
                          thin = thin)

print(stanfit_Hill_test_study,
      pars = c("a","b","c","g","sigma"))

# Validation of meta hierarchical model----

# use Allen data
# compare with figure S10 in Shao 2021 SI page 12
# because the data are available
# only left is to validate the model script

modelname <- "Hill_meta_fun_vec_rag"
dataName <- "meta"
df_input <- df_bladder_meta

# *Optional Specifications-------
# Testing on the order of study
OrderNew <- c("Bates 2004","Bates 1995","Wu 2013","Meliker 2010",
              "Steinmaus 2003","Steinmaus 2013",
              "Chen 2010","Sawada 2013M","Sawada 2013F")
OrderNew <- unique(df_input$Study)
df_input <- ReorderData(df_input = df_input,OrderNew = OrderNew)

# set initial values
inits <- list(
  list(c = 0.1,g = 10),
  list(c = 0.1,g = 10),
  list(c = 1, g = 2),
  list(c = 1, g = 2)
)

# *Hill Testing---------
# stan sampling and save samples in a Rdata file
stanfit <- get_metaStanfit(modelname = modelname,
                       df_input = df_input,
                       dataName = dataName)
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))

# stanfit diagnose
stan_diag(object = stanfit,
          information = c(
                          # "sample"
                          # "treedepth"
                          # "stepsize"
                          "divergence"
                          ),
          chain = 0)
stan_par(object = stanfit,
         par = "c",
         chain = 0)
pairs(stanfit,
      pars = c(
        # "a",
        "c",
        "g",
        "mu_a",
        "sigma_a",
        "mu_b",
        "sigma_b",
        "sigma"
      ))

# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(stanfit)

# when a and b with overarching lognormal dist
df_posteriors <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b),
  # c = rlnorm(iterEff, meanlog = mu_c, sdlog = sigma_c),
  # c = rlnorm(iterEff, meanlog = mu_c, sdlog = sigma_c),
  # g = rlnorm(iterEff, meanlog = mu_g, sdlog = sigma_g)
  # cast = exp(mu_c + sigma_c^2/2),
  aast = exp(mu_a + sigma_a^2/2),
  bast = exp(mu_b + sigma_b^2/2)
  )
write.csv(df_posteriors,
          file = paste0("posterior summary ",modelname,"_",dataName,".csv"),
          row.names = F)

# Determine which chains of posteriors to be plotted
ChainID <- 1:Nchain
df_extracted <- ExtractPosteriorbyChain(df_posteriors = df_posteriors,
                                         chains = ChainID)
write.csv(df_extracted,
          file = paste0("Posterior of chain ",ChainID,".csv"),
          row.names = F)
# examine the summary stats of extracted posterior
apply(df_extracted,MARGIN = 2, FUN = function(x)
  c(SD = sd(x), summary(x), Q = quantile(x,c(0.05,0.95))))

# Traceplot
plot(df_extracted$sigma,type = "l")
plot(df_extracted$c,type = "l")
plot(df_extracted$g,type = "l")
plot(df_extracted$mu_a,type = "l")
plot(df_extracted$sigma_a,type = "l")
plot(df_extracted$mu_b,type = "l")
plot(df_extracted$sigma_b,type = "l")

# display study specific curves
showcurve_specific_Hill(df_posteriors = df_extracted,df_input = df_input,
                        aname = "a",bname = "b", cname = "c", gname = "g",
                        xup = NULL,yup = NULL)

# overarching median curve
showcurve_overarching_Hill(df_posteriors = df_extracted,
                           aname = "a",bname = "b",cname = "c",gname = "g",
                      df_input = df_input,
                      xup = 10,yup = NULL)

# Hill BMD estimation

df_BMD <- getBMD_meta_Hill(BMR,Ref,
                           df_posteriors = df_extracted,
                           df_input = df_input)
# Note shown curves are fitted RR
# Have not been approximated by Hill curves
# Results replicate what in Shao's figures
# Model is validated


# Prior comparison----


stanfit_bladder_meta_Hill_NC1_dADD <- stan(model_code = modelstring_Hill_meta_NC1,
                                           # file =  "modelstring_Hill_meta_NC1.stan",
                                           data = append(dataAll,NC1prior),
                                           iter = iter,
                                           chains = 3,
                                           warmup = warmup * iter,
                                           thin = 1)
stanfit_bladder_meta_Hill_NC2_dADD <- stan(model_code = modelstring_Hill_meta_NC2,
                                           # file =  "modelstring_Hill_meta_NC2.stan",
                                           data = append(dataAll,NC2prior),
                                           iter = iter,
                                           chains = 3,
                                           warmup = warmup * iter,
                                           thin = 1)
stanfit_bladder_meta_Hill_CT_gLnorm_dADD <- stan(model_code = modelstring_Hill_meta_CT_gLnorm,
                                                 # file =  "modelstring_Hill_meta_CT_gLnorm.stan",
                                                 data = append(dataAll,CTgLnormprior),
                                                 iter = iter,
                                                 chains = 3,
                                                 warmup = warmup * iter,
                                                 thin = 1)
stanfit_bladder_meta_Hill_CT_gLnorm2_dADD <- stan(model_code = modelstring_Hill_meta_CT_gLnorm2,
                                                  # file =  "modelstring_Hill_meta_CT_gLnorm2.stan",
                                                  data = append(dataAll,CTgLnorm2prior),
                                                  iter = iter,
                                                  chains = 3,
                                                  warmup = warmup * iter,
                                                  thin = 1)
stanfit_bladder_meta_Hill_CT_gUnif_dADD <- stan(model_code = modelstring_Hill_meta_CT_gUnif,
                                                # file =  "modelstring_Hill_meta_CT_gUnif.stan",
                                                data = append(dataAll,CTgUnifprior),
                                                iter = iter,
                                                chains = 3,
                                                warmup = warmup * iter,
                                                thin = 1)


#* Hill all data----

modelname <- "Hill_meta_Shao_fd_old_ver2"
dataName <- "Allen"
df_input <- df_bladder_Allen

# stan sampling and save samples in a Rdata file
assign("stanfit",
       get_metaStanfit(modelname = modelname,
                       df_input = df_input,
                       dataName = dataName))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))


# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(stanfit)
df_posteriors <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b)
)
summary(df_posteriors[,c("a","b")])

# display study specific curves
showcurve_specific_Hill(df_posteriors = df_posteriors,
                   df_input = df_input)

# overarching median curve
showcurve_overarching_Hill(df_posteriors = df_posteriors,
                      df_input = df_input)
# Hill BMD estimation

df_BMD <- getBMD_meta_Hill(BMR,Ref,
                           df_posteriors = df_posteriors,
                           df_input = df_input)
#* Hill Case Control only----

modelname <- "Hill_meta_Shao_fd_old_ver2"
dataName <- "Allen"
df_input <- df_bladder_Allen

# stan sampling and save samples in a Rdata file
assign("stanfit",
       get_metaStanfit(modelname = modelname,
                       df_input = df_input,
                       dataName = dataName))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))


# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(stanfit)
df_posteriors <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b)
)
summary(df_posteriors[,c("a","b")])

# display study specific curves
showcurve_specific(df_posteriors = df_posteriors,
                   df_input = df_input)

# overarching median curve
showcurve_overarching(df_posteriors = df_posteriors,
                      df_input = df_input)
# Hill BMD estimation
BMR <- 0.01/100
Ref <- 0.03

df_BMD <- getBMD_meta_Hill(BMR,Ref,
                           df_posteriors = df_posteriors,
                           df_input = df_input)
#* Linear all data------

modelname <- "Linear_meta_fun_vec_rag"
dataName <- "meta"
df_input <- df_bladder_meta

# stan sampling and save samples in a Rdata file
assign("stanfit",
       get_metaStanfit(modelname = modelname,
                       df_input = df_input,
                       dataName = dataName))
# If stanfit data is available from storage, read it:
load(paste0("stanfit_",modelname,"_",dataName,".RData"))


# store posteriors in a dataframe and print summary stats
# overarching distribution parameter must be calculated mannually
fitsummary <- get_fitsummary(stanfit)
df_posteriors <- fitsummary %>% mutate(
  a = rlnorm(iterEff, meanlog = mu_a, sdlog = sigma_a),
  b = rlnorm(iterEff, meanlog = mu_b, sdlog = sigma_b)
)
summary(df_posteriors[,c("a","b")])

# display study specific curves
showcurve_specific_Linear(df_posteriors = df_posteriors,
                        df_input = df_input,
                        aname = "a", bname = "b")

# overarching median curve
showcurve_overarching_Linear(df_posteriors = df_posteriors,
                           df_input = df_input,
                           aname = "a" , bname = "b")

# BMD estimation

df_BMD <- getBMD_meta_Linear(BMR,Ref,
                   df_posteriors = df_posteriors,
                   df_input = df_input)
