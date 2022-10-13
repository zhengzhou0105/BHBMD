# Log Probability Density function for summary data---------
# Summarized Continuous data under lognormality assumption
# require homogenous variance, fd, sample size, summarized meanlog and sdlog 
"
  functions{
    // define probability function for summary continous response under lognormality assumption
    // allow vectorization
    real lognormal_summary_lpdf (
      vector ymeanL,                               // summary response mean at log scale
      vector ysdL,                                 // summary response sd at log scale
      vector Nsub,                                 // number of subjects per group
      vector fd,                                   // predicted response given dose response function
      real sigma){                                 // residual standard deviation
      ## vector muliplication .*
      return -sum(Nsub*log(sigma)+ ((Nsub-1).*square(ysdL)+(Nsub).*square(ymeanL-log(fd)))/(2*sigma^2));
    }
    
    // Power function vectorized
    vector pow_vec(vector x,real power){
      vector[dims(x)[1]] x_power;       // create a new vector x_power of the same length as x
      
      for(n in 1:dims(x)[1]){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
  }
"

# Hill Single study----

# works for summarized data given the current project context
modelstring_Hill_one = "

  functions{
    // define probability function for summary continous response under lognormality assumption
    // allow vectorization
    real lognormal_summary_lpdf (
      vector ymeanL,                               // summary response mean at log scale
      vector ysdL,                                 // summary response sd at log scale
      vector Nsub,                                 // number of subjects per group
      vector fd,                                   // predicted response given dose response function
      real sigma){                                 // residual standard deviation
      ## vector muliplication .*
      return -sum(Nsub*log(sigma)+ ((Nsub-1).*square(ysdL)+(Nsub).*square(ymeanL-log(fd)))/(2*sigma^2));
    }
    
    // Power function vectorized
    vector pow_vec(vector x,real power){
      vector[dims(x)[1]] x_power;       // create a new vector x_power of the same length as x
      
      for(n in 1:dims(x)[1]){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
  }
  
  data{
    int<lower=0> G;                   // number of dose groups
    vector<lower=0>[G] Nsub;          // sample size per dose group
    vector<lower=0>[G] ymeanL;        // summarized mean log per dose group
    vector<lower=0>[G] ysdL;          // summarized sdlog per dose group
    vector<lower=0>[G] dose;          // dose of each group
    
    // hyperparameters
    real prior_a[2];
    real prior_b[2];
    real prior_c[2];
    real prior_g[2];
    real prior_sigma[2];
  }
  
  parameters{
    real<lower=0> sigma;
    real a;
    real b;
    real<lower=0> c;
    real<lower=0,upper=18> g;
  }
  
  model{
    // Priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    a ~ uniform(prior_a[1],prior_a[2]);
    b ~ uniform(prior_b[1],prior_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // Likelihood
    target += lognormal_summary_lpdf(
        ymeanL |
        ysdL,
        Nsub,
        a + (b * pow_vec(dose,g))./(c^g + pow_vec(dose,g)),
        sigma
      ); 
  }
  
  generated quantities{
    vector[G] log_lik;
    for(i in 1:G){
      log_lik[i] = lognormal_summary_lpdf(
        ymeanL[i] |
        ysdL[i],
        Nsub[i],
        a + (b * dose[i]^g)/(c^g + dose[i]^g),
        sigma
      );
    }
  }

"

# Hill Partial Hierarchical model------

modelstring_Hill_meta_fun_vec_rag <- "
  //
  functions{
    // define probability function for summary continous response under lognormality assumption
    // allow vectorization
    real lognormal_summary_lpdf (
      vector ymeanL,                               // summary response mean at log scale
      vector ysdL,                                 // summary response sd at log scale
      vector Nsub,                                 // number of subjects per group
      vector fd,                                   // predicted response given dose response function
      real sigma){                                 // residual standard deviation
      return -sum(Nsub*log(sigma)+ ((Nsub-1).*square(ysdL)+(Nsub).*square(ymeanL-log(fd)))/(2*sigma^2));
    }
    
    // Power function vectorized
    vector pow_vec(vector x,real power){
      vector[dims(x)[1]] x_power;
      
      for(n in 1:dims(x)[1]){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
  }

  data{
    int<lower=0> N;                                        // number of observations
    int<lower=0> S;                                        // number of studies
    int<lower=0> G[S];                                     // number of dose groups in each study
    
    vector[N] dose;                                        // doses
    vector[N] Nsub;                                        // number of subjects per group
    vector[N] ymeanL;                                      // response mean at log scale
    vector[N] ysdL;                                        // response sd at log scale
    
    // hyperparameters
    real prior_sigma[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
    real prior_c[2];
    real prior_g[2];
    real g_constr;
  }
  
  //transformed data{
  //  vector[N] dosed;                                   // [0,1] standardized dose
  //  dosed = (dose - min(dose))/(max(dose) - min(dose));
  //}

  parameters{
    real<lower=0> sigma;
    real<lower=0> c;
    real<lower=g_constr> g;
    
    vector<lower=0>[S] a;
    vector<lower=0>[S] b;
    
    real mu_a;
    real<lower=0> sigma_a;
    
    real mu_b;
    real<lower=0> sigma_b;
  }

  model{
    int pos;
    
      // priors
    sigma ~ uniform(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

      // parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    pos = 1;
    for(s in 1:S){
      segment(ymeanL,pos,G[s]) ~ lognormal_summary(
        segment(ysdL,pos,G[s]),
        segment(Nsub,pos,G[s]),
        a[s]+(b[s]*segment(pow_vec(dose,g),pos,G[s]))./(c^g+segment(pow_vec(dose,g),pos,G[s])),
        sigma
      );
      pos = pos + G[s];
    }
  }
  
  

"

modelstring_Hill_meta_fun_vec_rag_test <- "
  //
  functions{
    // define probability function for summary continous response under lognormality assumption
    // allow vectorization
    real lognormal_summary_lpdf (
      vector ymeanL,                               // summary response mean at log scale
      vector ysdL,                                 // summary response sd at log scale
      vector Nsub,                                 // number of subjects per group
      vector fd,                                   // predicted response given dose response function
      real sigma){                                 // residual standard deviation
      return -sum(Nsub*log(sigma)+ ((Nsub-1).*square(ysdL)+(Nsub).*square(ymeanL-log(fd)))/(2*sigma^2));
    }
    
    // Power function vectorized
    vector pow_vec(vector x,real power){
      vector[dims(x)[1]] x_power;
      
      for(n in 1:dims(x)[1]){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
  }

  data{
    int<lower=0> N;                                        // number of observations
    int<lower=0> S;                                        // number of studies
    int<lower=0> G[S];                                     // number of dose groups in each study
    
    vector[N] dose;                                        // doses
    vector[N] Nsub;                                        // number of subjects per group
    vector[N] ymeanL;                                      // response mean at log scale
    vector[N] ysdL;                                        // response sd at log scale
    
    // hyperparameters
    real prior_sigma[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
    real prior_c[2];
    real prior_g[2];
    real g_constr;
  }
  
  //transformed data{
  //  vector[N] dosed;                                   // [0,1] standardized dose
  //  dosed = (dose - min(dose))/(max(dose) - min(dose));
  //}

  parameters{
    real<lower=0> sigma;
    real<lower=0> c;
    real<lower=g_constr> g;
    
    vector<lower=0>[S] a;
    vector<lower=0>[S] b;
    
    real mu_a;
    real<lower=0> sigma_a;
    
    real mu_b;
    real<lower=0> sigma_b;
  }

  model{
    int pos;
    
      // priors
    sigma ~ uniform(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

      // parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    pos = 1;
    for(s in 1:S){
    
      // segmentation of data
      vector[G[s]] seg_ymeanL;
      vector[G[s]] seg_ysdL;
      vector[G[s]] seg_Nsub;
      vector[G[s]] seg_dose;
      vector[G[s]] seg_fd;
      
      seg_ymeanL = segment(ymeanL,pos,G[s]);
      seg_ysdL = segment(ysdL,pos,G[s]);
      seg_Nsub = segment(Nsub,pos,G[s]);
      seg_dose = segment(dose,pos,G[s]);
      
      seg_fd = a[s]+(b[s]*pow_vec(seg_dose,g))./(c^g+pow_vec(seg_dose,g));
      
      seg_ymeanL ~ lognormal_summary(
        seg_ysdL,
        seg_Nsub,
        seg_fd,
        sigma
      );
      pos = pos + G[s];
    }
  }
  
  

"


# Hill Full Hierarchical-----
modelstring_Hill_meta_full_CT <- "
  data{
    // this structure requires data in long format
    int<lower=1> N;                 // number of data point
    int<lower=1> S;                 // number of study
    int<lower=0> Index[N];         // study Index
    int<lower=0> n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymean;       //  response SD at each dose level
    vector[N] ysd;         // response SD at each dose level
    
    //  hyperpriors
    vector[2] prior_sigma;
    vector[2] prior_mu_a;
    vector[2] prior_sigma_a;
    vector[2] prior_mu_b;
    vector[2] prior_sigma_b;
    vector[2] prior_mu_c;
    vector[2] prior_sigma_c;
    vector[2] prior_mu_g;
    vector[2] prior_sigma_g;
  }

  parameters{
    real<lower = 0> sigma;
    vector<lower = 0>[S] a;             // study specific a intercept
    vector<lower = 0>[S] b;             // study specific b slope
    vector<lower = 0>[S] c;            // study specific c half maximum
    vector<lower = 1, upper = 18>[S] g;           // study specific g   power
    
    real mu_a;               // overarching meanlog of a
    real<lower=0> sigma_a;            // overarching sdlog of a
    real mu_b;              // overarching meanlog of b
    real<lower=0> sigma_b;          // overarching sdlog of b
    real mu_c;
    real<lower=0> sigma_c;
    real mu_g;
    real<lower=0> sigma_g;
  }

  transformed parameters{
    vector<lower = 0>[N] fd;
    for(i in 1:N){
      fd[i] = a[Index[i]] + (b[Index[i]] * dose[i]^g[Index[i]]) / 
      (c[Index[i]]^g[Index[i]] + dose[i]^g[Index[i]]);
    }
  }

  model{
      // priors
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    mu_c ~ uniform(prior_mu_c[1],prior_mu_c[2]);
    sigma_c ~ uniform(prior_sigma_c[1],prior_sigma_c[2]);
    mu_g ~ uniform(prior_mu_g[1],prior_mu_g[2]);
    sigma_g ~ uniform(prior_sigma_g[1],prior_sigma_g[2]);

    
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);

    // parameters
    a ~ lognormal(mu_a, sigma_a);
    b ~ lognormal(mu_b, sigma_b);
    c ~ lognormal(mu_c, sigma_c);
    g ~ lognormal(mu_g, sigma_g);
    
    // likelihood
    for(i in 1:N){
      target += -(n[i] * log(sigma^2) / 2 + (n[i] * (ymean[i] - log(fd[i]))^2 + 
      (n[i] - 1) * ysd[i]^2) / (2 * sigma^2));

    }
  }

"

# Linear-------
modelstring_Linear_meta_fun_vec_rag <- "
// Author: Zheng Zhou
// Date: Jan 17 2022

  functions{
    // define probability function for summary continous response under lognormality assumption
    // allow vectorization
    real lognormal_summary_lpdf (
      vector ymeanL,                               // summary response mean at log scale
      vector ysdL,                                 // summary response sd at log scale
      vector Nsub,                                 // number of subjects per group
      vector fd,                                   // predicted response given dose response function
      real sigma){                                 // residual standard deviation
      return -sum(Nsub*log(sigma)+ ((Nsub-1).*square(ysdL)+(Nsub).*square(ymeanL-log(fd)))/(2*sigma^2));
    }
    
    // Power function vectorized
    vector pow_vec(vector x,real power){
      vector[dims(x)[1]] x_power;
      
      for(n in 1:dims(x)[1]){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
  }

  data{
    int<lower=0> N;                                        // number of observations
    int<lower=0> S;                                        // number of studies
    int<lower=0> G[S];                                     // number of dose groups in each study
    
    vector[N] dose;                                        // doses
    vector[N] Nsub;                                        // number of subjects per group
    vector[N] ymeanL;                                      // response mean at log scale
    vector[N] ysdL;                                        // response sd at log scale
    
    // hyperparameters
    real prior_sigma[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
  }

  parameters{
    real<lower=0> sigma;                    // residual standard deviation

    vector<lower=0>[S] a;                   // study specific intercept
    vector<lower=0>[S] b;                   // study specific slope
     
    real mu_a;                              // a's overarching dist- meanlog
    real<lower=0> sigma_a;                  //a's overarching dist- sdlog
    
    real mu_b;                              // b's overarching dist -meanlog
    real<lower=0> sigma_b;                  // b's overarching dist- sdlog
  }

  model{
    int pos;                                // position counter
    
      // priors
    sigma ~ uniform(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);

      // parameters
    a ~ lognormal(mu_a,sigma_a);              // vectorized
    b ~ lognormal(mu_b,sigma_b);

    pos = 1;                                   // ragged structure
    for(s in 1:S){                            // study specific loglik
      segment(ymeanL,pos,G[s]) ~ lognormal_summary(    // use self defined pdf
        segment(ysdL,pos,G[s]),
        segment(Nsub,pos,G[s]),
        // Insert dose response function
        // Linear equation
        a[s]+ b[s] * segment(dose,pos,G[s]),
        sigma                                     
      );
      pos = pos + G[s];
    }
  }
  
  

"
# Hill Meta Non-centered-----
# non centered models have not been fixed
# check out https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html

#* NC1 gUnif-------
modelstring_Hill_meta_NC_gUnif = "
  data {
    // this structure requires data in long format
    int N;                 // number of data point
    int S;                 // number of study
    int s[N];              // study Index
    int n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymeanL;       // sample mean at log scale
    vector[N] ysdL;         // sample SD at log scale
  }
  
  parameters {
    real<lower=0> sigma;
    // continuous Hill parameter BMDS ver. parameterization
    real<lower=0> k;                  // half-maximal dose
    real<lower=1,upper=18> p;// power, restricted; for unrestricted, lower = 0
    
    real<lower=0> g0;                 // overarching g;
    real<lower=0> v0;                 // overarching v;
    vector<lower=0>[S] g1;             // study specific g
    vector<lower=0>[S] v1;             // study specific v
  } 
  
  // transformed parameters{ 
  //   vector[N] fd;                     // hierarhcial Hill function
  //   vector<lower=0>[S] g_marg;        // marginalized g
  //   vector<lower=0>[S] v_marg;        // marginalized v
  //   for(i in 1:N){
  //       g_marg[s[i]] = g0 + g[s[i]] * gd;
  //       v_marg[s[i]] = v0 + v[s[i]] * vd;
  //       fd[i] = g_marg[s[i]] + (v_marg[s[i]]*dose[i]^p)/(k^p+dose[i]^p);   
  //   }
  // }
  
  model {
    sigma ~ cauchy(0,2.5);
    g0 ~ uniform(prior_g[1],prior_g[2]);
    v0 ~ uniform(prior_v[1],prior_v[2]);
    g1 ~ normal(0,10);
    v1 ~ normal(0,10);
    
    k ~ uniform(prior_k[1],prior_k[2]);
    p ~ uniform(prior_p[1],prior_p[2]); 
        
    for(i in 1:N){
      target += -(n[i]*log(sigma^2)/2 + ((n[i]-1)*ysdL[i]^2)/(2*sigma^2) +
      (n[i]*(ymeanL[i]-log((g0 + g1[s[i]]) + ((v0 + v1[s[i]])*dose[i]^p)/(k^p+dose[i]^p)))^2)/(2*sigma^2));
    }
  }
  

"



#* NC2 gLnorm-------
modelstring_Hill_meta_NC_gLnorm = "
  data {
    // this structure requires data in long format
    int N;                 // number of data point
    int S;                 // number of study
    int s[N];              // study Index
    int n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymeanL;       // sample mean at log scale
    vector[N] ysdL;         // sample SD at log scale
    
    // hierarchical priors
    vector[2] prior_sigma;
    vector[2] prior_g;
    vector[2] prior_v;
    vector[2] prior_p;
    vector[2] prior_k;
  }
  
  parameters {
    real<lower=0> sigma;
    // continuous Hill parameter BMDS ver. parameterization
    real<lower=0> k;                  // half-maximal dose
    real<lower=1,upper=18> p;// power, restricted; for unrestricted, lower = 0
    
    real<lower=0> g0;                 // overarching g;
    real<lower=0> v0;                 // overarching v;
    vector<lower=0>[S] g1;             // study specific g
    vector<lower=0>[S] v1;             // study specific v
  } 
  
  // transformed parameters{ 
  //   vector[N] fd;                     // hierarhcial Hill function
  //   vector<lower=0>[S] g_marg;        // marginalized g
  //   vector<lower=0>[S] v_marg;        // marginalized v
  //   for(i in 1:N){
  //       g_marg[s[i]] = g0 + g[s[i]] * gd;
  //       v_marg[s[i]] = v0 + v[s[i]] * vd;
  //       fd[i] = g_marg[s[i]] + (v_marg[s[i]]*dose[i]^p)/(k^p+dose[i]^p);   
  //   }
  // }
  
  model {
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    g0 ~ lognormal(prior_g[1],prior_g[2]);
    v0 ~ uniform(prior_v[1],prior_v[2]);
    g1 ~ normal(0,10);
    v1 ~ normal(0,10);
    
    k ~ uniform(prior_k[1],prior_k[2]);
    p ~ uniform(prior_p[1],prior_p[2]); 
        
    for(i in 1:N){
      target += -(n[i]*log(sigma^2)/2 + ((n[i]-1)*ysdL[i]^2)/(2*sigma^2) +
      (n[i]*(ymeanL[i]-log((g0 + g1[s[i]]) + ((v0 + v1[s[i]])*dose[i]^p)/(k^p+dose[i]^p)))^2)/(2*sigma^2));
    }
  }
  

"


#* Hill Meta Centered----
#* g follows lnorm----
modelstring_Hill_meta_CT_gLnorm = "
  data {
    // this structure requires data in long format
    int<lower=1> N;                 // number of data point
    int<lower=1> S;                 // number of study
    int<lower=0> Index[N];         // study Index
    int n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymean;       //  response SD at each dose level
    vector[N] ysd;         // response SD at each dose level
    
    //  hyperpriors
    vector[2] prior_sigma;
    vector[2] prior_g;
    vector[2] prior_c;
    vector[2] prior_mu_a;           // meanlog of a
    vector[2] prior_sigma_a;        // sdlog of a
    vector[2] prior_s_b;    // shape of b
    vector[2] prior_r_b;    // rate of b
  }
  
  parameters {
    real<lower=0> sigma;
    // continuous Hill parameter BMDS ver. parameterization
    real<lower=0> c;                  // half-maximal dose
    real<lower=1,upper=15> g;// power, restricted; for unrestricted, lower = 0
    vector<lower=0>[S] a;             // vector of study specific a
    vector<lower=0>[S] b;             // vector of study specific b
    real mu_a;                        // overarching meanlog of g
    real<lower=0> sigma_a;                     // overarching sdlog of g
    real<lower=0> s_b;                         // overarching shape of v
    real<lower=0> r_b;                         // overarching rate of v
  } 
  
  model {
    // priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ cauchy(prior_sigma_a[1],prior_sigma_a[2]);
    s_b ~ uniform(prior_s_b[1],prior_s_b[2]);
    r_b ~ uniform(prior_r_b[1],prior_r_b[2]);

    // partial hierarchical parameters
    a ~ lognormal(mu_a, sigma_a);
    b ~ gamma(s_b, r_b);
    
    // likelihood    
    for (i in 1:N) {
      target += -(n[i] * log(sigma^2) / 2 + ((n[i] - 1) * ysd[i]^2) /
      (2 * sigma^2) + (n[i] * (ymean[i] - log(a[Index[i]] + ( b[Index[i]] *
      dose[i]^g)/(c^g + dose[i]^g)))^2)/(2 * sigma^2));
    }
  }
  

"



#* g follows lnorm2 fixed mean-------
modelstring_Hill_meta_CT_gLnorm2 = "
  data {
    // this structure requires data in long format
    int<lower=1> N;                 // number of data point
    int<lower=1> S;                 // number of study
    int<lower=0> Index[N];         // study Index
    int n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymean;       //  response SD at each dose level
    vector[N] ysd;         // response SD at each dose level
    
    //  hyperpriors
    vector[2] prior_sigma;
    vector[2] prior_g;
    vector[2] prior_c;
    vector[2] prior_sigma_a;
    vector[2] prior_s_b;
    vector[2] prior_r_b;
  }
  
  parameters {
    real<lower=0> sigma;
    // continuous Hill parameter BMDS ver. parameterization
    real<lower=0> c;                  // half-maximal dose
    real<lower=1,upper=15> g;// power, restricted; for unrestricted, lower = 0
    vector<lower=0>[S] a;             // vector of study specific a
    vector<lower=0>[S] b;             // vector of study specific b
    real<lower=0> sigma_a;                     // overarching sdlog of g
    real<lower=0> s_b;                         // overarching shape of v
    real<lower=0> r_b;                         // overarching rate of v
  } 
  
  model {
    // priors
    sigma_a ~ cauchy(prior_sigma_a[1],prior_sigma_a[2]);
    s_b ~ uniform(prior_s_b[1],prior_s_b[2]);
    r_b ~ uniform(prior_r_b[1],prior_r_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // parameters
    a ~ lognormal(0, sigma_a);
    b ~ gamma(s_b, r_b);
        
    // likelihood
    for(i in 1:N){
      target += -(n[i] * log(sigma^2)/2 + ((n[i] - 1) * ysd[i]^2) / 
      (2 * sigma^2) + (n[i] * (ymean[i] - log(a[Index[i]] + ( b[Index[i]] *
      dose[i]^g) / (c^g + dose[i]^g)))^2) / (2 * sigma^2));
    }
  }
  

"


#* g follows uniform-------
modelstring_Hill_meta_CT_gUnif = "
  data {
    // this structure requires data in long format
    int<lower=1> N;                 // number of data point
    int<lower=1> S;                 // number of study
    int<lower=0> Index[N];         // study Index
    int n[N];              // number of subjects at each dose level
    vector[N] dose;        // dose
    vector[N] ymean;       //  response SD at each dose level
    vector[N] ysd;         // response SD at each dose level
    
    //  hyperparameters
    vector[2] prior_sigma;
    vector[2] prior_g;
    vector[2] prior_c;
    vector[2] prior_a_a;
    vector[2] prior_b_a;
    vector[2] prior_s_b;
    vector[2] prior_r_b;
  }
  
  parameters {
    real<lower=0> sigma;
    // continuous Hill parameter BBMD ver. parameterization
    real<lower=0> c;                  // half-maximal dose
    real<lower=1,upper=15> g;// power, restricted; for unrestricted, lower = 0
    vector<lower=0>[S] a;             // study specific a
    vector<lower=0>[S] b;             // study specific b
    real<lower=0> a_a;                        // overarching lower bound of a
    real<lower=0> b_a;                         // overarching upper bound of a
    real<lower=0> s_b;                         // overarching shape of b
    real<lower=0> r_b;                         // overarching rate of b
  } 
  
  model {
    // priors
    a_a ~ uniform(prior_a_a[1],prior_a_a[2]);
    b_a ~ cauchy(prior_b_a[1],prior_b_a[2]);
    s_b ~ uniform(prior_s_b[1],prior_s_b[2]);
    r_b ~ uniform(prior_r_b[1],prior_r_b[2]);
    
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // parameters
    a ~ uniform(a_a, b_a);
    b ~ gamma(s_b, r_b);
        
    // likelihood
    for(i in 1:N){
      target += -(n[i] * log(sigma^2)/2 + ((n[i] - 1) * ysd[i]^2) / 
      (2 * sigma^2) + (n[i] * (ymean[i] - log(a[Index[i]] + ( b[Index[i]] * 
      dose[i]^g) / (c^g + dose[i]^g)))^2) / (2 * sigma^2));
    }
  }
  

"


