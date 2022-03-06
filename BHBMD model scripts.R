# Model explanation----

# Hill_meta_Shao partial parameters a & b from Shao 2021 but vectorized
# Hill_meta_Shao_fd based on  Hill_meta_Shao, Hill function moved to transformed parameter block to display fitted f(d)
# Hill_meta_NC1 gUnif: partial parameters g&v with non-centered parameterization g = g0 + g[j], g0 ~ uniform
# Hill_meta_NC2 gLnorm: based on NC1 but g0 ~ Lnorm
# Hill_meta_CT_gLnorm: partial parameters a&b centere parameterization a[j] ~ Lnorm(mu_a,sigma_a), b[j]~gamma
# Hill_meta_CT_gLnorm2: based on CT_gLnorm but a[j]~Lnorm(0,sigma_a)
# Hill_meta_CT_gUnif: based on CT_gLnorm but a[j]~Unif(lower_a,upper_a)
# Hill_meta_full_CT

# Linear_meta full hierarchical structure

# Log density function---------
# Summarized Continuous data under lognormality assumption
"
functions{
  // Define log normal BMD log density function for summarized continuous response
  real SumlnormBMD_lpdf(real meanlog, real Nsub, real fd, real sdlog,real sigma){
    return -(
      Nsub/2*log(square(sigma)) + ( Nsub*square((meanlog-log(fd))) + (Nsub-1)*square(sdlog) )/(2*square(sigma))
      );
  }
}
"
# Individual Continous data under lognormality assumption
"
// Define log normal BMD log density function for individual continuous response
  real IndlnormBMD_ldf(real yL, real fd,real sigma){
    return 
  }
"
# Hill Single study----

# works for summarized data given the current project context
modelstring_Hill_one = "
  data {
    int G;                   // number of dose group
    vector[G] n;             // number of subjects per dose group
    vector[G] ymean;         // sample mean at log scale
    vector[G] ysd;           // sample SD at log scale
    vector[G] dose;          // dose
    
    // hyperpriors
    real prior_sigma[2];
    real prior_a[2];
    real prior_b[2];
    real prior_c[2];
    real prior_g[2];
  }
  
  parameters {
    real<lower=0> sigma;     // constant variance assumed for lnorm dist
    // continuous Hill parameter BMDS ver. parameterization
    real a;                  // intercept
    real<lower=0> b;                  // maximum change
    real<lower=0> c;         // half-maximal dose
    real<lower=1,upper=18> g;// power, restricted; for unrestricted, lower = 0
  }

  model {
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    a ~ uniform(prior_a[1],prior_a[2]);
    b ~ uniform(prior_b[1],prior_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // approximate MLE solution to loglikelihood when datadata 
    // are summarized in sample mean and SD. BMDS user manual
    for(i in 1:G){
          target += -(n[i] * log(sigma^2) / 2 + ((n[i] - 1) * (ysd[i])^2) /
          (2 * sigma^2) + (n[i] * (ymean[i] - log(a + ((b * dose[i]^g) / 
          (c^g + dose[i]^g))))^2) / (2 * sigma^2));
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


# Hill Meta Centered----
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



# Hill Shao 2021's model------

modelstring_Hill_meta_Shao_old ="
  data{
    int<lower=0> N;
    int<lower=0> S;
    int<lower=0,upper=S> Index[N];
    
    vector[N] ymeanL;
    vector[N] ysdL;
    vector[N] dose;
    row_vector<lower=0>[N] Nsub;
    
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
  
  transformed data{
    vector[N] dosed;
    for(n in 1:N){
      dosed[n] = (dose[n] - min(dose)) / (max(dose) / min(dose));
    }
  }
  
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
    vector[N] fd;
    
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    //for(n in 1:N){
    //  fd[n] = a[Index[n]] + (b[Index[n]] * dosed[n]^g) / (c^g + dosed[n]^g);
    //}
    //target += -( Nsub/2*log(square(sigma)) +  ( Nsub*square((ymeanL-log(fd))) + (Nsub-1)*square(ysdL) ) /
    //    (2*square(sigma)) );
    for(n in 1:N){
      target += -( Nsub[n]/2*log(square(sigma)) +  ( Nsub[n]*square((ymeanL[n]-log(
      // Hill equation
      a[Index[n]] + (b[Index[n]] * dosed[n]^g) / (c^g + dosed[n]^g)
      ))) + (Nsub[n]-1)*square(ysdL[n]) ) /   (2*square(sigma)) );
    }
  }

"
modelstring_Hill_meta_Shao_old_ver3 <- "
 data{
    int<lower=0> N;               // number of observation
    int<lower=0> S;               // number of studies
    int<lower=0> Gmax;            // largest number of dose groups
    int<lower=0,upper=Gmax> G[S]; // number of dose groups

    vector[N] ymeanL;             // ymean log
    vector[N] ysdL;               // ysd log
    vector<lower=0>[N] Nsub;      // number of subjects
    row_vector[N] dose;           // dose
    
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
  
  transformed data{
    row_vector[N] dosed;          // transformed dose
    dosed = ( dose - min(dose) ) / ( max(dose) - min(dose) );
  }
  
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
    vector[N] fd;              
    int pos;              

    
    // Priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // Parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    // Likelihood
    pos = 1;
    for(s in 1:S){
      for(gg in pos:(pos+G[s]-1)){
        target += -( Nsub[gg]/2*log(sigma^2) +  ( Nsub[gg]*(ymeanL[gg]-log(
          a[s] + ( b[s] * dosed[gg]^g ) / ( c^g + dosed[gg]^g )
        ))^2 + (Nsub[gg]-1)*ysdL[gg]^2 ) / (2*sigma^2) );
      }
    }
  }

"


modelstring_Hill_meta_Shao_fd_old_ver2 <- "
 data{
    int<lower=0> N;               // number of observation
    int<lower=0> S;               // number of studies
    int<lower=0> Gmax;            // largest number of dose groups
    int<lower=0,upper=Gmax> G[S]; // number of dose groups

    vector[N] ymeanL;             // ymean log
    vector[N] ysdL;               // ysd log
    row_vector<lower=0>[N] Nsub;      // number of subjects
    row_vector[N] dose;           // dose
    
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
  
  transformed data{
    row_vector[N] dosed;          // transformed dose
    for(n in 1:N){
      dosed[n] = (dose[n] - min(dose)) / (max(dose) / min(dose));
    }
  }
  
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
    vector[N] fd;              
    int pos;              

    
    // Priors
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // Parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    // Hill equation
    pos = 1;                          
    for(s in 1:S){
      for(gg in pos:(pos+G[s]-1)){
        fd[gg] = a[s] + ( b[s] * dosed[gg]^g ) / ( c^g + dosed[gg]^g );
      }
      pos = pos + G[s];
    }

    // Likelihood
      target += -( Nsub/2*log(square(sigma)) +  ( Nsub*square((ymeanL-log(fd))) + (Nsub-1)*square(ysdL) ) /
        (2*square(sigma)) );
  }
  
  
"
  
modelstring_Hill_meta_Shao_nonrag <- "
  data{
    int<lower=0> S;              // number of study
    int<lower=0> Gmax;              // number of dose group each study
    int<lower=0> N;              // number of records
    
    row_vector[Gmax] mtx_dose[S];                // dose
    vector[Gmax] mtx_Nsub[S];                // number of subjects
    vector[Gmax] mtx_ymeanL[S];              // log mean
    vector[Gmax] mtx_ysdL[S];                // log SD
    
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
  
  transformed data{
    row_vector[Gmax] mtx_dosed[S];
    matrix[S,Gmax] matrix_dose;
    for(s in 1:S){
      matrix_dose[s,] = mtx_dose[s];
    }
    for(s in 1:S){
      for(gg in 1:Gmax){
        mtx_dosed[s,gg] = (mtx_dose[s,gg] - min(matrix_dose))/(max(matrix_dose) - min(matrix_dose));
      }
    }  
  }
  
  parameters{
    real<lower=0> sigma;
    
    vector<lower=0>[S] a;
    vector<lower=0>[S] b;
    
    real mu_a;
    real<lower=0> sigma_a;
    real mu_b;
    real<lower=0> sigma_b;
    
    real<lower=0> c;
    real<lower=g_constr> g;
    
  }
  
  model{
    vector[Gmax] fd[S];
    
    // Priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

    // Parameters
    for(s in 1:S){
      a[s] ~ lognormal(mu_a,sigma_a);
      b[s] ~ lognormal(mu_b,sigma_b);
      for(gg in 1:Gmax){
        fd[s,gg] = a[s] + ( b[s] * mtx_dosed[s,gg]^g ) / ( c^g + mtx_dosed[s,gg]^g );

      }
    }
    
    // Likelihood
    for(s in 1:S){
      for(gg in 1:Gmax){
        target += -( mtx_Nsub[s,gg]/2*log(sigma^2) +  ( mtx_Nsub[s,gg]*(mtx_ymeanL[s,gg]-log(
          //a[s] + ( b[s] * mtx_dosed[s,gg]^g ) / ( c^g + mtx_dosed[s,gg]^g )
          fd[s,gg]
        ))^2 + (mtx_Nsub[s,gg]-1)*mtx_ysdL[s,gg]^2 ) / (2*sigma^2) );
      }
    }
  }
  
  
"
modelstring_Hill_meta_Shao_rag <- "
  data{
    int<lower=0> N;                  // number of data
    int<lower=0> S;                  // number of study
    int<lower=0> Gmax;               // largest number of dose groups
    int<lower=0,upper=Gmax> G[S];   // number of dose groups
    
    vector[N] dose;                 // dose
    vector<lower=0>[N] Nsub;        // number of subjects
    vector[N] ymeanL;                // RR mean at log scale
    vector[N] ysdL;                  // RR sd at log scale
    
    // hyperparameters
    real prior_sigma[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
    real prior_c[2];
    real prior_g[2];
    //real prior_mu_c[2];
    //real prior_sigma_c[2];
    real g_constr;
  }

  //transformed data{
  //  vector[N] dosed;
  //  for(n in 1:N){
  //    dosed[n] = (dose[n] - min(dose)) / (max(dose) - min(dose));
  //  }
  //}
  
  parameters{
    real<lower=0> sigma;
    
    vector<lower=0>[S] a;                  // study specific a
    real mu_a;                             // meanlog of a
    real<lower=0> sigma_a;                // sdlog of a
    vector<lower=0>[S] b;                 // study specific b
    real mu_b;                            // meanlog of b
    real<lower=0> sigma_b;                // sdlog of b
    
    real<lower=0> c;
    //vector<lower=0>[S] c;                    // threshold
    //real mu_c;
    //real<lower=0> sigma_c;
    
    real<lower=g_constr> g;             // g
  }
  
  model{
    int pos;
    vector[N] fd;
    
    // priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    //mu_c ~ uniform(prior_mu_c[1],prior_mu_c[2]);
    //sigma_c ~ uniform(prior_sigma_c[1],prior_sigma_c[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

    // parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    pos = 1;
    for(s in 1:S){
      for(gg in pos:(pos+G[s]-1)){
        fd[gg] = a[s] + (b[s] * dose[gg]^g) / (c^g + dose[gg]^g);;
      }
      pos = pos + G[s];
    }
    
    // likelihood
    for(n in 1:N){
        target += -( Nsub[n]/2*log(sigma^2) +  ( Nsub[n]*(ymeanL[n]-log(fd[n]))^2 + (Nsub[n]-1)*ysdL[n]^2 ) / (2*sigma^2) );
    }
  }

"

# 
# modelstring_Hill_meta_Shao_fd_old <- "
#   data {
#     // this structure requires data in long format
#     int<lower=1> N;                 // number of data point
#     int<lower=1> S;                 // number of study
#     int<lower=1,upper=S> Index[N];         // study Index
#     int N[N];              // number of subjects at each dose level
#     vector[N] dose;        // dose
#     vector[N] ymean;       //  response SD at each dose level
#     vector[N] ysd;         // response SD at each dose level
#     // hyperparameters
#     real prior_sigma[2];
#     real prior_mu_a[2];
#     real prior_sigma_a[2];
#     real prior_mu_b[2];
#     real prior_sigma_b[2];
#     real prior_c[2];
#     real prior_g[2];
#   }
#   
#   parameters {
#     real<lower = 0> sigma;
#     // continuous Hill parameter BMDS ver. parameterization
#     real<lower = 0> c;                  // half-maximal dose
#     real<lower = 0> g;// power, restricted; 
#     vector<lower = 0>[S] a;             // study specific a
#     vector<lower = 0>[S] b;             // study specific b
#     real mu_a;               // overarching meanlog of a
#     real<lower=0> sigma_a;            // overarching sdlog of a
#     real mu_b;              // overarching meanlog of b
#     real<lower=0> sigma_b;          // overarching sdlog of b
#   }
#   
#   transformed parameters{
#     real<lower = 0> fd[N];
#     for(i in 1:N){
#       fd[i] = a[Index[i]] + (b[Index[i]] * dose[i]^g) / (c^g + dose[i]^g);
#     }
#   }
#   
#   model {
#     // priors
#     mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
#     sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
#     mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
#     sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
#     sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
#     c ~ uniform(prior_c[1],prior_c[2]);
#     g ~ uniform(prior_g[1],prior_g[2]);
#     
#     // parameters
#     a ~ lognormal(mu_a, sigma_a);
#     b ~ lognormal(mu_b, sigma_b);
#     
#     // likelihood
#     for(i in 1:N){
#       target += -(n[i] * log(sigma^2) / 2 + (n[i] * (ymean[i] - log(fd[i]))^2 + 
#       (n[i] - 1) * ysd[i]^2) / (2 * sigma^2));
# 
#     }
#   }
#   
# 
# "

modelstring_Hill_meta_Shao_fdgen_old <- "
  data {
    // this structure requires data in long format
    int<lower=0> N;                 // number of data point
    int<lower=1> S;                 // number of study
    int<lower=1,upper=S> Index[N];         // study Index
    int Nsub[N];              // number of subjects at each dose level
    real dose[N];        // dose
    real ymeanL[N];       //  y mean log at each dose level
    real<lower=0> ysdL[N];         // y SD log at each dose level
    
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
  
  transformed data{
    real dosed[N];
    for(n in 1:N){
      dosed[n] = (dose[n] - min(dose)) / (max(dose) - min(dose));
    }
  }
  
  parameters {
    real<lower = 0> sigma;
    // continuous Hill parameter BMDS ver. parameterization
    real<lower = 0> c;                  // half-maximal dose
    real<lower = g_constr> g;// power, restricted; 
    real<lower = 0> a[S];             // study specific a
    real<lower = 0> b[S];             // study specific b
    real mu_a;               // overarching meanlog of a
    real<lower=0> sigma_a;            // overarching sdlog of a
    real mu_b;              // overarching meanlog of b
    real<lower=0> sigma_b;          // overarching sdlog of b
  } 
  
  model {
    // priors
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    
    // parameters
    a ~ lognormal(mu_a, sigma_a);
    b ~ lognormal(mu_b, sigma_b);

    // likelihood   
    for(i in 1:N){
      target += -(Nsub[i] * log(sigma^2) / 2 + ((Nsub[i] - 1) * ysdL[i]^2 + Nsub[i] * (ymeanL[i] - log(
      // Hill equation
      a[Index[i]] + (b[Index[i]] * dose[i]^g) / (c^g + dose[i]^g)
      ))^2) / (2 * sigma^2));
    }
  }

  generated quantities{
    real<lower = 0> fd[N];
    for(i in 1:N){
      fd[i] = a[Index[i]] + (b[Index[i]] * dose[i]^g) / (c^g + dose[i]^g);
    }
  }
  
  
"

# modelstring_Hill_meta_Shao_old1 <- "
#   data {
#     int<lower=0> N;                   // # of total records
#     int<lower=1> S;                   // # of studies
#     int<lower=1,upper=N> nstart[S];   // starting point of each study
#     int<lower=1,upper=N> nend[S];     // end point of each study
#     vector<lower=0>[N] dose;          // doses
#     vector<lower=0>[N] N;             // number of observatiosn at each dose group
#     vector[N] ymean;                  // ymean at each dose group
#     vector<lower=0>[N] ysd;           // ysd at each dose group
#     int<lower=1,upper=S> Index[N];    // study Index
# 
#       // hyperparameters
#       real prior_sigma[2];
#       real prior_mu_a[2];
#       real prior_sigma_a[2];
#       real prior_mu_b[2];
#       real prior_sigma_b[2];
#       real prior_c[2];
#       real prior_g[2];
# 
#   }
#   
#   parameters {
#     real<lower=0> sigma;              // error under homogenity assumption
#     vector<lower=0>[S] a;             // study-specific intercept a
#     vector<lower=0>[S] b;             // study-specific slope b
#     real<lower=0> c;                  // half maixum
#     real<lower=1,upper=15> g;         // power      
#     
#     real mu_a;                        // meanlog of a ~ Lnorm
#     real<lower=0> sigma_a;            // sdlog of a ~ Lnorm
#     real mu_b;                        // meanlog of b ~ Lnorm
#     real<lower=0> sigma_b;            // sdlog of b ~ Lnorm
#   }
#   
#   model {
#     // priors
      # mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
      # sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
      # mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
      # sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
      # sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
      # c ~ uniform(prior_c[1],prior_c[2]);
      # g ~ uniform(prior_g[1],prior_g[2]);
#     
#     // likelihood
#     a ~ lognormal(mu_a, sigma_a);
#     b ~ lognormal(mu_b, sigma_b);
# 
#     for(s in 1:S){
#       for(i in nstart[s]:nend[s]){
#         target += -(n[i] * log(sigma^2) / 2 + (n[i] * (ymean[i] - log(
#         // f(d) function
#         a[s] + (b[s] * dose[i]^g) / (c^g + dose[i]^g)
#         ))^2 + 
#         (n[i] - 1) * ysd[i]^2) / (2 * sigma^2))^2;
#       }
#     }
#   }
# 
# "
# 
# 
# modelstring_Hill_meta_Shao_fd_old1 <- "
#   data {
#     int<lower=0> N;                   // # of total records
#     int<lower=1> S;                   // # of studies
#     int<lower=1,upper=N> nstart[S];   // starting point of each study
#     int<lower=1,upper=N> nend[S];     // end point of each study
#     vector<lower=0>[N] dose;          // doses
#     vector<lower=0>[N] N;             // number of observatiosn at each dose group
#     vector[N] ymean;                  // ymean at each dose group
#     vector<lower=0>[N] ysd;           // ysd at each dose group
# 
#       // hyperparameters
#       real prior_sigma[2];
#       real prior_mu_a[2];
#       real prior_sigma_a[2];
#       real prior_mu_b[2];
#       real prior_sigma_b[2];
#       real prior_c[2];
#       real prior_g[2];
# 
#   }
#   
#   parameters {
#     real<lower=0> sigma;              // error under homogenity assumption
#     vector<lower=0>[S] a;             // study-specific intercept a
#     vector<lower=0>[S] b;             // study-specific slope b
#     real<lower=0> c;                  // half maixum
#     real<lower=1,upper=15> g;         // power      
#     
#     real mu_a;                        // meanlog of a ~ Lnorm
#     real<lower=0> sigma_a;            // sdlog of a ~ Lnorm
#     real mu_b;                        // meanlog of b ~ Lnorm
#     real<lower=0> sigma_b;            // sdlog of b ~ Lnorm
#   }
#   
#   transformed parameters {
#     vector<lower=0>[N] fd;            // f(d) at each dose group
#     //for(i in 1:N){
#     //  fd[i] = a[Index[i]] + (b[Index[i]] * dose[i]^g) / (c^g + dose[i]^g);
#     //}
#     for(s in 1:S){
#       for(i in nstart[s]:nend[s]){
#         fd[i] = a[s] + (b[s] * dose[i]^g) / (c^g + dose[i]^g);
#       }
#     }
#   }
#   
#   model {
#     // priors
#     mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
#     sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
#     mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
#     sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
#     sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
#     c ~ uniform(prior_c[1],prior_c[2]);
#     g ~ uniform(prior_g[1],prior_g[2]);
# 
#     // likelihood
#     a ~ lognormal(mu_a, sigma_a);
#     b ~ lognormal(mu_b, sigma_b);
#     
#     for(s in 1:S){
#       for(i in nstart[s]:nend[s]){
#         target += -(n[i] * log(sigma^2) / 2 + (n[i] * (ymean[i] - log(
#         // f(d) function
#         fd[i]
#         ))^2 + 
#         (n[i] - 1) * ysd[i]^2) / (2 * sigma^2));
#       }
#     }
#   }
# 
# "

modelstring_Hill_meta_Shao_fdgen_old1 <- "
  data {
    int<lower=0> N;                   // # of total records
    int<lower=0> S;                   // # of studies
    int<lower=1,upper=N> nstart[S];   // starting point of each study
    int<lower=1,upper=N> nend[S];     // end point of each study
    vector<lower=0>[N] dose;          // doses
    vector<lower=0>[N] Nsub;             // number of observatiosn at each dose group
    vector[N] ymeanL;                  // ymean log at each dose group
    vector<lower=0>[N] ysdL;           // ysd log at each dose group
    
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

  transformed data{
    real dosed[N];
    for(n in 1:N){
      dosed[n] = (dose[n] - min(dose)) / (max(dose) - min(dose));
    }
  }
  
  parameters {
    real<lower=0> sigma;              // error under homogenity assumption
    vector<lower=0>[S] a;             // study-specific intercept a
    vector<lower=0>[S] b;             // study-specific slope b
    real<lower=0> c;                  // half maixum
    real<lower=g_constr> g;         // power      
    
    real mu_a;                        // meanlog of a ~ Lnorm
    real<lower=0> sigma_a;            // sdlog of a ~ Lnorm
    real mu_b;                        // meanlog of b ~ Lnorm
    real<lower=0> sigma_b;            // sdlog of b ~ Lnorm
  }
  
  model {
    // priors
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

    // parameters
    a ~ lognormal(mu_a, sigma_a);
    b ~ lognormal(mu_b, sigma_b);
    
    // likelihood
    for(s in 1:S){
      for(i in nstart[s]:nend[s]){
        target += -(Nsub[i] * log(sigma^2) / 2 + (Nsub[i] * (ymeanL[i] - log(
        // f(d) function
        a[s] + (b[s] * dose[i]^g) / (c^g + dose[i]^g)
        ))^2 + 
        (Nsub[i] - 1) * ysdL[i]^2) / (2 * sigma^2));
      }
    }
  }
    
  generated quantities{
    vector<lower=0>[N] fd;            // f(d) at each dose group
    for(s in 1:S){
      for(i in nstart[s]:nend[s]){
        fd[i] = a[s] + (b[s] * dose[i]^g) / (c^g + dose[i]^g);
      }
    }
  }
    

"

modelstring_Hill_meta_Shao_matrix <- "
  data {
    // data must be in matrix forms- each group per row
    // data
    int<lower=0> S;                          // # of studies
    int<lower=1> Gmax;                     // largest # of dose groups
    int<lower=1,upper=Gmax> G[S];          // # of dose group in each study
    matrix[S,Gmax] mtx_dose;             // matrix of dose
    matrix[S,Gmax] mtx_ymeanL;             // matrix of y mean at log scale
    matrix<lower=0>[S,Gmax] mtx_ysdL;     // matrix of y SD at log scale
    matrix<lower=0>[S,Gmax] mtx_N;        // matrix of # of subjects
    
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
  
  transformed data{
    matrix[S,Gmax] mtx_doseD;            // matrix of standaridzed dose
    real mindose;                         // min of dose
    real maxdose;                        // max of dose
    mindose = min(mtx_dose);
    maxdose = max(mtx_dose);
    mtx_doseD = (mtx_dose - mindose) / (maxdose - mindose);
  }
  
  parameters {
    real<lower=0> sigma;        // SD of noise distribution
    
    vector<lower=0>[S] a;       // study specific a, RR>0
    vector<lower=0>[S] b;       // study specific b, monotonic increasing
    real<lower=0> c;            // half maximal dose >0
    real<lower=g_constr> g;   // power, restricted bayesian
    
    // Overarching parameters based on lognomral prior distributions
    real mu_a;                  // mean log of a
    real<lower=0> sigma_a;      // sd log of a
    real mu_b;                  // mean log of b
    real<lower=0> sigma_b;      // sd log of b
  }
  
  model {
    // priors
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

    // parameters
    a ~ lognormal(mu_a,sigma_a);      // vectorized
    b ~ lognormal(mu_b,sigma_b);      // vectorized
    
    // likelihood function
    for(s in 1:S){
      for(i in 1:G[s]){
        target += -(mtx_N[s,i] * log(sigma^2) / 2 + (mtx_N[s,i] * (mtx_ymeanL[s,i] - log(
        // f(d) function
        a[s] + (b[s] * mtx_doseD[s,i]^g) / (c^g + mtx_doseD[s,i]^g)
        ))^2 + 
        (mtx_N[s,i] - 1) *  mtx_ysdL[s,i]^2) / (2 * sigma^2));
      }
    }
  }

"

modelstring_Hill_meta_fun_index <- "
  functions{
    // non vectorized loglik function under lognormality assumption
    real lognormal_summary_lpdf(
      real ymeanL,
      real ysdL,
      real Nsub,
      real fd,
      real sigma
    ){
      return -Nsub*log(sigma)-(Nsub*(ymeanL-fd)^2+(Nsub-1)*ysdL^2)/(2*sigma^2);
    }
  }
  
  data{
    int<lower=0> N;           // number of dose groups
    int<lower=0> S;           // number of studies
    int<lower=0,upper=S> Index[N];   // study index
    vector[N] dose;           // doses
    vector[N] Nsub;           // number of subjects
    vector[N] ymeanL;         // summary mean log
    vector[N] ysdL;           // summary sd log
    
    // hyperparameters
    real g_constr;
    real prior_sigma[2];
    real prior_c[2];
    real prior_g[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
  }
  
  //transformed data{
  //  real Ntotal;
  //  Ntotal = sum(Nsub);
  //}
  
  parameters{
    real<lower=0> sigma;
    real<lower=0> c;
    real<lower=g_constr> g;
    
    vector<lower=0>[S] a; 
    real mu_a;
    real<lower=0> sigma_a;
    
    vector<lower=0>[S] b;
    real mu_b;
    real<lower=0> sigma_b;
  }
  
  model{
    // priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    
    // Parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    // likelihood
    for(n in 1:N){
      ymeanL[n] ~ lognormal_summary_lpdf(ysdL[n],Nsub[n],
      // insert Hill equation
        log(a[Index[n]] + (b[Index[n]] * dose[n]^g)/(c^g + dose[n]^g)),
      sigma);
    }
  }

"
modelstring_Hill_meta_fun_index_fd <- "

  functions{
    // non vectorized loglik function under lognormality assumption
    real lognormal_summary_lpdf(
      real ymeanL,
      real ysdL,
      real Nsub,
      real fd,
      real sigma
    ){
      return -Nsub*log(sigma)-(Nsub*(ymeanL-fd)^2+(Nsub-1)*ysdL^2)/(2*sigma^2);
    }
  }
  
  data{
    int<lower=0> N;           // number of dose groups
    int<lower=0> S;           // number of studies
    int<lower=0,upper=S> Index[N];   // study index
    vector[N] dose;           // doses
    vector[N] Nsub;           // number of subjects
    vector[N] ymeanL;         // summary mean log
    vector[N] ysdL;           // summary sd log
    
    // hyperparameters
    real g_constr;
    real prior_sigma[2];
    real prior_c[2];
    real prior_g[2];
    real prior_mu_a[2];
    real prior_sigma_a[2];
    real prior_mu_b[2];
    real prior_sigma_b[2];
  }
  
  //transformed data{
  //  real Ntotal;
  //  Ntotal = sum(Nsub);
  //}
  
  parameters{
    real<lower=0> sigma;
    real<lower=0> c;
    real<lower=g_constr> g;
    
    vector<lower=0>[S] a; 
    real mu_a;
    real<lower=0> sigma_a;
    
    vector<lower=0>[S] b;
    real mu_b;
    real<lower=0> sigma_b;
  }
  
  model{
    vector[N] fd;
    
    // priors
    sigma ~ cauchy(prior_sigma[1],prior_sigma[2]);
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);
    mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
    sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
    mu_b ~ uniform(prior_mu_b[1],prior_mu_b[2]);
    sigma_b ~ uniform(prior_sigma_b[1],prior_sigma_b[2]);
    
    // Parameters
    a ~ lognormal(mu_a,sigma_a);
    b ~ lognormal(mu_b,sigma_b);
    
    for(n in 1:N){
      fd[n] = a[Index[n]] + (b[Index[n]] * dose[n]^g)/(c^g + dose[n]^g);
    }
    // likelihood
    for(n in 1:N){
      ymeanL[n] ~ lognormal_summary_lpdf(ysdL[n],Nsub[n],log(fd[n]),sigma);
    }
  }
  

"

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