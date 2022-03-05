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
    real prior_c[2];
    real prior_g[2];
    real g_constr;
  }

  parameters{
    real<lower=0> sigma;                    // residual standard deviation
    real<lower=0> c;                        // Hill half maximal
    real<lower=g_constr> g;                 // Hill power
    
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
    c ~ uniform(prior_c[1],prior_c[2]);
    g ~ uniform(prior_g[1],prior_g[2]);

      // parameters
    a ~ lognormal(mu_a,sigma_a);              // vectorized
    b ~ lognormal(mu_b,sigma_b);

    pos = 1;                                   // ragged structure
    for(s in 1:S){                            // study specific loglik
      segment(ymeanL,pos,G[s]) ~ lognormal_summary(    // use self defined pdf
        segment(ysdL,pos,G[s]),
        segment(Nsub,pos,G[s]),
        // Insert dose response function
        // Hill equation
        a[s]+(b[s]*segment(pow_vec(dose,g),pos,G[s]))./(c^g+segment(pow_vec(dose,g),pos,G[s])),
        sigma                                     
      );
      pos = pos + G[s];
    }
  }
  
  