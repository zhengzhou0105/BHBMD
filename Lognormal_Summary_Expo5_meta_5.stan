functions{
  real Normal_Summary_lpdf(real[] ymean,real[] ysd,real[] Nsub,real[] fd,real sigma){
    return -sum(to_vector(Nsub)*log(sigma)+ ((to_vector(Nsub)-1).*square(to_vector(ysd))+(to_vector(Nsub)).*square(to_vector(ymean)-to_vector(fd)))/(2*sigma^2));
  }
  
      // Power function vectorized
    real[] pow_vec(real[] x,real power){
      real x_power[size(x)];
      
      for(n in 1:size(x)){
        x_power[n] = pow(x[n],power);
      }
      return x_power;
    }
}

data{
  int<lower=0> N;                           // number of observation
  int<lower=0> S;                           // number of studies
  int<lower=0> G[S];                        // number of dose groups per study
  int<lower=1,upper=S> Index[N];            // study indexes
  
  real ymeanL[N];                           // ymean at log scale per dose group
  real ysdL[N];                             // ysd at log scale per dose group
  real Nsub[N];                             // number of subjects per dose group
  real doseZ[N];                             // standardized dose per dose group
  // real dose[N];                             // dose per dose group
  
  real prior_sigma;
  real prior_mu_a[2];
  real prior_sigma_a;
  // real prior_a_b[2];
  // real prior_b_b[2];
  real prior_scale_b;
  // real prior_mu_c[2];
  // real prior_sigma_c[2];
  real prior_scale_c;
  real c_lower;
  real prior_mu_d[2];
  real prior_sigma_d;
  real d_constr[2];
}

parameters{
  real<lower=0> sigma;
  real<lower=0> a[S];
  real<lower=0> b[S];
  real<lower=c_lower> c[S];
  real<lower=d_constr[1],upper=d_constr[2]> d[S];
  
  real mu_a;
  real<lower=0> sigma_a;
  // real<lower=0> a_b;
  // real<lower=0> b_b;
  real<lower=0> scale_b;
  // real mu_c;
  // real<lower=0> sigma_c;
  real<lower=0> scale_c;
  real mu_d;
  real<lower=0> sigma_d;
}


model{
  int pos;

  // sigma ~ uniform(prior_sigma[1],prior_sigma[2]);
  sigma ~ exponential(prior_sigma);
  // print("sigma = ",target());
  
  mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
  // print("mu_a = ",target());
  
  // sigma_a ~ uniform(prior_sigma_a[1],prior_sigma_a[2]);
  sigma_a ~ exponential(prior_sigma_a);
  // print("sigma_a = ",target());
  
  // a_b ~ uniform(prior_a_b[1],prior_a_b[2]);
  // print("sigma = ",target());
  
  // b_b ~ uniform(prior_b_b[1],prior_b_b[2]);
  // print("sigma = ",target());
  
  scale_b ~ exponential(prior_scale_b);
  // print("scale_b = ",target());
  
  // mu_c ~ uniform(prior_mu_c[1],prior_mu_c[2]);
  // print("sigma = ",target());
  
  // sigma_c ~ uniform(prior_sigma_c[1],prior_sigma_c[2]);
  // print("sigma = ",target());

  scale_c ~ exponential(prior_scale_c);
  // print("scale_c = ",target());
  
  mu_d ~ uniform(prior_mu_d[1],prior_mu_d[2]);
  // print("mu_d = ",target());
  
  // sigma_d ~ uniform(prior_sigma_d[1],prior_sigma_d[2]);
  sigma_d ~ exponential(prior_sigma_d);
  // print("sigma_d = ",target());
  
  a ~ lognormal(mu_a,sigma_a);
  // print("a = ",target());
  
  b ~ cauchy(0,scale_b);
  // print("b = ",target());
  
  c ~ cauchy(0,scale_c);
  // print("c = ",target());
  
  d ~ lognormal(mu_d,sigma_d);
  // print("d = ",target());

  pos = 1;
  
  for(s in 1:S){
    real seg_ymeanL[G[s]];
    real seg_ysdL[G[s]];
    real seg_Nsub[G[s]];
    real seg_fdL[G[s]];
    real seg_doseZ[G[s]];
    // real seg_dose[G[s]];
    
    seg_ymeanL = segment(ymeanL,pos,G[s]);
    seg_ysdL = segment(ysdL,pos,G[s]);
    seg_Nsub = segment(Nsub,pos,G[s]);
    seg_doseZ = segment(doseZ,pos,G[s]);
    // seg_dose = segment(dose,pos,G[s]);
    
    seg_fdL = to_array_1d(log( a[s] * (c[s] - (c[s] -1) * exp(-1 * to_vector(pow_vec(to_array_1d(b[s] * to_vector(seg_doseZ)),d[s])))) ));  // use fdL instead of log(fd)
    // seg_fdL = to_array_1d(log( a[s] * (c[s] - (c[s] -1) * exp(-1 * to_vector(pow_vec(to_array_1d(b[s] * to_vector(seg_doseZ)),d[s])))) ));  // use fdL instead of log(fd)
    target += Normal_Summary_lpdf(seg_ymeanL | seg_ysdL, seg_Nsub, seg_fdL, sigma);
    
    pos = pos + G[s];
  }
}

generated quantities{
  real log_lik[N];
  real fd[N];

  for(n in 1:N){
    fd[n] = a[Index[n]] *(c[Index[n]] - (c[Index[n]]-1) * exp(-(b[Index[n]]*doseZ[n])^d[Index[n]]));                          // simulate fd for double check
    // fd[n] = a[Index[n]] *(c[Index[n]] - (c[Index[n]]-1) * exp(-(b[Index[n]]*doseZ[n])^d[Index[n]]));                          // simulate fd for double check
    log_lik[n] = Normal_Summary_lpdf({ymeanL[n]} | {ysdL[n]},{Nsub[n]},{log(fd[n])},sigma);
  }
}
