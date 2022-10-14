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
  
  real doseZ[N];                             // dose per dose group
  // real dose[N];
  real ymeanL[N];                           // ymean at log scale per dose group
  real ysdL[N];                             // ysd at log scale per dose group
  real Nsub[N];                             // number of subjects per dose group
  
  // real prior_sigma[2];
  real prior_sigma;
  real prior_mu_a[2];
  // real prior_sigma_a[2];
  real prior_sigma_a;
  // real prior_a_b[2];
  // real prior_b_b[2];
  real prior_scale_b;
}

parameters{
  real<lower=0> sigma;
  real<lower=0> a[S];
  real<lower=0> b[S];
  
  real mu_a;
  real<lower=0> sigma_a;
  // real<lower=0> a_b;
  // real<lower=0> b_b;
  real<lower=0> scale_b;
}

model{
  int pos;

  // sigma ~ lognormal(prior_sigma[1],prior_sigma[2]);
  sigma ~ exponential(prior_sigma);
  // print("target sigma =",target());
  
  mu_a ~ uniform(prior_mu_a[1],prior_mu_a[2]);
  // print("target mu_a =",target());
  
  // sigma_a ~ lognormal(prior_sigma_a[1],prior_sigma_a[2]);
  sigma_a ~ exponential(prior_sigma_a);
  // print("target sigma_a =",target());
  
  // a_b ~ uniform(prior_a_b[1],prior_a_b[2]);
  // print("target a_b =",target());
  
  // b_b ~ uniform(prior_b_b[1],prior_b_b[2]);
  // print("target b_b =",target());
  
  scale_b ~ exponential(prior_scale_b);
  
  a ~ lognormal(mu_a,sigma_a);
  // print("target a =",target());
  
  // b ~ inv_gamma(a_b,b_b);
  // print("target b =",target());
  
  b ~ cauchy(0,scale_b);
  
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
    
    seg_fdL = to_array_1d(log(a[s] + b[s] * to_vector(seg_doseZ)));  // use fdL instead of log(fd)
    // seg_fdL = to_array_1d(log(a[s]+b[s]*to_vector(seg_dose)));     // use fdL instead of log(fd)
  
    target += Normal_Summary_lpdf(seg_ymeanL | seg_ysdL, seg_Nsub, seg_fdL, sigma);
    
    pos = pos + G[s];
  }
}

generated quantities{
  real log_lik[N];
  real fd[N];

  // Method 1: double index
  for(n in 1:N){
    fd[n] = a[Index[n]] + b[Index[n]] * doseZ[n];                          // simulate fd for double check
    // fd[n] = a[Index[n]] + b[Index[n]] * dose[n];
    log_lik[n] = Normal_Summary_lpdf({ymeanL[n]} | {ysdL[n]},{Nsub[n]},{log(fd[n])},sigma);
  }
}
