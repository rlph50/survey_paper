data {
  int<lower=1> N;            // num observations
  real y[N];                 // observed outputs
}
parameters {
  real<lower=-1,upper=1> phi;                  // autoregression coeff
  real<lower=0.0,upper=0.5> d;
  real<lower=0.0,upper=0.5> f;
}
model {
  vector[N] nu;
  vector[N] err;
  real g[N+1];
  real u;
  
  u=cos(2.0*pi()*f);
  
  g[1] = 1.0;                // this and next few lines to derive the Ggbr Coefficients.
  g[2] = 2.0*u*d;
  for (j in 3:(N+1)) 
  {
    real rj1;
    real temp;
    
    rj1 = j-1;
    temp = (d-1.0)/rj1;
    g[j]=(2.0*u*(rj1+d-1.0)*g[j-1]-(rj1+2.0*d-2.0)*g[j-2])/rj1;
  }

  nu[1] = 0.0;
  err[1] = y[1] - nu[1];
  for (t in 2:N) {
    real sum_g;
    sum_g = 0.0;
    for (j in 1:(t-1)) sum_g += g[j+1]*err[t-j];

    nu[t] = phi * y[t-1] + sum_g;
    err[t] = y[t] - nu[t];
  }
  phi ~ normal(0.0, 0.5);
  d   ~ normal(0.0, 0.25);
  f   ~ lognormal(0.0, 0.2);
  err ~ normal(0.0, 1.0);
}
