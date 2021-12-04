library(rstan)
library(rstantools)
library(bayesplot)
options(mc.cores=parallel::detectCores())
setwd("./sim_1000")

synProcess_22_40_80<-readRDS("synProcess_22_40_80.RDS")
res <- data.frame(d_mean=numeric(1000), f_mean=numeric(1000), phi_mean=numeric(1000),
                  d_neff=numeric(1000), f_neff=numeric(1000), phi_neff=numeric(1000),
                  d_rhat=numeric(1000), f_rhat=numeric(1000), phi_rhat=numeric(1000), t=numeric(1000))
#res <- readRDS('bayes_est_1000.rds')
    
for (i in 1:1000) {
  cat(sprintf("sim %d\n",i))
  y<-synProcess_22_40_80[,paste0("y",i)]
  data<-list(N=length(y),y=y)
  tt<-system.time(fit1 <- stan(
    file = "bayes_est.stan", # Stan program
    data = data,              # named list of data
    chains = 4,               # number of Markov chains
    warmup = 500,             # number of warmup iterations per chain
    iter = 2000,              # total number of iterations per chain
    cores = 4,                # number of cores (could use one per chain)
    seed=1798213,
    control=list(max_treedepth=12,adapt_delta=0.9),
    #control=list(adapt_delta=0.9),
    #control=list(max_treedepth=14,adapt_delta=0.99),
    refresh = 0             # no progress shown
  ))
  xx<-summary(fit1)$summary[,c(1,9,10)]
  print(xx)
  res$d_mean[i]<-xx['d','mean']
  res$f_mean[i]<-xx['f','mean']
  res$phi_mean[i]<-xx['phi','mean']
  res$d_neff[i] <- xx['d','n_eff']
  res$f_neff[i] <- xx['f','n_eff']
  res$phi_neff[i] <- xx['phi','n_eff']
  res$d_rhat[i] <- xx['d','Rhat']
  res$f_rhat[i] <- xx['f','Rhat']
  res$phi_rhat[i] <- xx['phi','Rhat']
  res$t[i]<-tt['elapsed']
  saveRDS(res,'bayes_est_1000.rds')
}


mse <- function(x, x_true) {return(mean((x-x_true)^2))}
res %>%
  filter(abs(d_rhat-1.0)<0.01) %>%
  summarise(d_m=mean(d_mean),d_mse=mse(d_mean,0.22),f_m=mean(f_mean),f_mse=mse(f_mean,0.40),phi_m=mean(phi_mean), phi_mse=mse(phi_mean,0.80), n = n())

saveRDS(res,'bayes_est_1000.rds')
    

i<-1
y<-synProcess_22_40_80[,paste0("y",i)]
data<-list(N=length(y),y=y)
tt<-system.time(fit1 <- stan(
  file = "sims_gar1b1.stan", # Stan program
  data = data,            # named list of data
  chains = 4,             # number of Markov chains
  warmup = 500,          # number of warmup iterations per chain
  iter = 2000,           # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  seed=1798213,
  control=list(max_treedepth=12,adapt_delta=0.9)
  #control=list(max_treedepth=14,adapt_delta=0.99),
  #refresh = 0             # no progress shown
))
fit1

traceplot(fit1) + ggtitle('stan Traceplot for Series 1')
pairs(fit1,main='stan pairs plot for series 1')

x<-extract(fit1)
mean(x$d)
mean(x$d[x$d>0.1])
mean(x$f)
mean(x$f[x$f>0.1])
mean(x$phi)
mean(x$phi[x$f>0.1])

