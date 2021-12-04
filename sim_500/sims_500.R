library(waveslim)
library(nloptr)
library(fracdiff)
library(pracma)
library(FKF)
library(Rsolnp)
setwd("./sim_500")

# Calculate Gegenbauer coeficients
ggbr.coef<-function(n,d,eta) {
  cf<-c(1,2*d*eta,2*d*(d+1)*eta^2-d)
  for (j in 3:(n-1)) cf[j+1]<-2*eta*((d-1)/j+1)*cf[j]-(2*(d-1)/j+1)*cf[j-1]
  return(cf)
}

initial.est<-function(y) {
  # Initial estimate for GARMA model with f,d,AR1
  # first identify the peak
  ssx<-spectrum(y,plot=F)
  startIdx<-as.integer(length(ssx$spec)/10)+1
  f<-ssx$freq[f_idx<-which.max(ssx$spec[startIdx:length(ssx$spec)])+startIdx-1]
  # next, estimate d
  m<-as.integer(length(y)^0.5)
  v<-log(1:m)
  v_mean<-sum(v)/m
  v<-v-v_mean
  denom<-2*sum(v^2)
  numer<-0
  for (j in 1:m) {
    idx1<-f_idx+j
    if (idx1>length(ssx$spec)) idx1<-length(ssx$spec)-(idx1-length(ssx$spec))
    idx2<-f_idx-j
    if (idx2<1) idx2<-abs(idx2)+1
    numer<-numer+v[j]*(log(ssx$spec[idx1])+log(ssx$spec[idx2]))
  }
  d<- (-0.5)*numer/denom
  if (d>0.5) d<-0.25
  if (d<0.0) d<-0.25
  x<-LMFilter(y,d,cos(2*pi*f))
  fit<-arima(x,order=c(1,0,0))
  phi<-fit$coef
  return(c(d=d,f=f-0.005,phi[1]))
}

# generic optimize function
  sim_optim<-function(par, objective, lower, upper, ...) {
    # check pars
    for (i in 1:length(par))
      if (par[i]<=lower[i]|par[i]>=upper[i]) par[i]<- (lower[i]+upper[i])/2   # if below lower bound then set to middle value
    
    res <- solnp(par, objective, LB=lower, UB=upper, control=list(tol=1e-12,trace=0), ...)
    res$par <- res$pars
    return(res)
  }
  

# Objective function for CSS method - GAR(1)
css.obj<-function(par,y) {
  theta1<-par[1]
  theta2<-par[2]
  theta3<-par[3]
  d<-theta1
  sigma2<-1 #exp(par[3])
  f<-theta2
  eta<-cos(2*pi*f)
  phi<-theta3
  n<-length(y)
  xx<-Toeplitz(ggbr.coef(n,d,eta),c(1,rep(0,n-1)))
  eps<-solve(xx)%*%(y-phi*c(0,y[1:(n-1)]))
  return(sum(eps^2)/(2*sigma2))
}

# Find CSS estimates for parameters for GAR(1)
css.est<-function(y) {
  fit<-sim_optim(initial.est(y), css.obj, lower=c(0,0,-1), upper=c(0.5,0.5,1), y=y)
  return(c(d=unname(fit$par[1]),f=unname(fit$par[2]),phi=unname(fit$par[3]),convergence=fit$convergence))
}


## Whittle Estimate
# Objective Function for Whittle method - GAR(1)
whittle.obj<-function(theta,ss) {
  fd<-theta[1]
  f<-theta[2]
  phi<-theta[3]
  u <- cos(2*pi*f)
  cos_2_pi_f <- cos(2.0*pi*ss$freq)
  mod_phi <- (1+phi^2-2*phi*cos_2_pi_f)
  spec_den_inv <- 2.0*pi * (4.0*((cos_2_pi_f-u)^2))^fd * mod_phi   # Inverse of spectral density
  
  spec_den_inv[is.infinite(spec_den_inv)] <- NA
  spec_den_inv[spec_den_inv<=0] <- .Machine$double.eps
  I_f <- ss$spec*spec_den_inv
  res <- sum(I_f-log(spec_den_inv),na.rm=TRUE)
  return(res)
}

# Find Whittle Estimates for parameters GAR(1)
whittle.est<-function(y) {
  ss<-spectrum(y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)
  fit<-sim_optim(initial.est(y), whittle.obj, lower=c(0,0,-1), upper=c(0.5,0.5,1), ss=ss)
  return(c(d=unname(fit$par[1]),f=unname(fit$par[2]),phi=unname(fit$par[3]),convergence=fit$convergence))
}

wll.ggbr.obj<-function(par,ss) {
  fd    <- par[1]
  f     <- par[2]
  phi   <- par[3]
  
  cos_2_pi_f <- cos(2.0*pi*ss$freq)
  u <- cos(2*pi*f)
  mod_phi <- (1+phi^2-2*phi*cos_2_pi_f)
  spec_den_inv <- 2.0*pi * (4.0*((cos_2_pi_f-u)^2))^fd * mod_phi   # Inverse of spectral density
  
  spec_den_inv[is.infinite(spec_den_inv)] <- NA
  spec_den_inv[spec_den_inv<=0] <- .Machine$double.eps
  I_f <- ss$spec*spec_den_inv
  res <- sum((log(I_f))^2,na.rm=TRUE)
  return(res)
}

wll.est<-function(y) {
  ss<-spectrum(y,plot=FALSE,detrend=FALSE,demean=FALSE,method='pgram',taper=0,fast=FALSE)
  lb<-c(0,0.05,-0.99)
  ub<-c(0.5,0.45,0.99)
  pars<-initial.est(y)
  fit <- sim_optim(pars, wll.ggbr.obj, lower=lb, upper=ub, ss=ss)
  return(c(d=unname(fit$par[1]),f=unname(fit$par[2]),phi=unname(fit$par[3]),convergence=fit$convergence))
}

#QMLE
# build state-space matricies for QML for GAR(1)
# With thanks to Hau Wu, U Syd.
qml<-function(d,u,phi) {
  si=as.vector(ggbr.coef(51,d,u))
  # add in an AR(1) factor
  si2<-as.vector(stats::filter(si,phi,method="recursive",sides=1))
  si_sq=si2 %*% t(si2)
  Tt <- rbind(cbind(rep(0,50),diag(50)),rep(0,51))
  Zt <- matrix(c(1, rep(0,50)), ncol = 51)
  ct <- matrix(0)
  dt <- matrix(0, nrow = 51)
  GGt <- matrix(0)
  H <- matrix(c(si2[1:51]), nrow = 51)
  HHt <- H %*% t(H)
  a0 <- c(rep(0,51))
  entry<-vector()
  for(j in 1:49){
    entry[j]=sum(diag(si_sq[-1:-j,-(52-j):-51]))
  }
  P0 <- toeplitz(c(sum(diag(si_sq)),entry[1:49],si_sq[1,51]))
  
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt, HHt = HHt))
}

# Objective Function for QML for GAR(1)
qml.obj <- function(theta, yt) {
  sp <- qml(theta[1], cos(2*pi*theta[2]),theta[3])
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
  return(-ans$logLik)
}

# Find QML Estimates for GAR(1)
qml.est<-function(y) {
  # fit<-lbfgs(initial.est(y), qml.obj, lower=c(0,0,-1), upper=c(0.5,0.5,1), control=list(maxeval=500), yt=rbind(y))
  fit<-sim_optim(initial.est(y), qml.obj, lower=c(0,0,-1), upper=c(0.5,0.5,1), yt=rbind(y))
  return(c(d=unname(fit$par[1]),f=unname(fit$par[2]),phi=unname(fit$par[3]),convergence=fit$convergence))
}

# take a g-process and remove the Gegenbauer component to leave the short memory process
LMFilter<-function(y,d,u) {
  n<-length(y)
  xx<-Toeplitz(ggbr.coef(n,d,u),c(1,rep(0,n-1)))
  x<-solve(xx)%*%y
  return(x)
}

# Estimate AR parameter given d, u
estPhi<-function(y,d,u) {
  x<-LMFilter(y,d,u)
  return(a<-arima(x,order=c(1,0,0),include.mean=FALSE)$coef[1])
}

# Method of Arteche & Robinson (2000) to estimate fractional differencing.
# use max of periodogram as per Yajima (1996) as Gegenbauer frequency
# from short memory component use standard R "arima" function to estimate phi.
ar1.est<-function(y) {
  # first identify the peak
  ssx<-spectrum(y,plot=F)
  startIdx<-as.integer(length(ssx$spec)/10)+1
  f<-ssx$freq[f_idx<-which.max(ssx$spec[startIdx:length(ssx$spec)])+startIdx-1]
  # next, estimate d
  m<-as.integer(length(y)^0.7)
  v<-log(1:m)
  v_mean<-sum(v)/m
  v<-v-v_mean
  denom<-2*sum(v^2)
  numer<-0
  for (j in 1:m) {
    idx1<-f_idx+j
    if (idx1>length(ssx$spec)) idx1<-length(ssx$spec)-(idx1-length(ssx$spec))
    idx2<-f_idx-j
    if (idx2<1) idx2<-abs(idx2)+1
    numer<-numer+v[j]*(log(ssx$spec[idx1])+log(ssx$spec[idx2]))
    #numer<-numer+v[j]*log(ssx$spec[idx1])
  }
  d<- (-0.5)*numer/denom
  phi<-estPhi(y,d,cos(2*pi*f))
  return(c(d=d,f=f,phi=unname(phi)))
}

# "local Whittle" Method of Arteche & Robinson (2000) to estimate fractional differencing.
# use max of periodogram as per Yajima (1996) as Gegenbauer frequency
# from short memory component use standard R "arima" function to estimate phi.
lw.est<-function(y) {
  # first identify the peak
  ssx<-spectrum(y,plot=F)
  startIdx<-as.integer(length(ssx$spec)/10)+1
  f<-ssx$freq[f_idx<-which.max(ssx$spec[startIdx:length(ssx$spec)])+startIdx-1]
  # next estimate d via local-whittle technique
  l<-0
  lw.r<-function(d,y,f_idx) {
    m<-as.integer(length(y)^0.7)
    sum1<-sum2<-sumlog<-0
    for (j in (l+1):m) {
      idx1<-f_idx+j
      if (idx1>length(ssx$spec)) idx1<-length(ssx$spec)-(idx1-length(ssx$spec))
      idx2<-f_idx-j
      if (idx2<1) idx2<-abs(idx2)+1
      sum1<-sum1+((2*pi*ssx$freq[j])^(2*d))*(ssx$spec[idx1])
      sum2<-sum2+((2*pi*ssx$freq[j])^(2*d))*(ssx$spec[idx2])
      sumlog<-sumlog+log(2*pi*ssx$freq[j])
    }
    lw.c1<-sum1/(m-l)
    lw.c2<-sum2/(m-l)
    return(log(lw.c1)-2*d/(m-l)*sumlog)
  }
  # now minimise lw.r
  fit<-sim_optim(c(0.1), lw.r, lower=c(0), upper=c(0.5), y=y, f_idx=f_idx)
  phi<-estPhi(y,fit$par[1],cos(2*pi*f))
  return(c(d=unname(fit$par[1]),f=f,phi=unname(phi)))
}

# Wavelets
# first we re-define the standard function to force it to ignore the low-frequency peak from the AR(1) component
# the standard version gets confused by this and reports the g-frequency to be low.
spp.mle.mod <- function(y, wf, J=log(length(y),2)-1, p=0.01, frac=1)
{
  ##
  ##  g a u s s g k . R  Adapitve Gauss-Kronrod
  ##
  
  # We re-define this funtion as the default can generate many errors around the unbounded peak in the spectral density.
  quadgk <- function(f, a, b, tol = .Machine$double.eps^0.5, ...) {
    stopifnot(is.numeric(a), length(a) == 1,
              is.numeric(b), length(b) == 1)
    eps <- .Machine$double.eps
    
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    
    if (a == b)     return(0)
    else if (a > b) return(-1 * quadgk(f, b, a, tol = tol))
    
    # Nodes and weights for Gauss-Kronrod (7, 15)
    n15 <- c(-0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
             -0.7415311855993944, -0.5860872354676911, -0.4058451513773972,
             -0.2077849550078985,  0.0,                 0.2077849550078985,
             0.4058451513773972,  0.5860872354676911,  0.7415311855993944,
             0.8648644233597691,  0.9491079123427585,  0.9914553711208126)
    n7  <- c(-0.9491079123427585, -0.7415311855993944, -0.4058451513773972,
             0.0,
             0.4058451513773972, 0.7415311855993944,  0.9491079123427585)
    
    w15 <- c(0.02293532201052922, 0.06309209262997855,  0.1047900103222502,
             0.1406532597155259,  0.1690047266392679,   0.1903505780647854,
             0.2044329400752989,  0.2094821410847278,   0.2044329400752989,
             0.1903505780647854,  0.1690047266392679,   0.1406532597155259,
             0.1047900103222502,  0.06309209262997855,  0.02293532201052922)
    w7  <- c(0.1294849661688697,  0.2797053914892767,   0.3818300505051189,
             0.4179591836734694,
             0.3818300505051189,  0.2797053914892767,   0.1294849661688697)
    
    .gkadpt <- function(f, a, b, tol = tol) {
      # use nodes and weights from the environment
      x15 <- 0.5 * ((b - a) * n15 + b + a)
      x7  <- 0.5 * ((b - a) * n7  + b + a)
      
      f7<-f(x7)
      if (length(which(is.infinite(f7)))) {
        n<-which(is.infinite(f7))
        for (nn in n) if (nn==7) f7[7]<-f7[6] else f7[nn]<-f7[nn+1]
      }
      f15<-f(x15)
      if (length(which(is.infinite(f15)))) {
        n<-which(is.infinite(f15))
        for (nn in n) if (nn==15) f15[15]<-f15[14] else f15[nn]<-f15[nn+1]
      }
      
      Q7  <- sum(w7  * f7)  * (b-a)/2
      Q15 <- sum(w15 * f15) * (b-a)/2
      
      if (!is.finite(Q7) || !is.finite(Q15)) {
        #warning("Infinite or NA function value encountered.")
        return(Q15)
      } else if (abs(Q15 - Q7) < tol) {
        return(Q15)
      } else if (abs(b-a) < 16*eps) {
        #warning("Minimum step size reached; singularity possible.")
        return(Q15)
      } # else
      
      Q2 <- .gkadpt(f, (a+b)/2, b, tol = tol)
      Q1 <- .gkadpt(f, a, (a+b)/2, tol = tol)
      
      return(Q1 + Q2)
    }
    
    # start the recursive procedure
    .gkadpt(f, a, b, tol = tol)
  }
  spp.ar1.sdf <- function(freq, ll)
    {return(abs(2 * (cos(2*pi*freq) - cos(2*pi*ll$fG)))^(-2*ll$d)/(1-2*ll$phi*cos(2*pi*freq)+ll$phi^2))}
  bandpass.spp.ar1 <- function(a, b, d, fG, phi) {
    # We change this to use Gauss-Konrad quadrature rather than the standard R "integrate" function
    # which had problems with the unbounded peak.
    if(fG > a && fG < b) {
      result1 <- quadgk(spp.ar1.sdf, a, fG, ll=list(d=d, fG=fG, phi=phi))
      result2 <- quadgk(spp.ar1.sdf, fG, b, ll=list(d=d, fG=fG, phi=phi))
    }
    else {
      result1 <- quadgk(spp.ar1.sdf, a, b, ll=list(d=d, fG=fG, phi=phi))
      result2 <- 0
    }
    return(2*(result1 + result2))
  }
  sppLL <- function(x, y) {
    delta <- x[1]
    fG <- x[2]
    phi <- x[3]
    y.dwpt <- y[[1]]
    y.basis <- y[[2]]
    n <- y[[3]]
    J <- y[[4]]
    
    ## Establish the limits of integration for the band-pass variances
    a <- unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 2^(rep(1:J, 2^(1:J))) / 2
    b <- unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) / 2^(rep(1:J, 2^(1:J))) / 2
    
    ## Define some useful parameters for the wavelet packet tree
    length.jn <- n / rep(2^(1:J), 2^(1:J))
    scale.jn <- rep(2^(1:J+1), 2^(1:J))
    
    ## Initialize various parameters for the reduced LL
    Basis <- (1:length(y.basis))[y.basis]
    bp.var <- numeric(length(Basis))
    delta.n <- 100
    
    ## Compute the band-pass variances according to \delta and f_G
    omega.diag <- NULL
    for(i in 1:sum(y.basis)) {
      jn <- Basis[i]
      bp.var[i] <- bandpass.spp.ar1(a[jn], b[jn], delta, fG, phi)
      # If we get infinities, replace with a large number
      if (is.infinite(bp.var[i])) bp.var[i]<- 1e15
      omega.diag <- c(omega.diag, scale.jn[jn] * rep(bp.var[i], length.jn[jn]))
    }

    ## Compute reduced log-likelihood 
    rLL <- n * log(1/n * sum(y.dwpt^2 / omega.diag, na.rm=TRUE)) + sum(length.jn[y.basis] * log(scale.jn[y.basis] * bp.var))
    rLL
  }
  
  n <- length(y)
  x0 <- numeric(3)
  
  ## Perform discrete wavelet packet transform (DWPT) on Y
  y.dwpt <- dwpt(y, wf, n.levels=J)
  n <- length(y)
  if(frac < 1) {
    for(i in 1:length(y.dwpt)) {
      vec <- y.dwpt[[i]]
      ni <- length(vec)
      j <- rep(1:J, 2^(1:J))[i]
      vec[trunc(frac * n/2^j):ni] <- NA
      y.dwpt[[i]] <- vec
    }
  }
  y.basis <- as.logical(ortho.basis(portmanteau.test(y.dwpt, p)))
  y.dwpt <- as.matrix(unlist(y.dwpt[y.basis]))
  
  ## Compute initial estimate of the Gegenbauer frequency
  y.per <- per(y - mean(y))
  y.per1<-y.per
  y.per1[1:as.integer(n/10)]<-0
  x0[2] <- (0:(n/2)/n)[max(y.per1) == y.per1]
  
  ## Compute initial estimate of the fractional difference parameter
  muJ <- (unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0, by=1)) / 2^(rep(1:J, 2^(1:J))) + 
            unlist(apply(matrix(2^(1:J)), 1, seq, from=1, by=1)) / 2^(rep(1:J, 2^(1:J)))) / 4
  y.modwpt <- modwpt(y, wf=wf, n.levels=J)
  y.varJ <- rep(2^(1:J), 2^(1:J)) * unlist(lapply(y.modwpt, FUN=function(x)sum(x*x,na.rm=TRUE)/length(x[!is.na(x)])))
    
  lb <- c(0.001,0.050,-1)
  ub <- c(0.499,0.45,1)
  x0[1] <- min(-0.5 * lsfit(log(abs(muJ[y.basis] - x0[2])), log(y.varJ[y.basis]))$coef[2], 0.49)
  if (x0[1]<=0) x0[1]<- (lb[1]+ub[1])/2   # if below lower bound then set to middle value
  
  x0[3] <- 0    # Initial estimate for phi
  result <- sim_optim(x0, sppLL, lower=lb, upper=ub, y=list(y.dwpt, y.basis, n, J))
  return(result)
}

# Find Wavelet Estimates for GAR(1)
wavelet.est<-function(y,filter){
  flag<-TRUE
  fit<-list(par=c(d=NA,f=NA,phi=NA),convergence=1)
  tryCatch(fit<-spp.mle.mod(y, filter), error = function(e) flag<-FALSE)
  return(c(d=unname(fit$par[1]),f=unname(fit$par[2]),phi=unname(fit$par[3]),convergence=fit$convergence))
}

# Generate Realizations
# Theoretical Spectral Density
GgbrAR1SpecDen<-function(x,ll) {return (cos(ll$k*x)*(4*(cos(x)-cos(2*pi*ll$f))^2)^(-ll$d)/(1-2*ll$phi*cos(x)+ll$phi^2));}

# Get acf of theoretical spectral density using numerical integration
ggbr.sim<-function(n,d,f,phi,sigma2) {
  # get acf by numerical integration
  g<-rep(0,n)
  for (k in 1:n) g[k]<-quadgk(GgbrAR1SpecDen, 0, 2*pi, ll=list(d=d,f=f,phi=phi,k=k-1,sigma2=sigma2))
  g<-unlist(g)
  g<-g/max(g)
  return(g)
}

# data frame to store the simulated processes
synProcess_22_40_80<-data.frame(obs=1:500)
#periodograms for each of the above
sd_22_40_80<-data.frame(obs=1:250)
i<-1
while (i<=1000) {
  # from "waveslim" package, use method of Hosking (1984) to generate a realization
  y2<-hosking.sim(1000,g)
  synProcess_22_40_80[,paste0("y",i)]<-y2[501:1000]
  ss<-spectrum(y2[501:1000],plot=F,detrend=F)
  # check the peak in the spectral density is high enough; if not throw this one away.
  if (max(ss$spec[50:250])>15) {
    if (i==1) sd_22_40_80$freq<-ss$freq
    sd_22_40_80[,paste0("spec",i)]<-ss$spec
    i<-i+1
  }
}
  
synProcess_22_40_80<-readRDS("synProcess_22_40_80.RDS")
results_22_40_80<-data.frame(css_d=numeric(1000),css_f=numeric(1000),css_phi=numeric(1000),css_t=numeric(1000), css_conv=numeric(1000),
                             whittle_d=numeric(1000),whittle_f=numeric(1000),whittle_phi=numeric(1000),
                             whittle_t=numeric(1000),whittle_conv = numeric(1000),
                             wll_d=numeric(1000),wll_f=numeric(1000),wll_phi=numeric(1000),wll_t=numeric(1000),wll_conv = numeric(1000),
                             qml_d=numeric(1000),qml_f=numeric(1000),qml_phi=numeric(1000), qml_t=numeric(1000),qml_conv = numeric(1000),
                             mb8_d=numeric(1000),mb8_f=numeric(1000),mb8_phi=numeric(1000), mb8_t=numeric(1000),mb8_conv = numeric(1000),
                             lw_d=numeric(1000),lw_f=numeric(1000),lw_phi=numeric(1000),lw_t=numeric(1000),lw_conv = numeric(1000),
                             ar1_d=numeric(1000),ar1_f=numeric(1000),ar1_phi=numeric(1000),ar1_t=numeric(1000),ar1_conv = numeric(1000)
)

for (i in 1:1000) {
  y<-synProcess_22_40_80[,paste0("y",i)]
  tt<-system.time(css<-css.est(y))
  results_22_40_80$css_d[i]<-css[['d']]
  results_22_40_80$css_f[i]<-css[['f']]
  results_22_40_80$css_phi[i]<-css[['phi']]
  results_22_40_80$css_conv[i]<-css[['convergence']]
  results_22_40_80$css_t[i]<-tt[1]
  tt<-system.time(qml1<-qml.est(y))
  results_22_40_80$qml_d[i]<-qml1[['d']]
  results_22_40_80$qml_f[i]<-qml1[['f']]
  results_22_40_80$qml_phi[i]<-qml1[['phi']]
  results_22_40_80$qml_conv[i]<-qml1[['convergence']]
  results_22_40_80$qml_t[i]<-tt[1]
  tt<-system.time(whittle<-whittle.est(y))
  results_22_40_80$whittle_d[i]<-whittle[['d']]
  results_22_40_80$whittle_f[i]<-whittle[['f']]
  results_22_40_80$whittle_phi[i]<-whittle[['phi']]
  results_22_40_80$whittle_conv[i]<-whittle[['convergence']]
  results_22_40_80$whittle_t[i]<-tt[1]
  tt<-system.time(wll<-wll.est(y))
  results_22_40_80$wll_d[i]<-wll[['d']]
  results_22_40_80$wll_f[i]<-wll[['f']]
  results_22_40_80$wll_phi[i]<-wll[['phi']]
  results_22_40_80$wll_conv[i]<-wll[['convergence']]
  results_22_40_80$wll_t[i]<-tt[1]
  yd<-c(y,rep(0,12))
  tt<-system.time(mb8<-wavelet.est(yd, "mb8"))
  results_22_40_80$mb8_d[i]<-mb8[['d']]
  results_22_40_80$mb8_f[i]<-mb8[['f']]
  results_22_40_80$mb8_phi[i]<-mb8[['phi']]
  results_22_40_80$mb8_conv[i]<-mb8[['convergence']]
  results_22_40_80$mb8_t[i]<-tt[1]
  tt<-system.time(x1<-ar1.est(y))
  results_22_40_80$ar1_d[i]<-x1[['d']]
  results_22_40_80$ar1_f[i]<-x1[['f']]
  results_22_40_80$ar1_phi[i]<-x1[['phi']]
  results_22_40_80$ar1_t[i]<-tt[1]
  results_22_40_80$ar1_conv[i]<-0
  tt<-system.time(x2<-lw.est(y))
  results_22_40_80$lw_d[i]<-x2[['d']]
  results_22_40_80$lw_f[i]<-x2[['f']]
  results_22_40_80$lw_phi[i]<-x2[['phi']]
  results_22_40_80$lw_t[i]<-tt[1]
  results_22_40_80$lw_conv[i]<-0
  if ((i%%25)==0) {
    cat(paste(i,"\n"))
    saveRDS(results_22_40_80,"results_22_40_80.RDS")
  }
}
summary(results_22_40_80)

saveRDS(results_22_40_80,"results_22_40_80.RDS")
#   results_22_40_80<-readRDS("results_22_40_80.RDS")

mse <- function(x, x_true) {return(mean((x-x_true)^2))}

methods<-c("css","whittle","wll","qml","mb8","ar1","lw")
mse<-function(x,true_value) {sum((x-true_value)^2,na.rm=T)/length(!is.na(x))}
res<-sapply(methods,function(x) sprintf("%10s %9.6f %9.6f",x,mean(results_22_40_80[,paste0(x,"_d")],na.rm=T)-0.22,mse(results_22_40_80[,paste0(x,"_d")],0.22)))
cat(c("     d est      mean       mse",res),sep="\n")
res<-sapply(methods,function(x) sprintf("%10s %9.6f %9.6f",x,mean(results_22_40_80[,paste0(x,"_f")],na.rm=T)-0.4,mse(results_22_40_80[,paste0(x,"_f")],0.40)))
cat(c("     f est      mean       mse",res),sep="\n")
res<-sapply(methods,function(x) sprintf("%10s %9.6f %9.6f",x,mean(results_22_40_80[,paste0(x,"_phi")],na.rm=T)-0.8,mse(results_22_40_80[,paste0(x,"_phi")],0.80)))
cat(c("   phi est      mean       mse",res),sep="\n")
res<-sapply(methods,function(x) sprintf("%10s %9.6f %9.6f",x,mean(results_22_40_80[,paste0(x,"_t")],na.rm=T),mse(results_22_40_80[,paste0(x,"_t")],0)))
cat(c("  time est      mean       mse",res),sep="\n")

length(which(is.na(results_22_40_80$mb8_d)))
