##
## Script name: patent_lin.R
##
## Purpose of script: estimates BDCTM_lin and plots rootograms and quantile residuals
##
## Author: Manuel Carlan
##
## Date Created: 2021-10-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

set.seed(42)

library(tidyverse)

# rm(list=ls())

source("code/nuts/adnuts_helper.R")


packages <- c("MASS","BayesX", "Rcpp", "RcppArmadillo", "splines",  "Matrix",  "sdPrior",
              "R2BayesX", "tidyverse",  "numDeriv", "tictoc",  "scales",  "scales", "coda",
              "mlt", "doParallel", "scam", "mvtnorm", "MCMCpack", "mgcv", "mcmcplots", "cotram", "bamlss", "RhpcBLASctl")# library(doParallel)

lapply(packages, library, character.only = TRUE) %>%  invisible()

data <-read.table("raw_data/patent.raw", header=T)
sourceCpp("code/rcpp/log_disc_xx2.cpp")

source("code/nuts/nuts_patent.R")


# scale continuous variables, create and scale dummys
data$nclaimss <- scale(data$nclaims)
data$years <- scale(data$year)
data$ncountrys <- scale(data$ncountry)

data$opps <- ifelse(data$opp == 0, 1, -1)
data$biopharms <- ifelse(data$biopharm  == 0, 1, -1)
data$ustwins  <- ifelse(data$ustwin  == 0, 1, -1)
data$patuss  <- ifelse(data$patus  == 0, 1, -1)
data$patgsgrs     <- ifelse(data$patgsgr     == 0, 1, -1)

data$opps <- scale(data$opps)
data$biopharms <- scale(data$biopharms)
data$ustwins  <- scale(data$ustwins)
data$patuss  <- scale(data$patuss)
data$patgsgrs     <- scale(data$patgsgrs)


# apply log(y+1) transformation or not
log_plus <-T
m <- 6

# indices of exponentiated coefs
ind_exp=2:(m+2)

ystore <- y <- data$ncit

# construct monotonic trafo
sspec <- mgcv::s(y,k=m+2,bs="ps")
sspec$mono <- 1

yfun <- function(x) x
if(log_plus) yfun <- function(x) log(x+1)

y <- yfun(ystore)
sm0 <- smoothCon(sspec,data=data.frame(y=y), scale.penalty=F, absorb.cons=F ,knots=NULL)[[1]]
B0 <- sm0$X
K0 <- sm0$S[[1]]

# construct lagged basis
B0l <-Predict.matrix(sm0, data=data.frame(y=yfun((ystore-1))))

# linear shift effects
x_lin <- data %>% dplyr::select(ncountrys, years, nclaimss)

p0 <- length(x_lin)

# precision matrices
Klist <- list(K0)
S <- as.matrix(bdiag(K0, matrix(0,p0, p0)))


trafo <- "tram"
X <- cbind(B0, -x_lin) 
Xl <- cbind(B0l, -x_lin)
if(trafo=="tram"){
  
  X <- cbind(B0, -x_lin)
  Xl <- cbind(B0l, -x_lin)

} else if(trafo =="hurdle"){
  p0 <- 2*p0 +1
  S <- as.matrix(bdiag(K0, matrix(0,p0, p0)))
  y_ind <- ifelse(ystore==0, 0, 1) #%>% scale
  y_ind_lag <- ifelse(ystore-1==0, 0, 1)# %>% scale
  
  X <- as.matrix(data.frame( B0=B0, -x_lin, y_ind , -y_ind*x_lin))
  Xl <- as.matrix(data.frame( B0l=B0l, -x_lin, y_ind_lag , -y_ind_lag*x_lin))
  
}

p <-ncol(X)

# design matrices used in posterior calculation
Xzero <- X[ystore ==0,  ] %>% as.matrix
Xrest <- X[ystore >0,  ] %>% as.matrix
Xlrest <- Xl[ystore >0,  ] %>% as.matrix

# starting values
start <- rep(0,p)

hyperparams <- list(a=1, b=0.001)
method <- "logit"

Xl[ystore==0,] = 0

X <- X %>% as.matrix
Xl <- Xl %>% as.matrix
# list including model components passed to rcpp

Sb = 
xx <- list(Xzero=Xzero, Xrest=Xrest, Xlrest=Xlrest,S= S, Klist=Klist, zero_ind = which(ystore ==0)-1, X=X, Xl=Xl,
           ind_exp=ind_exp-1, p=p, p0=p0, pnlin=length(Klist)-1, hyperparams=hyperparams)

if(method=="logit"){
  cdfun <- plogis
  posterior = posterior_logit
  gradf = gradf_logit
  ll <- ll_logit
} else if(method=="probit"){
  cdfun <- pnorm
  posterior = posterior_probit
  gradf = gradf_probit
  ll <- ll_probit
} else if(method == "cloglog"){
  cdfun <- pevd
  posterior = posterior_cll
  gradf = gradf_cll
  # start <- runif(length(start), min=-2, max=-1)
  start <- c(-10, rep(-1,p-1))
  ll <- ll_cll
}


# NUTS settings
settings <- list(max_treedepth=12, adapt_delta=0.8)
its <- 2000
warmup <- burnin <- its/2


fit <-NUTS(n_iter=its, xx=xx, f=posterior, gr=gradf, ll=ll, start=start, warmup=warmup, nuts_settings=settings)

# save(fit, file="processed_data/patent_reg.RData")
load("processed_data/patent_reg.RData")

beta_samples <- fit$beta[1001:2000, 1:ncol(X)] 

# Elapsed Time: 76.8 seconds (Warmup)
# Elapsed Time: 66.1 seconds (Sampling)
# Elapsed Time: 142.9 seconds (Total)

effectiveSize(beta_samples %>% as.mcmc)
# var1     var2     var3     var4     var5     var6     var7     var8     var9    var10    var11 
# 332.5035 302.6217 426.4560 498.3336 671.2074 789.1774 818.1849 760.1287 704.1059 880.5740 844.1630 

beta <- beta_samples %>% colMeans
bt <- beta
bt[ind_exp] <- exp(bt[ind_exp])


X <- as.matrix(X)
Xl <- as.matrix(Xl)

# Rootograms
counts <- 15
ypred <- 0:counts 

y <- data$ncit

B0pred <-Predict.matrix(sm0, data=data.frame(y=ypred %>% yfun))
B0lpred <-Predict.matrix(sm0, data=data.frame(y=(ypred-1) %>% yfun))


probs <-c(cdfun(B0pred[1,]%*%bt[1:(m+2)]  ),
          cdfun(B0pred[2:(counts+1),]%*%bt[1:(m+2)])- cdfun(B0lpred[2:(counts+1),]%*%bt[1:(m+2)])) 


library(countreg)
m_p <-    glm(ncit ~ opp+ biopharm  +patus+ patgsgr + ncountryn + yearn + nclaimsn, data=data, family="poisson")
m_nb <- glm.nb(ncit  ~ opp + biopharm + ustwin + patus +patgsgr, data = data)

par(mfrow=c(1,3))
rootogram(table(y)[1:(counts+1)], probs*length(y), style="hanging", main="Regular BDCTM")
rootogram(m_p, main="Poisson")
rootogram(m_nb, main="Negative Binomial")



# Quantile Residuals
lq <- B0l %*% bt[1:(m+2)]
lq[ystore==0] <- 0
lq <- plogis(lq)

rq <- B0 %*% bt[1:(m+2)]
rq[ystore==0] <- 0
rq <- plogis(rq)

u <-runif(nrow(X), min=lq, max=rq)



qres <- rnorm(u)
# save(qres, file="processed_data/qr_bdctm_reg.RData")

ggpubr::ggqqplot(data.frame(qres=qres+1), x = "qres", title="Regular BDCTM")

