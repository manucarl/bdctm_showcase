##
## Script name: patent_hurdle.R
##
## Purpose of script: estimates BDCTM_hurdle and plots quantile residuals
##
## Author: Manuel Carlan
##
## Date Created: 2021-10-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

library(tidyverse)
set.seed(42)

# rm(list=ls())
source("code/nuts/adnuts_helper.R")


packages <- c("MASS","BayesX", "Rcpp", "RcppArmadillo", "splines",  "Matrix",  "sdPrior",
              "R2BayesX", "tidyverse",  "numDeriv", "tictoc",  "scales",  "scales", "coda",
              "mlt", "doParallel", "scam", "mvtnorm", "MCMCpack", "mgcv", "mcmcplots", "cotram", "bamlss", "RhpcBLASctl")# library(doParallel)

lapply(packages, library, character.only = TRUE) %>%  invisible()


data <-read.table("raw_data/patent.raw", header=T)
sourceCpp("code/rcpp/log_disc_xx.cpp")

source("code/nuts/nuts_patent_hurdle.R")


# scale continuous variables, create and scale dummys
data$nclaimss <- scale(data$nclaims)
data$years <- scale(data$year)
data$ncountrys <- scale(data$ncountry)

data$opps <- ifelse(data$opp == 0, 1, -1)
data$biopharms <- ifelse(data$biopharm  == 0, 1, -1)
data$ustwins  <- ifelse(data$ustwin  == 0, 1, -1)
data$patuss  <- ifelse(data$patus  == 0, 1, -1)
data$patgsgrs     <- ifelse(data$patgsgr     == 0, 1, -1)

data$opps <- scale(data$opp)
data$biopharms <- scale(data$biopharm)
data$ustwins  <- scale(data$ustwin)
data$patuss  <- scale(data$patus)
data$patgsgrs     <- scale(data$patgsgr)


# apply log(y+1) transformation or not
log_plus <-T
plus_diag <- 0
m <- 6


# indices of exponentiated params
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

# lagged trafo basis
B0l <-Predict.matrix(sm0, data=data.frame(y=yfun((ystore-1))))

# linear shift and hurdle effects
x_lin <- data %>% dplyr::select(ncountrys, years, nclaimss)

p0 <- length(x_lin)

# precision matrix
Klist <- list(K0)
S <- as.matrix(bdiag(K0, matrix(0,p0, p0)))


trafo <- "hurdle"
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

# adapted bases used in posterior calculation
Xzero <- X[ystore ==0,  ] %>% as.matrix
Xrest <- X[ystore >0,  ] %>% as.matrix
Xlrest <- Xl[ystore >0,  ] %>% as.matrix


# starting values
start <- rep(0,p)

hyperparams <- list(a=1, b=0.001)
method <- "logit"

# list including model components passed to rcpp
xx <- list(Xzero=Xzero, Xrest=Xrest, Xlrest=Xlrest,S= S, Klist=Klist,
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

its <- 2000
warmup <- burnin <- its/2
settings <- list(max_treedepth=12, adapt_delta=0.95)

fit <-NUTS(n_iter=its, xx=xx, f=posterior, gr=gradf, ll=ll, start=start, warmup=warmup, nuts_settings=settings)

# Final acceptance ratio=0.96, and target=0.95
# Final step size=0.053; after 1000 warmup iterations
# Elapsed Time: 353.3 seconds (Warmup)
# Elapsed Time: 299.0 seconds (Sampling)
# Elapsed Time: 652.3 seconds (Total)

# save(fit, file="processed_data/patent_hurdle.RData")
load("processed_data/patent_hurdle.RData")

beta_samples <- fit$beta
effectiveSize(beta_samples %>% as.mcmc)
# var1      var2      var3      var4      var5      var6      var7      var8      var9     var10     var11     var12     var13     var14     var15 
# 464.9498  421.4687  468.3144  811.4905  831.8899  839.5818 1144.8623  758.8530 1482.0083  996.4343 1697.7438  494.8198 1359.5784 1834.5187 1549.5993 


# load("processed_data/patent_hurdle.RData")

beta <- fit$beta[1001:2000, 1:ncol(X)] %>% colMeans
bt <- beta
bt[ind_exp] <- exp(bt[ind_exp])

# quantile residuals

h0_ind <- m+3
h0l <- h0 <- rep(0, nrow(X))
h0[ystore == 0] <- bt[h0_ind]
h0l[(ystore-1) == 0] <- bt[h0_ind]



lq <- B0l %*% bt[1:(m+2)]  + h0l
lq[ystore==0] <- 0
lq <- plogis(lq)

rq <- B0 %*% bt[1:(m+2)]  + h0
rq[ystore==0] <- 0
rq <- plogis(rq)

u <-runif(nrow(X), min=lq, max=rq)



qres <- rnorm(u)
# save(qres, file="processed_data/qr_bdctm_hurdle.RData")
ggpubr::ggqqplot(data.frame(qres=qres), x = "qres", title="BDCTM_hurdle")
