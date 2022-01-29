##
## Script name: patent_nonlin.R
##
## Purpose of script: estimates BDCTM_nl and plots nonlinear effects and quantile residuals
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

# functions from package adnuts
source("code/nuts/adnuts_helper.R")


#needed packages
packages <- c("MASS","BayesX", "Rcpp", "RcppArmadillo", "splines",  "Matrix",  "sdPrior",
              "R2BayesX", "tidyverse",  "numDeriv", "tictoc",  "scales",  "scales", "coda",
              "mlt", "doParallel", "scam", "mvtnorm", "MCMCpack", "mgcv", "mcmcplots", "cotram", "bamlss", "RhpcBLASctl")# library(doParallel)

lapply(packages, library, character.only = TRUE) %>%  invisible()


data <-read.table("raw_data/patent.raw", header=T)
sourceCpp("code/rcpp/log_disc_xx.cpp")

source("code/nuts/nuts_patent.R")


# scale continuous variables, create and scale dummys
data$nclaimss <- scale(data$nclaims)
data$years <- scale(data$year)
data$ncountrys <- scale(data$ncountry)

data$opps <- ifelse(data$opp == 0, 1, -1) %>% scale
data$biopharms <- ifelse(data$biopharm  == 0, 1, -1) %>% scale
data$ustwins  <- ifelse(data$ustwin  == 0, 1, -1)%>% scale
data$patuss  <- ifelse(data$patus  == 0, 1, -1)%>% scale
data$patgsgrs     <- ifelse(data$patgsgr     == 0, 1, -1)%>% scale




# apply log(y+1) transformation or not
log_plus <- T

# number of knots for monotonic spline
m <- 6

# indices of exponentiated parameters
ind_exp=2:(m+2)
ystore <- y <- data$ncit


# construct trafo
sspec <- s(y,k=m+2,bs="ps")
sspec$mono <- 1
yfun <- function(x) x
if(log_plus) yfun <- function(x) log(x+1)

sm0 <- smoothCon(sspec,data=data.frame(y=yfun(ystore)), scale.penalty=F, absorb.cons=F ,knots=NULL)[[1]]
B0 <- sm0$X
K0 <- sm0$S[[1]]

# lagged basis
B0l <-PredictMat(sm0, data=data.frame(y=yfun((ystore-1))))

# list of precision matrices
Klist <- list(K0)


# construct nonlinear bases
sm1 <- smoothCon(s(ncountryn, k=10+2,bs="ps", m=2), absorb.cons=T, data=data, sparse.cons = 1, scale.penalty=F)
B1 <- sm1[[1]]$X
K1 <- sm1[[1]]$S[[1]]

sm2 <- smoothCon(s(yearn,k=10+2,bs="ps", m=2), absorb.cons=T, data=data, sparse.cons = 1, scale.penalty=F)
B2 <- sm2[[1]]$X
K2 <- sm2[[1]]$S[[1]]

sm3 <- smoothCon(s(nclaimsn,k=10+2,bs="ps", m=2), absorb.cons=T, data=data, sparse.cons =1, scale.penalty=F)
B3 <- sm3[[1]]$X
K3 <- sm3[[1]]$S[[1]]


x <- cbind(B1, B2, B3)
S <- as.matrix(bdiag(K0, K1, K2, K3))



trafo <- "tram"

if(trafo=="tram"){
  
  
  X <- cbind(B0, -x)
  Xl <- cbind(B0l, -x)

} else if(trafo =="hurdle"){
  S2 <- as.matrix(bdiag(0,S[-(1:(m+2)), -(1:(m+2))]))
  
  X <- as.matrix(data.frame( B0=B0, -x, ystore==0 , -(ystore==0)*x))
  Xl <- as.matrix(data.frame( B0l=B0l, -x, (ystore-1)==0 , -((ystore-1)==0)*x))
  
  S <- as.matrix(bdiag(S, S2))
  Klist <- c(Klist, Klist[-1])
}

p <-ncol(X)
p0 <- 0

Xzero <- X[ystore ==0,  ]
Xrest <- X[ystore >0,  ]
Xlrest <- Xl[ystore >0,  ]


start <- rep(0,p)
hyperparams=list(a=1, b=0.001)

Klist <- list(K0, K1, K2, K3)



method <- "logit"
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

# list including model components passed to rcpp
Xl[ystore==0,] = 0

X <- X %>% as.matrix
Xl <- Xl %>% as.matrix
# list including model components passed to rcpp
xx <- list(Xzero=Xzero, Xrest=Xrest, Xlrest=Xlrest,S= S, Klist=Klist, zero_ind = which(ystore ==0)-1, X=X, Xl=Xl,
           ind_exp=ind_exp-1, p=p, p0=p0, pnlin=length(Klist)-1, hyperparams=hyperparams)



# NUTS settings
settings <- list(max_treedepth=10, adapt_delta=0.99)
its <- 2000
warmup <- burnin <- its/2

fit <-NUTS(n_iter=its, xx=xx, f=posterior, gr=gradf, ll=ll, start=start, warmup=warmup, nuts_settings=settings)


load("processed_data/patent_nl.RData")

# Final acceptance ratio=0.98, and target=0.99
# Final step size=0.022; after 1000 warmup iterations
# Elapsed Time: 2811.9 seconds (Warmup)
# Elapsed Time: 2077.4 seconds (Sampling)
# Elapsed Time: 4889.3 seconds (Total)
beta_samples <- fit$beta

effectiveSize(beta_samples %>% as.mcmc)
# > effectiveSize(beta_samples %>% as.mcmc)
# var1      var2      var3      var4      var5      var6      var7      var8      var9     var10     var11     var12     var13     var14     var15     var16     var17 
# 852.4153 1130.5178 1290.3243 1608.3513  954.5306 1312.7861 1245.9771 1669.6004  885.5888  948.0742  907.8785  651.5658 1102.3253 1086.3215  229.1450  465.6823  329.9278 
# var18     var19     var20     var21     var22     var23     var24     var25     var26     var27     var28     var29     var30     var31     var32     var33     var34 
# 271.4043  363.5001 1063.7991 1028.3985  749.3656  477.9426  649.7952  819.4769  414.4136  274.1306  692.8603  530.1630  660.1992  405.3505  374.7174  481.0817  720.5784 
# var35     var36     var37     var38     var39     var40     var41 
# 521.9760  443.5388  434.3296  355.0840  351.7907  391.7712  427.9160 

beta <- fit$beta[1001:2000, 1:ncol(X)] %>% colMeans


bt <- beta
bt[ind_exp] <- exp(bt[ind_exp])


X <- as.matrix(X)
Xl <- as.matrix(Xl)

# Quantile Residuals

u <- runif(nrow(X), plogis(ifelse(B0l %*% bt[1:(m+2)] %>%  is.nan, 0, B0l %*% bt[1:(m+2)])), plogis(B0 %*% bt[1:(m+2)]))

# u <- runif(nrow(X), plogis(Xl%*%bt), plogis(X%*%bt))

qres <- rnorm(u)
# save(qres, file="processed_data/qr_bdctm_nl.RData")

ggpubr::ggqqplot(data.frame(qres=qres+1), x = "qres", title="BDCTM_nl")




# nonlinear effects on log odds ratio
betas <- fit$beta[1001:2000, 1:ncol(X)]

B1pred <- Predict.matrix(sm1[[1]], data=data.frame(ncountryn=data$ncountry))
B2pred <- Predict.matrix(sm2[[1]], data=data.frame(yearn=data$year))
B3pred <- Predict.matrix(sm3[[1]], data=data.frame(nclaimsn=data$nclaims))

gg_dat <- bind_rows(
ncountry = data.frame(x=data$ncountry, t(apply(betas[,14:24]%*%t(B1), 2, quantile,  probs=c(0.025, 0.5, 0.975)))),
year = data.frame(x=data$year, t(apply(betas[,20:30]%*%t(B2), 2, quantile,  probs=c(0.025, 0.5, 0.975)))),
nclaims = data.frame(x=data$nclaims, t(apply(betas[,31:41]%*%t(B3), 2, quantile,  probs=c(0.025, 0.5, 0.975)))), .id = "effect"
) %>% set_names(c("effect", "x", "lower", "median", "upper"))


p1 <- gg_dat %>% ggplot(aes(x=x, y=median)) + geom_line()+ facet_wrap(~effect, scales="free") + geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.1) + ylab("h") + xlab("")
p1

# ggsave("figures/patent_nl_effects.png", plot=p1, width=12, height=4)
