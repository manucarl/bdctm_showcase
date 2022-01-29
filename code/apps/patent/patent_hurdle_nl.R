##
## Script name: patent_hurdle_nl.R
##
## Purpose of script: estimates BDCTM_hurdle_nl and  quantile residuals
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

source("code/nuts/nuts_patent_hurdle_nl.R")


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
log_plus <- TRUE
m <- 6

# construct trafo
ind_exp=2:(m+2)
ystore <- y <- data$ncit


# construct monotonic transformation
sspec <- s(y,k=m+2,bs="ps")
sspec$mono <- 1
yfun <- function(x) x
if(log_plus) yfun <- function(x) log(x+1)

sm0 <- smoothCon(sspec,data=data.frame(y=yfun(ystore)), scale.penalty=F, absorb.cons=F ,knots=NULL)[[1]]
B0 <- sm0$X
K0 <- sm0$S[[1]]
B0l <-Predict.matrix(sm0, data=data.frame(y=yfun((ystore-1))))

colnames(B0) <- colnames(B0l) <- paste0("h(",all.vars(formula)[1L], ").", 1:(m+2 ))


# nonlinear shift effects
sm1 <- smoothCon(s(ncountryn, k=10+2,bs="ps"), absorb.cons=TRUE, data=data)
B1 <- sm1[[1]]$X
K1 <- sm1[[1]]$S[[1]]

sm2 <- smoothCon(s(yearn,k=10+2,bs="ps"), absorb.cons=TRUE, data=data)
B2 <- sm2[[1]]$X
K2 <- sm2[[1]]$S[[1]]

sm3 <- smoothCon(s(nclaimsn,k=10+2,bs="ps"), absorb.cons=TRUE, data=data)
B3 <- sm3[[1]]$X
K3 <- sm3[[1]]$S[[1]]


x <- cbind(B1, B2, B3)
S <- as.matrix(bdiag(K0, K1, K2, K3))

# precision matrix
Klist <- list(K0,  K1, K2, K3)

trafo <- "hurdle"

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
plin <- 1

Xzero <- X[ystore ==0,  ]
Xrest <- X[ystore >0,  ]
Xlrest <- Xl[ystore >0,  ]


start <- rep(0,p)
hyperparams=list(a=1, b=0.001)



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
xx <- list(y=ystore, X=X, Xl=Xl, Xzero=Xzero, Xrest=Xrest, Xlrest=Xlrest, 
           S= S,
           Klist=Klist,  
           ind_exp=ind_exp-1, p=p, pnlin=length(Klist), p0=1, hyperparams=hyperparams)


settings <- list(max_treedepth=12, adapt_delta=0.8)
its <- 2000
warmup <- burnin <- its/2

fit <-NUTS(n_iter=its, xx=xx, f=posterior, gr=gradf, ll=ll, start=start, warmup=warmup, nuts_settings=settings)
# save(fit, file="processed_data/patent_hurdle_nl.RData")

# Final acceptance ratio=0.95, and target=0.95
# Final step size=0.037; after 1000 warmup iterations
# Elapsed Time: 2059.7 seconds (Warmup)
# Elapsed Time: 1571.0 seconds (Sampling)
# Elapsed Time: 3630.7 seconds (Total)


load("processed_data/patent_hurdle_nl.RData")

beta_samples <- fit$beta
effectiveSize(beta_samples %>% as.mcmc)
# var1      var2      var3      var4      var5      var6      var7      var8      var9     var10     var11     var12     var13     var14     var15     var16     var17 
# 776.7715  607.7740  927.6414 1402.7601 1450.6606 1804.1709 1307.0275 1464.6537  868.0617  504.4394  532.0313  687.5327  861.0906  690.0902  416.4634  530.8789  406.1412 
# var18     var19     var20     var21     var22     var23     var24     var25     var26     var27     var28     var29     var30     var31     var32     var33     var34 
# 642.9123  898.3313  756.7210 1578.5280  791.8142  792.1035  946.2911  653.7324  523.4653  833.4510  830.9004  563.2838 1223.1192 1371.4679  799.6989  780.4607 1076.6510 
# var35     var36     var37     var38     var39     var40     var41     var42     var43     var44     var45     var46     var47     var48     var49     var50     var51 
# 1104.7309 1038.0987 1192.1717 1054.3736  989.9367  925.8560  781.3127  804.6572  971.2084  850.3077  826.7423  854.9895 1058.8050  888.5885  852.7303  706.0320  691.1924 
# var52     var53     var54     var55     var56     var57     var58     var59     var60     var61     var62     var63     var64     var65     var66     var67     var68 
# 707.4885  870.4873  714.9014 1266.5771 1161.9330  887.2175  853.1411  785.3888  740.9553  741.0475  641.3655  962.1623 1752.4439 1599.3745 1137.0214 1185.5976 1261.2040 
# var69     var70     var71     var72     var73     var74     var75 
# 1316.9055 1217.0827 1316.1331 1253.6892  936.1005  829.0624  902.3165

beta <- fit$beta[1001:2000, 1:ncol(X)] %>% colMeans
bt <- beta
bt[ind_exp] <- exp(bt[ind_exp])


X <- as.matrix(X)
Xl <- as.matrix(Xl)

# Quantile Residuals

h0_ind <- 42
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


# save(qres, file="processed_data/qr_bdctm_hurdle_nl.RData")
ggpubr::ggqqplot(data.frame(qres=qres+1), x = "qres", title="BDCTM_hurdle_nl")

