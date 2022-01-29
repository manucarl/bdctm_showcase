##
## Script name: forest.R
##
## Purpose of script: estimates discrete BDCTM in form of a partial proportional odds model for forest health data with spatial tensor spline and random effectn
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

rm(list=ls())

packages <- c("MASS","BayesX", "Rcpp", "RcppArmadillo", "splines", "mgcv", "Matrix",
              "R2BayesX", "tidyverse", "profvis", "numDeriv", "scales", 
              "doParallel", "scam", "mvtnorm", "MCMCpack", "rlist", "RhpcBLASctl")# library(doParallel)

setwd("D:/OneDrive/bdctm_paper/")
lapply(packages, library, character.only = TRUE) %>%  invisible()


source("code/nuts/adnuts_helper.R")
source("code/nuts/nuts_forest_te_omega_ig.R")
sourceCpp("code/rcpp/log_disc_xx2.cpp")





## load forest health data and tree location map
data <-read.table("raw_data/beach.raw", header=T)

#################################
# prepare data
#################################

data <-
  data %>% as_tibble  %>% mutate(
    defol3 = case_when(
      defol %in% 0 ~ 0,
      defol %in% c(12.5, 25, 37.5) ~ 1,
      defol >= 50 ~ 2
    )) %>% as.data.frame

data$x_s <- scale(data$x)
data$y_s <- scale(data$y)
data$canopyd_s <- scale(data$canopyd)

n <- nrow(data)

# nonlinear covariate smooths
smooths <- NULL
Slist <-NULL


#################################
# construct trafos
#################################


ystore <- y <- data$defol3
ymax <- max(ystore)

# no of linear effects
p0 <- 0
yfun <- function(x) x


ystore <- y <- data$defol3
ymax <- max(ystore)

# build discrete basis and its lagged basis
y <- as.numeric(y)
B <-  matrix(0, nrow=n, ncol=length(unique(y)))
Bl <- matrix(0, nrow=n, ncol=length(unique(y)))

for(i in 1:n){
  B[i,y[i]+1 ] <- 1
  Bl[i,y[i]] <-1
}


B <- B[,-3]
Bl <- Bl[, -3]
# stair matrix
Sigma <- lower.tri(diag(ncol(B)), diag = T)

B0 <- B%*%Sigma
B0l <-  Bl%*%Sigma


# precision matrix
K0 <- matrix(0, ncol(B0), ncol(B))

X <- B0
Xl <- B0l


# # # construct non-proportional effect
sspec <- s(canopyd_s, k=6+2,bs="ps", m=1)
sm <- smoothCon(sspec,data=data,
                scale.penalty=F, absorb.cons=T, knots=NULL)[[1]]  

SigmaStar <- Sigma %x% diag(ncol(sm$X))
B01 <- tensor.prod.model.matrix(list(B, sm$X)) %*% SigmaStar
B01l <- tensor.prod.model.matrix(list(Bl, sm$X)) %*% SigmaStar

colnames(B01) <- paste0("vcm(canopyd_s).", 1:ncol(B01))

# penalty matrix
K01_2 <- diag(2)%x%sm$S[[1]]


# # random effect
Br <- model.matrix(~as.factor(data$id)-1)[,-1] # %>% scale
colnames(Br) <- paste0("re_",sort(unique(data$id))[-1])
Kr <- diag(ncol(Br))




# spatial tensor spline

# basis dimension of marginal effects
m_no <- 10

sspec <- te(x_s, y_s, k=c(m_no, m_no), m=list(c(2,1),c(2,1)),  bs="ps")
smspat <- smoothCon(sspec,data=data, sparse.cons= -1,
                    scale.penalty=F, absorb.cons=T, knots=NULL)[[1]]
#

margins <- smooth.construct(sspec,data=data,  knots=NULL)[[1]]
Bxy1 <- margins[[1]]$X
Bxy2 <- margins[[2]]$X
drop <- m_no

Bxy <- tensor.prod.model.matrix(list(Bxy1, Bxy2))

C <- matrix(colMeans(Bxy), 1, ncol(Bxy))


qrc <- c(drop, as.numeric(C)[-drop])
Bxy <- Bxy[, -drop, drop = FALSE] - matrix(qrc[-1], 
                                           nrow(Bxy), ncol(Bxy) - 1, byrow = TRUE)

Kxy_1 <- (margins[[1]]$S[[1]] %x% diag(ncol(Bxy1)))[-drop, -drop]
Kxy_2 <- (diag(ncol(Bxy1)) %x% margins[[2]]$S[[1]])[-drop, -drop]


colnames(Bxy) <- paste0("te(x,y).", 1:ncol(Bxy))


#################################
# model transformation matrices
#################################

X <- cbind(B0, B01, -Br, - Bxy) %>% as.matrix
Xl <- cbind(B0l,B01l, - Br, -Bxy) %>% as.matrix

# # List of Precision matrices
Klist <- list(K0, K01_2, Kr, list(Kxy_1, Kxy_2))

# diagonal precision matrix
S <- bdiag(K0, K01_2, Kr, Kxy_1 + Kxy_2)  %>% as.matrix

ind_re <- which(startsWith(colnames(X), "re"))
ind_npo <- which(startsWith(colnames(X), "vcm("))

ind_spat <- which(startsWith(colnames(X), "te("))

# indices of exponentiated coefs
ind_exp <- c(2, ind_npo[-1:-7])
p <- ncol(X)

beta <- runif(p)
bt <- beta
bt[ind_exp] <- exp(beta[ind_exp])


zero_ind <- which(ystore==0)
max_ind <- which(y==max(y))

Xnew <- X
Xlnew <- Xl

Xzero <- X[ystore ==0,  ]
Xrest <- X[ystore>0,  ]
Xlrest <- as.matrix(Xl[ystore>0,  ])
Xzerof <- Xzero[,-c(3) ]
Xrestf <- Xrest[,-c(3)]
Xlrestf <- Xlrest[,-c(3)]

Xlone <- Xl[ystore==2,]
Xone <- matrix(1,  nrow(Xlone), 1)

hyperparams=list(a=0.001, b=0.001)

library(MASS)
xx <- list(y=ystore, X=X, Xl=Xl,
           S= S,
           Xzero=Xzero, Xrest=Xrest, Xlrest=Xlrest, Xzerof=Xzerof, Xrestf=Xrestf, Xlrestf=Xlrestf,
           Xnew = Xnew, Xlnew=Xlnew, zero_ind=zero_ind-1, max_ind = max_ind-1,
           Klist=Klist,  ymax=ymax, 
           ind_exp=ind_exp-1, p=p, p0=p0, pnlin=length(Klist), hyperparams=hyperparams)

settings=list(max_treedepth=12, adapt_delta=0.95, verbose=TRUE)

its=10000
warmup=its/2
burnin=its/2

# starting values
start <- rep(0, p)
fit <-NUTS(n_iter=its, xx=xx, f=posterior_logitd3, gr=gradf_logitd3, ll=ll_logitd3, fixed=fixed, start=start, warmup=warmup, nuts_settings=settings)

# first order differences:
# Final acceptance ratio=0.88, and target=0.9
# Final step size=0.033; after 5000 warmup iterations
# Elapsed Time: 9498.3 seconds (Warmup)
# Elapsed Time: 8489.1 seconds (Sampling)
# Elapsed Time: 17987.4 seconds (Total)

#second order differences
# Final acceptance ratio=0.71, and target=0.95
# Final step size=0.01; after 5000 warmup iterations
# Elapsed Time: 26770.0 seconds (Warmup)
# Elapsed Time: 19220.0 seconds (Sampling)
# Elapsed Time: 45990.0 seconds (Total)


# save(fit, file="processed_data/forest_npo_re_spat.RData")

# load("processed_data/forest_npo_re_spat.RData")
load("processed_data/forest_npo_id_re_spat_omega_ig_m2_10000it.RData")

library(mcmcplots)

beta_samples <- fit$beta[1000:10000,]
colnames(beta_samples) <- colnames(xx$Xnew)
# mcmcplot(beta_samples)
tau2_samples <- fit$tau[,1:3]
colnames(tau2_samples) <- c("tau2_npo_effect", "tau2_random_effect", "tau2_spatial_effect")
# mcmcplot((tau2_samples))
# effectiveSize(beta_samples)
betas <- beta_samples#[500:1000,]
bts <- betas
bts[,ind_exp] <- exp(bts[,ind_exp])

beta <- colMeans(betas)
colnames(betas) <- colnames(xx$Xnew)
bt <- beta
bt[ind_exp] <- exp(bt[ind_exp])



bt_npo <- bt[ind_npo]


#############################
# spatial effect
############################

# granularity of grid
n_grid <- 200
pred_grid <- expand.grid(x=seq(min(data$x), max(data$x), length=n_grid), y=seq(min(data$y), max(data$y), length=n_grid))
pred_grid <- expand.grid(x_s=seq(min(data$x_s), max(data$x_s), length=n_grid), y_s=seq(min(data$y_s), max(data$y_s), length=n_grid))


beta_spat <- betas %>% as_tibble %>% dplyr::select(starts_with("te(x,y)")) %>% colMeans
pred_spat <- PredictMat(smspat, pred_grid)%*% beta_spat


gg_spat <- data.frame(pred_spat, pred_grid)

spat_plot <- ggplot( data=gg_spat) + geom_tile(aes(x=x_s, y=y_s, fill=pred_spat))  + 
  geom_point(aes(x_s, y_s), data=data, shape=2, size=5) +
  scale_fill_viridis_c(option="A", name="h(y,x)", direction=-1, begin=0.1, end=1) +
  theme_void() +
  expand_limits(x=0, y=0) 

spat_plot


# ggsave("figures/forest_spat_m2.png", plot=spat_plot, width=10)


#############################
# random effect
############################

beta_rand <- betas  %>% as_tibble %>% dplyr::select(starts_with("re"))

# > paste0(colMeans(beta_rand)) %>% as.numeric %>% var
# [1] 2.416654
# > paste0(colMeans(beta_rand)) %>% as.numeric %>% sd
# [1] 1.554559

gg_rand <- beta_rand %>% pivot_longer(everything())


library(ggmcmc)
colnames(gg_rand) <- c("Parameter", "value") 
gg_rand$Parameter <- as.factor(gg_rand$Parameter)

re_plot <- ggs_caterpillar(gg_rand, horizontal=F) + ylab("id") +xlab("")
re_plot
# ggsave("figures/forest_random.png", plot=re_plot, height=4)


#############################
# nonlinear npo effect plot
############################
sspec <- s(canopyd_s, k=6+2,bs="ps", m=1)
sm <- smoothCon(sspec,data=data,
                scale.penalty=F, absorb.cons=T, knots=NULL)[[1]] 
Sigma <- lower.tri(diag(2), diag = T)
SigmaStar <- Sigma %x% diag(ncol(sm$X))

B01pred <- matrix(c(1,0),nrow=n, ncol=2, byrow=T)%.%sm$X%*%SigmaStar
B02pred <- matrix(c(0,1),nrow=n, ncol=2, byrow=T)%.%sm$X%*%SigmaStar

effquants <- function(x, ind, betatilde_samples){
  effquan <-apply(x%*%t(betatilde_samples[,ind]), 1, quantile, probs= c(0.025, 0.975))
  
  res <-tibble(q025=effquan[1,], med =(x%*%apply(betatilde_samples[,ind], 2, median)),  q975= effquan[2,])
  return(res)
}

gg_data <- bind_rows(canopy= bind_rows(tibble(x=data$canopyd_s, effquants(B01pred, ind_npo, bts), cat=1),
                                       tibble(x=data$canopyd_s,effquants(B02pred, ind_npo,bts), cat=2)),
                     .id="effect") 

npo_plot <- gg_data %>% 
  ggplot(aes(x=x, y=med[,1])) + geom_line() + facet_wrap(cat~effect, scale="free", ncol=3) +
geom_ribbon(aes(ymin=q025,ymax=q975),alpha=0.3)
npo_plot

# ggsave("figures/forest_npo.png", plot=npo_plot, height=4)
