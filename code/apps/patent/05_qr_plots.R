##
## Script name: patent_qr_plots.R
##
## Purpose of script: produces QR plots of all patent models
##
## Author: Manuel Carlan
##
## Date Created: 2021-10-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------
library(tidyverse)
library(ggpubr)
library(countreg)

data <-read.table("raw_data/patent.raw", header=T)


load("processed_data/qr_bdctm_reg.RData")
qres_lin <- qres
load("processed_data/qr_bdctm_nl.RData")
qres_nl <- qres

load("processed_data/qr_bdctm_hurdle.RData")
qres_hurdle <- qres
load("processed_data/qr_bdctm_hurdle_nl.RData")
qres_hurdle_nl <- qres


qres_bdctm <- bind_rows(BDCTM_lin=qres_lin, BDCTM_nl=qres_nl, BDCTM_hurdle=qres_hurdle, BDCTM_hurdle_nl=qres_hurdle_nl, .id="model")


ggpubr::ggqqplot(data.frame(qres=qres+1), x = "qres", title="BDCTM_hurdle")


# get poisson, negative binomial, zero-inflated poisson and zero-inflated negative binomial data
f1 <-ncit ~ opp+ biopharm  +patus+ patgsgr + ncountryn + yearn + nclaimsn

m.p <-glm(f1, data=data, family="poisson")
m.nb <-glm.nb(f1, data=data)
m.zip <-countreg::zeroinfl(f1, data=data, dist="poisson")
m.zinb <-countreg::zeroinfl(f1, data=data, dist="negbin")

qres <- function(model) {
  data.frame(
    sample = qresiduals(mc, nsim = 3),
    median = qresiduals(mc, type = "quantile"),
    mean100 = rowMeans(qresiduals(mc, nsim = 100)),
    range = qresiduals(mc, type = "quantile", prob = c(0, 1))
  )
}

glm_models <- list(Poisson=m.p, "Negative Binomial"=m.nb, "Zero-inflated Poisson"=m.zip, "Zero-inflated Negative Binomial"=m.zinb)
qres_glm <-glm_models %>% lapply(qresiduals) %>% bind_rows() #%>% pivot_longer(cols=1:4,names_to = "model", values_to = "qr") 

p + stat_qq() + stat_qq_line()

gg_dat <- bind_cols(qres_bdctm, qres_glm)  %>% pivot_longer(cols=everything()) %>% mutate(row = ifelse(name %in% c("Poisson", "Negative Binomial", "BDCTM_lin", "BDCTM_nl"), 1, 2))

p1 <- gg_dat %>% ggplot(aes(sample=value)) + stat_qq() + stat_qq_line() +
  facet_wrap(~ factor(name, levels = c("Poisson", "Negative Binomial", "BDCTM_lin", "BDCTM_nl", "Zero-inflated Poisson", "Zero-inflated Negative Binomial", "BDCTM_hurdle", "BDCTM_hurdle_nl") ), nrow=2)
p1

ggsave("figures/qq_plots.png", plot=p1, height=8)
