library(here,usethis)
library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(simcausal)
library(survtmle2)
library(SuperLearner)
library(MOSSATE)
library(MOSS)




get.data <- function(iti,samplesize, conmode, endtime,ratDiv){
  D <- DAG.empty()
  D <- D +
    node("W1", distr ="rbinom", prob = .5,size=1)+
    node("W2", distr ="runif", min = 0, max = 1)+
    node("W3", distr ="rbinom", prob = .5,size=1)+
    node("W4", distr ="runif", min = 0, max = 1)+
    node("W5", distr ="rbinom", prob = .5,size=1)+
    node("W6", distr ="runif", min = 0, max = 1)+
    node("W7", distr ="rbinom", prob = .5,size=1)+
    node("W8", distr ="runif", min = 0, max = 1)+
    node("W9", distr ="runif", min = 0, max = 1)+
    node("W10", distr ="runif", min = 0, max = 1)

  if(conmode == "scenario 3"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4+log(1.2)*W5+log(1.2)*W6+log(2)*W7+log(2)*W8)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "no"){
    D <- D+ node("odds",distr = "rconst", const = 1)+
      node("A", distr = "rbinom", size = 1, prob = .5) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }

  D <- D+
    node("rate1", distr = "rconst", const = rate/ratDiv) +
    node("Trexp", distr = "rexp", rate = rate1) +
    node("Cweib", distr = "rweibull", shape = .8 - .1*W1, scale = 20) +
    node("T", distr = "rconst", const = round(Trexp,0)) +
    node("C", distr = "rconst", const = round(Cweib,0)) +
    #node("C", distr = "rconst", const = T+1) +
    node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
    node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
  setD <- set.DAG(D)

  dat <- sim(setD,n=samplesize,rndseed= iti)

  data_out <- dat[,names(dat) %in% c("ID",wnames,"A","T.tilde","Delta" )]

  return(list(datt_out = as.data.frame(data_out)))
}

#simulate data
dat <-  get.data(12345,2000,"scenario 3",endtime =12,ratDiv=450)
data_out <- dat$datt_out
data_out <- data_out[data_out$T.tilde<150 & data_out$T.tilde>0,]
data_out <- data_out[complete.cases(data_out),]
table(data_out$T.tilde)
#Show the event rate before time 12
table(data_out$Delta[data_out$T.tilde<12])/nrow(data_out)
compute_true_effect <- function(x,ratDiv){
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+2*(x[,1]+x[,2]+x[,3]+x[,4]))
  s_diff_true_1 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+x[,1]+x[,2]+x[,3]+x[,4])
  s_diff_true_0 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  return(list(effect = s_diff_true_1-s_diff_true_0,
              treated = s_diff_true_1,
              control = s_diff_true_0))
}
True.eff <-compute_true_effect(data_out[,grep("W",names(data_out),value = T)],ratDiv=450)
treated.eff.true.avg <- colMeans(True.eff$treated)
control.eff.true.avg <- colMeans(True.eff$control)


endtime <- 12
cores <- parallel::detectCores()
options(mc.cores = cores)
set.seed(25431, "L'Ecuyer-CMRG")
try(detach('package:mgcv',unload=TRUE),silent = T)

#Test estimation of one survival curve
ATEfit <- MOSSATE::MOSS$new(data_out, dW = 1, verbose = T, epsilon.step = 1e-3, max.iter = 50)
ATEfit$initial_fit()
ATEfit$onestep_curve(g.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         Delta.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         ht.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         env = parent.frame())
ATEfit$Psi.hat[1:endtime]
colMeans(ATEfit$Qn.A1.t_initial)[1:endtime]

#Experiments at dw=1
paste(round((treated.eff.true.avg-colMeans(ATEfit$Qn.A1.t_initial)[1:endtime])/treated.eff.true.avg,6)*100,"%")
#Bias before adjustment
# "1.2284 %" "1.126 %"  "1.0288 %" "0.9368 %" "0.85 %"   "0.7686 %" "0.6924 %" "0.6215 %" "0.5561 %" "0.496 %"  "0.4413 %" "0.3921 %"
paste(round((treated.eff.true.avg-ATEfit$Psi.hat[1:endtime])/treated.eff.true.avg,6)*100,"%")
#Bias after adjustment
#"0.2657 %" "0.2687 %" "0.2789 %" "0.2963 %" "0.3211 %" "0.3532 %" "0.3927 %" "0.4396 %" "0.4941 %" "0.5562 %" "0.6258 %" "0.7032 %"


#Experiments at dw=0
paste(round((control.eff.true.avg-colMeans(ATEfit$Qn.A1.t_initial)[1:endtime])/treated.eff.true.avg,6)*100,"%")
#Bias before adjustment
#"1.0953 %" "1.3083 %" "1.5229 %" "1.739 %"  "1.9567 %" "2.176 %"  "2.3968 %" "2.6192 %" "2.8432 %" "3.0688 %" "3.296 %"  "3.5247 %"
paste(round((control.eff.true.avg-ATEfit$Psi.hat[1:endtime])/treated.eff.true.avg,6)*100,"%")
#Bias after adjustment
#"1.1368 %" "1.5428 %" "1.9534 %" "2.3686 %" "2.7885 %" "3.213 %"  "3.6423 %" "4.0763 %" "4.515 %"  "4.9587 %" "5.4072 %" "5.8605 %"



MOSS_fit <- MOSS::MOSS$new(dat = data_out, dW = 1, epsilon.step = 1e-3, max.iter = 50, verbose = T)
MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                       Delta.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                       ht.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"))

paste(round((treated.eff.true.avg-MOSS_fit$Psi.hat[1:endtime])/treated.eff.true.avg,6)*100,"%")

#Experiments at dw=1
#Bias before adjustment (by setting epsilon.step = 1e-10, max.iter = 1)
# "0.1814 %" "0.1888 %" "0.2024 %" "0.2224 %" "0.2486 %" "0.2812 %" "0.3203 %" "0.3657 %" "0.4177 %" "0.4762 %" "0.5413 %" "0.613 %"
#Bias after adjustment
#"0.169 %"  "0.1812 %" "0.1996 %" "0.2242 %" "0.2552 %" "0.2924 %" "0.3361 %" "0.3862 %" "0.4427 %" "0.5058 %" "0.5754 %" "0.6515 %"

#Experiments at dw=0
paste(round((control.eff.true.avg-MOSS_fit$Psi.hat[1:endtime])/treated.eff.true.avg,6)*100,"%")
#Bias before adjustment
# "0.2022 %" "0.5621 %" "0.929 %"  "1.303 %"  "1.6841 %" "2.0724 %" "2.468 %"  "2.8708 %" "3.281 %"  "3.6986 %" "4.1237 %" "4.5564 %"
#Bias after adjustment
# "0.1695 %" "0.5597 %" "0.9567 %" "1.3608 %" "1.7718 %" "2.19 %"   "2.6153 %" "3.0479 %" "3.4876 %" "3.9347 %" "4.3892 %" "4.851 %



