
library(here,usethis)
library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(simcausal)
library(survtmle2)
library(SuperLearner)
library(MOSSATE)



get.data <- function(iti,samplesize, conmode, endtime=50,ratDiv){
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

  if (conmode == "scenario 1"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 2"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 3"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4+log(1.2)*W5+log(1.2)*W6+log(2)*W7+log(2)*W8)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "no"){
    D <- D+ node("odds",distr = "rconst", const = 1)+
      node("A", distr = "rbinom", size = 1, prob = .5) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 3.1"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4+log(1.2)*W5+log(1.2)*W6+log(2)*W7+log(2)*W8)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario HTE 1"){
    for (i in 10:30){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='runif', min = 0, max = 1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A+A)")))
    wnames <- grep('W',names(D),value = T)

  }else if(conmode == "scenario HTE 2"){
    for (i in 10:30){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='rbinom', prob = .5,size=1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A*log(1.2)+A)")))

  }else if(conmode == "scenario HTE 3"){
    for (i in 10:60){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='runif', min = 0, max = 1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A*log(1.2)+A)")))
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
surv_simulate <- function(data_out,size,nsim){
  truefit <- MOSSATE::MOSS$new(data_out, dW = 1,
                               verbose = TRUE, epsilon.step = 1e-1, max.iter = 0)
  truefit$onestep_curve(g.SL.Lib = c("SL.gam","SL.glmnet","SL.mean","SL.earth"),
                        Delta.SL.Lib = c("SL.mean","SL.earth","SL.glmnet","SL.gam"),
                        ht.SL.Lib = c("SL.gam","SL.mean","SL.earth","SL.gam_.5","SL.glmnet"),
                        env = parent.frame())
  truefit.1 <- MOSSATE::MOSS$new(data_out, dW = 0, verbose = TRUE, epsilon.step = 1e-1,
                                 max.iter = 0,pred_data = T,
                                 ftimeMod = truefit$ftimeMod,
                                 ctimeMod = truefit$ctimeMod,
                                 trtMod = truefit$trtMod)
  truefit.1$onestep_curve( env = parent.frame())
  truefit.0 <- MOSSATE::MOSS$new(data_out, dW = 1, verbose = TRUE, epsilon.step = 1e-1,
                                 max.iter = 0,pred_data = T,
                                 ftimeMod = truefit$ftimeMod,
                                 ctimeMod = truefit$ctimeMod,
                                 trtMod = truefit$trtMod)
  truefit.0$onestep_curve( env = parent.frame())
  curve0 <- truefit.0$Qn.A1.t_initial
  curve1 <- truefit.1$Qn.A1.t_initial
  T.uniq <- truefit.1$T.uniq
  simulation <- get_data(data_out,size,nsim,T.uniq,curve1,curve0,mode="true_censor")
  return(simulation)
}
get_data <- function(one_sample,size,nsim,T.uniq,curve1,curve0,mode=NULL){
  one_sample$ID <- seq(1,nrow(one_sample),1)

  interval <- c(diff(T.uniq),1)
  #curve_individual_1 <- apply(curve1, 1, function(x) c(1, rep(x,interval)[1:max(T.uniq)])) %>% t()
  #curve_individual_0 <- apply(curve0, 1, function(x) x +pexp(seq(0,max(T.uniq)-1,1)/100, 2)/10%>% t()) %>% t()
  curve_individual_1 <- curve1
  curve_individual_0 <- curve0

  if(mode == "true_censor") curve_c <- curve
  curve_individual_c_1 <- curve1
  curve_individual_c_0 <- curve0


  n <- nrow(curve_individual_0)
  ids <- tnew <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
  data_out <- list()

  for(sim in 1:nsim) {
    idxs <- sample(n, size, replace = TRUE)
    ids[,sim] <- one_sample$ID[idxs]
    # event time
    u <- runif(size, 0, 1)
    # the first time survival drops below u
    stime <-  ifelse(one_sample$A[idxs]==1,
                     apply(curve_individual_1[idxs,] < u, 1, function(x) which(x)[1]),
                     apply(curve_individual_0[idxs,] < u, 1, function(x) which(x)[1]))
    w <- ifelse(one_sample$A[idxs]==1,
                curve_individual_1[idxs,length(T.uniq)] > u,
                curve_individual_0[idxs,length(T.uniq)] > u)
    stime <- T.uniq[stime]
    stime[w] <- max(T.uniq) + 1

    w

    # censoring time
    # u <- runif(size, 0, 1)
    # ctime <-  ifelse(one_sample$A[idxs]==1,
    #                  apply(curve_individual_c_1[idxs,] < u, 1, function(x) which(x)[1]),
    #                  apply(curve_individual_c_0[idxs,] < u, 1, function(x) which(x)[1]))
    # w <- ifelse(one_sample$A[idxs]==1,
    #             curve_individual_c_1[idxs,length(T.uniq)] > u,
    #             curve_individual_c_0[idxs,length(T.uniq)] > u)
    # ctime <- T.uniq[ctime]
    # ctime[w] <- max(T.uniq)
    if(mode == "true_censor") ctime <-ifelse(one_sample$Delta[idxs]==0, one_sample$T.tilde[idxs],max(T.uniq))

    # put it together
    tnew[,sim] <- pmin(stime, ctime)
    names(tnew) <- paste("T.tilde", 1:nsim, sep = "")
    ynew[,sim] <- stime == tnew[,sim]
    names(ynew) <- paste("Delta", 1:nsim, sep = "")
    data_out[[sim]] <- data.frame(one_sample[ids[,sim],!names(one_sample) %in% c("T.tilde","Delta")],
                                  Delta =as.numeric(ynew[,sim]), T.tilde = tnew[,sim])
    data_out[[sim]]$ID <- seq(1:nrow(data_out[[sim]]))
    rownames(data_out[[sim]]) <- NULL

  }

  to.return <- list(data_out = data_out,
                    deltas = ynew,
                    T.tildes = tnew,
                    ids =ids
  )
  return(to.return)
}
SL.gam_.5 <- function (... , deg.gam = .5) {
  SL.gam (... , deg.gam = deg.gam)
}
SL.gam_2 <- function (... , deg.gam = 2) {
  SL.gam (... , deg.gam = deg.gam)
}



#simu data
dat <-  get.data(12345,2000,"scenario 3",endtime =12,ratDiv=200)
data_out <- dat$datt_out
data_out <- data_out[data_out$T.tilde<150 & data_out$T.tilde>0,]
data_out <- data_out[complete.cases(data_out),]
table(data_out$T.tilde)
table(data_out$Delta[data_out$T.tilde<12])/nrow(data_out)
table(data_out$A)/nrow(data_out)
compute_true_effect <- function(x,ratDiv=70){
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+2*(x[,1]+x[,2]+x[,3]+x[,4]))
  s_diff_true_1 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+x[,1]+x[,2]+x[,3]+x[,4])
  s_diff_true_0 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  return(list(effect = s_diff_true_1-s_diff_true_0,
              treated = s_diff_true_1,
              control = s_diff_true_0))
}
True.eff <-compute_true_effect(data_out[,grep("W",names(data_out),value = T)],ratDiv=200)



endtime <- 12
tpts <- seq(1:endtime)

cores <- parallel::detectCores()
options(mc.cores = cores)
set.seed(25431, "L'Ecuyer-CMRG")
try(detach('package:mgcv',unload=TRUE),silent = T)
onestepfit <- MOSSATE::MOSS$new(data_out, dW = 1, verbose = T, epsilon.step = 1e-2, max.iter = 50)
onestepfit$initial_fit()
onestepfit$onestep_curve(g.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         Delta.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         ht.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                         env = parent.frame())
onestepfit$Psi.hat[1:endtime]
colMeans(onestepfit$Qn.A1.t_initial)[1:endtime]

MOSS_fit <- MOSS::MOSS$new(dat = data_out, dW = 1, epsilon.step = 1e-3, max.iter = 0, verbose = T)
MOSS_fit$onestep_curve(g.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                       Delta.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"),
                       ht.SL.Lib = c("SL.mean", "SL.glm", "SL.earth"))

MOSS_fit$Psi.hat[1:endtime]


treated.eff.true.avg <- colMeans(True.eff$treated)
treated.eff.true.avg
control.eff.true.avg <- colMeans(True.eff$control)






onestepdiff <- MOSSATE::MOSS$new(data_out, verbose = TRUE, epsilon.step = 1e-2,
               max.iter = 10,
               ftimeMod = onestepfit$ftimeMod,
               ctimeMod = onestepfit$ctimeMod,
               trtMod = onestepfit$trtMod)

onestepdiff <- MOSSATE::MOSS_difference$new(data_out, verbose = TRUE, epsilon.step = 1e-2,
                                            max.iter = 10,
                                            ftimeMod = onestepfit$ftimeMod,
                                            ctimeMod = onestepfit$ctimeMod,
                                            trtMod = onestepfit$trtMod)

onestepfit$initial_fit()

onestepdiff$initial_fit(env = parent.frame())
onestepdiff$onestep_diff_curve()

sd_EIC <- (onestepdiff$D1.t)[,1:endtime]


A <- onestepdiff$dat$A
controls <- onestepdiff$dat[,grep("W",names(onestepdiff$dat),value = T)]

#Get the true effect by treatment
testdat <- data.frame(A=A,controls)

#Plug in the estimator
Y.hat.1 <- onestepdiff$MOSS_A1$Qn.A1.t_initial[,1:endtime]
Y.hat.0 <- onestepdiff$MOSS_A0$Qn.A1.t_initial[,1:endtime]
prediction <- Y.hat.1-Y.hat.0
prediction_tmle <- onestepdiff$MOSS_A1$Qn.A1.t_full[,1:endtime]-
  onestepdiff$MOSS_A0$Qn.A1.t_full[,1:endtime]






