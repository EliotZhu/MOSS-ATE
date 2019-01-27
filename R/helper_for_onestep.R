#' row mulitiplication
#' @export
multiply_vector_to_matrix <- function(M, V) {
  t(t(M) * V)
}

#' generate step function object using y value and corresponding jump points
#'
#' @param y.vec a vector of step function values
#' @param t.vec a vector for all jump locations of step function.
#' NOTE: the first element value of y.vec is the flat part LEFT of the first jump point specified in t.vec
#'
#' @return sfun a function object, that can perform any mapping
#' @export
create_step_func <- function(y.vec, t.vec) {
  if (length(y.vec) != (length(t.vec) + 1)) {
    warning("the legnth of input vectors incorrect!")
  }
  # before the first jump point, the step function is 0 value
  sfun <- stepfun(t.vec, y.vec, f = 0)
  return(sfun)
}

#' compute \eqn{I\{T.tilde >= t\}}
#'
#' loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
create_Yt_vector <- function(Time, t.vec) {
  (Time >= t.vec) + 0
}

#' compute \eqn{I\{T.tilde == t, Delta = 1\}}
#'
#' loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param Delta length n vector of censoring indicator
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
create_Yt_vector_with_censor <- function(Time, Delta, t.vec) {
  ((Time == t.vec) & (Delta == 1)) + 0
}


#' compute cumulative distribution function of a step-shaped (empirical) density
#'
#' @param pdf.mat if input vector = compute cdf for a step-function pdf;
#'              if input matrix = compute cdf for several step-function pdf with same jump points
#' @param t.vec unique jump points of step function
#' @param start if -Inf = from left to right; if Inf = from right to left.
#'
#' @return vector of cdf value
#' @export
compute_step_cdf <- function(pdf.mat, t.vec, start = -Inf) {
  interval.size <- diff(t.vec)
  # interval.size <- c(0, interval.size)
  interval.size <- c(interval.size, 0)

  # compute the mass
  if (is.matrix(pdf.mat)) { # if input with multi-sample
    mass.by.interval <- sweep(pdf.mat, MARGIN = 2, interval.size, `*`)
    # multiplies the interval length to each row of the y-values
    # the result is a matrix, each row is a single pdf, and entries are the mass
  } else { # if input with one-sample
    mass.by.interval <- pdf.mat * interval.size
  }

  if (is.infinite(start) & (start < 0)) { # start from -Inf
    if (is.matrix(pdf.mat)) { # if input with multi-sample
      cdf.by.interval <- t(apply(mass.by.interval, 1, cumsum)) # cumsum of mass for each row, from left to right
    } else { # if input with one-sample
      cdf.by.interval <- cumsum(mass.by.interval)
    }
  } else { # start from +Inf
    if (is.matrix(pdf.mat)) { # if input with multi-sample
      cdf.by.interval <- t(apply(mass.by.interval, 1, function(obj) rev(cumsum(rev(obj)))))
    } else { # if input with one-sample
      cdf.by.interval <- rev(cumsum(rev(mass.by.interval)))
    }
  }
  return(cdf.by.interval)
}


#' compute l2 inner product of two step functions
#'
#' f and g
#'
#' @param f.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param g.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param T.grid shared jump points
#'
#' @return scalar
#' @export
l2_inner_prod_step <- function(f.step, g.step, T.grid) {
  if (is.vector(f.step) & is.vector(g.step)) {
    # both f and g are one sample
    f.times.g <- f.step * g.step
  }
  if (!is.vector(f.step) & is.vector(g.step)) {
    # f: multi-sample
    # g: one-sample
    f.times.g <- sweep(f.step, MARGIN = 2, g.step, `*`) # multiply g to each row of f.
  }
  if (is.vector(f.step) & !is.vector(g.step)) {
    # f: one-sample
    # g: multi-sample
    f.times.g <- sweep(g.step, MARGIN = 2, f.step, `*`) # multiply f to each row of g.
  }
  if (!is.vector(f.step) & !is.vector(g.step)) {
    # both f and g are multi-sample of same sample size
    if (nrow(f.step) != nrow(g.step)) stop("f and g have different sample size!")
    f.times.g <- f.step * g.step
  }
  result <- compute_step_cdf(f.times.g, T.grid)
  if (!is.vector(f.step) | !is.vector(g.step)) {
    # there is multi-sample
    result <- apply(result, 1, function(obj) tail(obj, 1))
  } else {
    # both f and g are one-sample
    result <- tail(result, 1)
  }
  return(result)
}


#' Perform one-step TMLE update of survival curve
#'
#' @param D1.t.func.prev n*p matrix of previous influence curve
#' @param Pn.D1.func.prev p vector of previous mean influence curve
#' @param dat input data.frame
#' @param T.uniq grid of unique event times
#' @param W_names vector of the names of baseline covariates
#' @param dW dynamic intervention
#'
#' @return vector of length n_sample
#' @export
#' @importFrom dplyr left_join
compute_onestep_update_matrix <- function(
  D1.t.func.prev, Pn.D1.func.prev, dat, T.uniq, W_names, dW) {
  # calculate the number inside exp{} expression in submodel
  # each strata of Q is updated the same
  # numerator <- l2_inner_prod_step(-abs(Pn.D1.func.prev), D1.t.func.prev, T.uniq)
  # numerator <- l2_inner_prod_step(abs(Pn.D1.func.prev), D1.t.func.prev, T.uniq) # wrong?
  # numerator <- sweep(D1.t.func.prev, MARGIN=2, -abs(Pn.D1.func.prev),`*`) # mat useful
  numerator <- sweep(D1.t.func.prev, MARGIN = 2, abs(Pn.D1.func.prev), `*`) # colsum useful
  # numerator <- sweep(D1.t.func.prev, MARGIN=2, abs(Pn.D1.func.prev),`*`) # WROOOOOONG
  result <- numerator /
    sqrt(l2_inner_prod_step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

  return(result)
}


#' Estimation for the Method of Cause-Specific Hazards
#'

estimateCensoring <- function(dataList,
                              adjustVars,
                              t0,
                              SL.ctime = NULL,
                              returnModels = FALSE,
                              verbose = TRUE,
                              gtol = 1e-3,
                              env,
                              testdat,
                              ...) {
  include <- !(dataList[[1]]$t == dataList[[1]]$ftime &
                 dataList[[1]]$C != 1 &
                 dataList[[1]]$t < t0)&
    !(dataList[[1]]$t == dataList[[1]]$ftime &
        dataList[[1]]$C == 1 &
        dataList[[1]]$t == t0)


  if (!all(dataList[[1]]$C == 0)) {
    fit.dat.1 <- as.data.frame(dataList[[1]][include & dataList[[1]]$trt==1,])
    message(paste0("Get censor model 1..."))
    ctimeMod.1 <- SuperLearner::mcSuperLearner(
      Y = fit.dat.1$C,
      X = fit.dat.1[, c("t",names(adjustVars)),drop=F],
      id = fit.dat.1$id,
      family = "binomial",
      SL.library = SL.ctime,
      verbose = verbose,
      env = env
    )

    fit.dat.0 <- as.data.frame(dataList[[1]][include & dataList[[1]]$trt==0,])
    message(paste0("Get censor model 0..."))
    ctimeMod.0 <- SuperLearner::mcSuperLearner(
      Y = fit.dat.0$C,
      X = fit.dat.0[, c("t",names(adjustVars)),drop=F],
      id = fit.dat.0$id,
      family = "binomial",
      SL.library = SL.ctime,
      verbose = verbose,
      env = env
    )

  } else {
    dataList <- lapply(dataList, function(x) {
      x$G_dC <- 1
      x
    })
    ctimeMod.1 <- "No censoring observed"
    class(ctimeMod.1) <- "noCens"
  }

  if (class(ctimeMod.1) != "noCens") {
    get.pred.censor <- function(x,ctimeMod) {
      g_dC <- rep(1, nrow(x))
      if (t0 != 1) {
        # temporarily replace time with t-1
        # NOTE: this will fail if t enters model as a factor
        x$t <- x$t - 1
        g_dC <-
          suppressWarnings(
            1 - predict(ctimeMod,newdata = x[, c("t",names(adjustVars))], onlySL = TRUE)[[1]]
          )
        # put time back to normal
        x$t <- x$t + 1
        # replace any observations with t = 1
        # to avoid extrapolation at t = 0
        g_dC[x$t == 1] <- 1
      }
      x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
      return(x)
    }
    message("Finalizing censoring model...")

    dataList[[2]] <- get.pred.censor(dataList[[2]] , ctimeMod.0)
    dataList[[3]] <- get.pred.censor(dataList[[3]] , ctimeMod.1)

  } else {
    dataList <- lapply(dataList, function(x) {
      x$G_dC <- 1
      x
    })
  }

  # truncate small propensities at gtol

  dataList$`0`$G_dC[dataList$`0`$G_dC < gtol] <- gtol
  dataList$`1`$G_dC[dataList$`1`$G_dC < gtol] <- gtol
  if (class(ctimeMod.1) != "noCens") {
    out <- list(
      dataList = dataList,
      ctimeMod = if (returnModels) { list(ctimeMod.1,ctimeMod.0,get.pred.censor) } else {NULL}
    )
  }else{
    out <- list(
      dataList = dataList
    )
  }

  return(out)
}

estimateHazards <- function(dataList,
                            J,
                            adjustVars,
                            SL.ftime = NULL,
                            returnModels,
                            bounds,
                            verbose,
                            env,
                            ...) {
  ftimeMod.1 <- vector(mode = "list", length = length(J))
  names(ftimeMod.1) <- paste0("J", J)
  ftimeMod.0 <- vector(mode = "list", length = length(J))
  names(ftimeMod.0) <- paste0("J", J)

  for (j in J) {
    # add all events less than current j to see who to include in regression
    NlessthanJ <- rep(0, nrow(dataList[[1]]))
    for (i in J[J < j]) {
      NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
    }

    fit.dat.1 <- as.matrix(dataList[[1]][NlessthanJ == 0 & dataList[[1]]$trt==1,]) %>% as.tibble(.)
    message(paste0("Get model 1.",j))
    Qj_mod.1 <- SuperLearner::mcSuperLearner(
      Y = fit.dat.1[[paste0("N", j)]],
      X = fit.dat.1[,c("t", names(adjustVars))],
      id = fit.dat.1$id,
      family = stats::binomial(),
      SL.library = SL.ftime,
      env =env
    )

    fit.dat.0 <- as.matrix(dataList[[1]][NlessthanJ == 0 & dataList[[1]]$trt==0,]) %>% as.tibble(.)
    message(paste0("Get model 0.",j))
    Qj_mod.0 <- SuperLearner::mcSuperLearner(
      Y = fit.dat.0[[paste0("N", j)]],
      X = fit.dat.0[,c("t", names(adjustVars))],
      id = fit.dat.0$id,
      family = stats::binomial(),
      SL.library = SL.ftime,
      env=env
    )

    ftimeMod.1[[paste0("J", j)]] <- Qj_mod.1
    ftimeMod.0[[paste0("J", j)]] <- Qj_mod.0

    # get predictions back
    get.pred.hazard <- function(x, j, Qj_mod) {
      suppressWarnings(
        x[[paste0("Q", j, "PseudoHaz")]] <- predict(
          Qj_mod,
          onlySL = TRUE,
          newdata = x[, c("t", names(adjustVars))]
        )[[1]]
      )
      if (j != min(J)) {
        x[[paste0("hazLessThan", j)]] <- rowSums(cbind(
          rep(0, nrow(x)),
          x[, paste0("Q", J[J < j], "Haz")]
        ))
        x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]] *
          (1 - x[[paste0("hazLessThan", j)]])
      } else {
        x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
        x[[paste0("hazLessThan", j)]] <- 0
      }
      return(x)
    }
    message("Finalizing hazard model...")
    dataList[[2]] <- get.pred.hazard( dataList[[2]],j = j,Qj_mod.0)
    dataList[[3]] <- get.pred.hazard( dataList[[3]],j = j,Qj_mod.1)
  }

  out <- list()
  out$dataList <- dataList
  out$ftimeMod <- list(ftimeMod.0,ftimeMod.1,get.pred.hazard)
  return(out)
}


estimateTreatment <- function(dat, adjustVars, glm.trt = NULL, SL.trt = NULL,
                              returnModels = FALSE, verbose = FALSE,
                              gtol = 1e-3, ...) {
  if (length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else {
    # binarize the outcome
    thisY <- as.numeric(dat$trt == max(dat$trt))

    # fit Super Learner
    if (!is.null(SL.trt)) {
      if (class(SL.trt) != "SuperLearner") {
        trtMod <- SuperLearner::mcSuperLearner(
          Y = thisY, X = adjustVars,
          newX = adjustVars,
          SL.library = SL.trt,
          id = dat$id, verbose = verbose,
          family = "binomial"
        )
      } else {
        trtMod <- SL.trt
      }
      dat[[paste0("g_", max(dat$trt))]] <- trtMod$SL.predict
      dat[[paste0("g_", min(dat$trt))]] <- 1 - trtMod$SL.predict
    } else if (!is.null(glm.trt) & is.null(SL.trt)) {
      # set up model formula and data for the treatment regression
      trt_form <- paste("thisY", "~", glm.trt, sep = " ")
      trt_data_in <- as.data.frame(cbind(adjustVars, thisY))

      # fit GLM if Super Learner not requested
      if (!("glm" %in% class(glm.trt)) & !("speedglm" %in% class(glm.trt))) {
        # fit the treatment model
        trtMod <- fast_glm(
          reg_form = stats::as.formula(trt_form),
          data = trt_data_in,
          family = stats::binomial()
        )
      } else {
        trtMod <- glm.trt
      }
      suppressWarnings(
        pred <- predict(trtMod, newdata = trt_data_in, type = "response")
      )
      dat[[paste0("g_", max(dat$trt))]] <- pred
      dat[[paste0("g_", min(dat$trt))]] <- 1 - pred
    }
  }

  # truncate propensities
  eval(parse(text = paste0(
    "dat$g_", min(dat$trt), "[dat$g_", min(dat$trt),
    "< gtol]<- gtol"
  )))
  eval(parse(text = paste0(
    "dat$g_", max(dat$trt), "[dat$g_", max(dat$trt),
    "< gtol]<- gtol"
  )))
  out <- list()
  out$dat <- dat
  out$trtMod <- NULL
  if (returnModels) out$trtMod <- trtMod
  return(out)
}



makeDataList <- function (dat, J, ntrt, uniqtrt, t0, bounds = NULL, ...)
{
  n <- nrow(dat)
  dataList <- vector(mode = "list", length = ntrt + 1)
  rankftime <- match(dat$ftime, sort(unique(dat$ftime)))
  dataList[[1]] <- dat[rep(1:nrow(dat), rankftime), ]
  for (j in J) {
    dataList[[1]][[paste0("N", j)]] <- 0
    dataList[[1]][[paste0("N", j)]][cumsum(rankftime)] <- as.numeric(dat$ftype ==
                                                                       j)
  }
  dataList[[1]]$C <- 0
  dataList[[1]]$C[cumsum(rankftime)] <- as.numeric(dat$ftype == 0)
  n.row.ii <- nrow(dataList[[1]])
  uniqftime <- unique(dat$ftime)
  orduniqftime <- uniqftime[order(uniqftime)]
  row.names(dataList[[1]])[row.names(dataList[[1]]) %in% paste(row.names(dat))] <- paste0(row.names(dat),
                                                                                          ".0")
  dataList[[1]]$t <- orduniqftime[as.numeric(paste(unlist(strsplit(row.names(dataList[[1]]),
                                                                   ".", fixed = TRUE))[seq(2, n.row.ii * 2, 2)])) + 1]
  if (!is.null(bounds)) {
    boundFormat <- data.frame(t = bounds$t)
    for (j in J) {
      if (paste("l", j, sep = "") %in% colnames(bounds)) {
        boundFormat[[paste0("l", j)]] <- bounds[, paste0("l",
                                                         j)]
      }
      else {
        boundFormat[[paste0("l", j)]] <- 0
      }
      if (paste("u", j, sep = "") %in% names(bounds)) {
        boundFormat[[paste0("u", j)]] <- bounds[, paste0("u",
                                                         j)]
      }
      else {
        boundFormat[[paste0("u", j)]] <- 1
      }
    }
    suppressMessages(dataList[[1]] <- plyr::join(x = dataList[[1]],
                                                 y = boundFormat, type = "left"))
    for (j in J) {
      tmp <- is.na(dataList[[1]][, paste0("l", j)])
      dataList[[1]][tmp, paste0("l", j)] <- 0
      tmp <- is.na(dataList[[1]][, paste0("u", j)])
      dataList[[1]][tmp, paste0("u", j)] <- 1
    }
  }
  else {
    for (j in J) {
      dataList[[1]][[paste0("l", j)]] <- 0
      dataList[[1]][[paste0("u", j)]] <- 1
    }
  }
  for (i in seq_len(ntrt)) {
    dataList[[i + 1]] <- dat[sort(rep(1:nrow(dat), t0)),
                             ]
    dataList[[i + 1]]$t <- rep(1:t0, n)
    for (j in J) {
      typejEvents <- dat$id[which(dat$ftype == j)]
      dataList[[i + 1]][[paste0("N", j)]] <- 0
      dataList[[i + 1]][[paste0("N", j)]][dataList[[i +
                                                      1]]$id %in% typejEvents & dataList[[i + 1]]$t >=
                                            dataList[[i + 1]]$ftime] <- 1
    }
    censEvents <- dat$id[which(dat$ftype == 0)]
    dataList[[i + 1]]$C <- 0
    dataList[[i + 1]]$C[dataList[[i + 1]]$id %in% censEvents &
                          dataList[[i + 1]]$t >= dataList[[i + 1]]$ftime] <- 1
    dataList[[i + 1]]$trt <- uniqtrt[i]
    dataList[[i + 1]]$ftime <- t0
    if (!is.null(bounds)) {
      suppressMessages(dataList[[i + 1]] <- plyr::join(x = dataList[[i +
                                                                       1]], y = boundFormat, type = "left"))
      for (j in J) {
        tmp <- is.na(dataList[[i + 1]][, paste0("l",
                                                j)])
        dataList[[i + 1]][tmp, paste0("l", j)] <- 0
        tmp <- is.na(dataList[[i + 1]][, paste0("u",
                                                j)])
        dataList[[i + 1]][tmp, paste0("u", j)] <- 1
      }
    }
    else {
      for (j in J) {
        dataList[[i + 1]][[paste0("l", j)]] <- .Machine$double.eps
        dataList[[i + 1]][[paste0("u", j)]] <- 1 - .Machine$double.eps
      }
    }
  }
  names(dataList) <- c("obs", uniqtrt)
  return(dataList)
}




predictTreatment <- function(dat, adjustVars, trtMod,gtol) {
  if (length(unique(dat$trt)) == 1) {
    eval(parse(text = paste0("dat$g_", unique(dat$trt), "<- 1")))
  } else {
    # binarize the outcome
    thisY <- as.numeric(dat$trt == max(dat$trt))
    if (!is.null(trtMod)) {
      suppressWarnings(pred <- predict(trtMod,dat))
      dat[[paste0("g_", max(dat$trt))]] <- pred$pred
      dat[[paste0("g_", min(dat$trt))]] <- 1 - pred$pred
    }
  }
  # truncate propensities
  eval(parse(text = paste0(
    "dat$g_", min(dat$trt), "[dat$g_", min(dat$trt),
    "< gtol]<- gtol"
  )))
  eval(parse(text = paste0(
    "dat$g_", max(dat$trt), "[dat$g_", max(dat$trt),
    "< gtol]<- gtol"
  )))
  return(dat)
}

predictHazards <- function(dataList,
                            J,
                            adjustVars,
                            ftimeMod,
                            env,
                            ...) {
  for (j in J) {
    # add all events less than current j to see who to include in regression
    NlessthanJ <- rep(0, nrow(dataList[[1]]))
    for (i in J[J < j]) {
      NlessthanJ <- NlessthanJ + dataList[[1]][[paste0("N", i)]]
    }

    Qj_mod.1 <- ftimeMod[[1]][[paste0("J", j)]]
    Qj_mod.0 <- ftimeMod[[2]][[paste0("J", j)]]

    # get predictions back
    get.pred.hazard <- function(x, j, Qj_mod) {
      suppressWarnings(
        x[[paste0("Q", j, "PseudoHaz")]] <- predict(
          Qj_mod,
          onlySL = TRUE,
          newdata = x[, c("t", names(adjustVars))]
        )[[1]]
      )
      if (j != min(J)) {
        x[[paste0("hazLessThan", j)]] <- rowSums(cbind(
          rep(0, nrow(x)),
          x[, paste0("Q", J[J < j], "Haz")]
        ))
        x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]] *
          (1 - x[[paste0("hazLessThan", j)]])
      } else {
        x[[paste0("Q", j, "Haz")]] <- x[[paste0("Q", j, "PseudoHaz")]]
        x[[paste0("hazLessThan", j)]] <- 0
      }
      return(x)
    }
    message("Predict hazard...")
    dataList[[2]] <- get.pred.hazard( dataList[[2]],j = j,Qj_mod.0)
    dataList[[3]] <- get.pred.hazard( dataList[[3]],j = j,Qj_mod.1)
  }
  return(dataList)
}

predictCensoring <- function(dataList,
                              adjustVars,
                              t0,
                              gtol = 1e-3,
                              env,
                              ctimeMod,
                              ...) {
  include <- !(dataList[[1]]$t == dataList[[1]]$ftime &
                 dataList[[1]]$C != 1 &
                 dataList[[1]]$t < t0)&
    !(dataList[[1]]$t == dataList[[1]]$ftime &
        dataList[[1]]$C == 1 &
        dataList[[1]]$t == t0)

  ctimeMod.1 <- ctimeMod[[1]]
  ctimeMod.0 <- ctimeMod[[2]]
  if (class(ctimeMod.1) != "noCens") {
    get.pred.censor <- function(x,ctimeMod) {
      g_dC <- rep(1, nrow(x))
      if (t0 != 1) {
        # temporarily replace time with t-1
        # NOTE: this will fail if t enters model as a factor
        x$t <- x$t - 1
        g_dC <-
          suppressWarnings(
            1 - predict(ctimeMod,newdata = x[, c("t",names(adjustVars))], onlySL = TRUE)[[1]]
          )
        # put time back to normal
        x$t <- x$t + 1
        # replace any observations with t = 1
        # to avoid extrapolation at t = 0
        g_dC[x$t == 1] <- 1
      }
      x$G_dC <- as.numeric(unlist(by(g_dC, x$id, FUN = cumprod)))
      return(x)
    }
    message("Predict censoring...")
    dataList[[2]] <- get.pred.censor(dataList[[2]] , ctimeMod.0)
    dataList[[3]] <- get.pred.censor(dataList[[3]] , ctimeMod.1)
  } else {
    dataList <- lapply(dataList, function(x) {
      x$G_dC <- 1
      x
    })
  }

  # truncate small propensities at gtol
  dataList$`0`$G_dC[dataList$`0`$G_dC < gtol] <- gtol
  dataList$`1`$G_dC[dataList$`1`$G_dC < gtol] <- gtol

  return(dataList)
}

