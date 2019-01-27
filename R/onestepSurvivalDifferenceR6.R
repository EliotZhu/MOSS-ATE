require("R6")
require("SuperLearner")

#' onestep TMLE of difference survival curve
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' # MOSS_difference$new(dat, epsilon.step = 1e-5, max.iter = 1e3, tol = 1/nrow(dat), T.cutoff = NULL, verbose = FALSE)
#' @field dat data.frame
#' @field epsilon.step float
#' @field max.iter int
#' @field tol float
#' @field T.cutoff int
#' @field verbose bool
#' @section Methods:
#' @export
MOSS_difference <- R6Class("MOSS_difference",
  public = list(
    iticount = NULL,
    dat = NULL,
    epsilon.step = NULL,
    max.iter = NULL,
    tol = NULL,
    T.cutoff = NULL,
    verbose = NULL,
    ftimeMod = NULL,
    trtMod = NULL,
    ctimeMod  = NULL,

# single MOSS object
    MOSS_A1 = NULL,
    MOSS_A0 = NULL,

# targeting
    D1.t = NULL,
    Pn.D1.t = NULL,
    stopping_history = numeric(),
    Psi.hat = NULL,
    sd_EIC = NULL,
    upper_CI = NULL,
    lower_CI = NULL,
    # simultaneous CI
    simul_CI = NULL,
    initialize = function(dat,
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          T.cutoff = NULL,
                          verbose = FALSE,
                          ftimeMod,
                          trtMod,
                          ctimeMod) {
      self$ftimeMod <- ftimeMod
      self$trtMod <- trtMod
      self$ctimeMod <- ctimeMod

      self$MOSS_A1 <- MOSSATE::MOSS$new(data_out, dW = 0,
                                         epsilon.step = epsilon.step,
                                         max.iter = max.iter,pred_data = T,
                                         ftimeMod = self$ftimeMod,
                                         ctimeMod = self$ctimeMod,
                                         trtMod = self$trtMod)
      self$MOSS_A0 <- MOSSATE::MOSS$new(data_out, dW = 1,
                                        epsilon.step = epsilon.step,
                                        max.iter = max.iter,pred_data = T,
                                        ftimeMod = self$ftimeMod,
                                        ctimeMod = self$ctimeMod,
                                        trtMod = self$trtMod)

      self$dat <- self$MOSS_A1$dat
      self$epsilon.step <- self$MOSS_A1$epsilon.step
      self$max.iter <- self$MOSS_A1$max.iter
      self$tol <- self$MOSS_A1$tol
      self$T.cutoff <- self$MOSS_A1$T.cutoff
      self$verbose <- self$MOSS_A1$verbose
    },

    initial_fit = function(env){

      self$MOSS_A1$surv.predict(env = env)
      self$MOSS_A1$transform_failure_hazard_to_survival()
      self$MOSS_A1$transform_failure_hazard_to_pdf()
      self$MOSS_A1$compute_EIC()

      self$MOSS_A0$surv.predict(env = env)
      self$MOSS_A0$transform_failure_hazard_to_survival()
      self$MOSS_A0$transform_failure_hazard_to_pdf()
      self$MOSS_A0$compute_EIC()
    },

    onestep_diff_curve = function(){
      self$D1.t <- self$MOSS_A1$D1.t - self$MOSS_A0$D1.t
      self$Pn.D1.t <- colMeans(self$D1.t)
      self$MOSS_A1$D1.t <- self$D1.t
      self$MOSS_A0$D1.t <- self$D1.t
      self$MOSS_A1$Pn.D1.t <- self$Pn.D1.t
      self$MOSS_A0$Pn.D1.t <- self$Pn.D1.t

      iter_count <- 0
      stopping_prev <- Inf

      stopping <- self$MOSS_A1$compute_stopping()
      while ((stopping >= self$tol) & (iter_count <= self$max.iter)) {
      # while ((stopping >= self$tol) & (iter_count <= self$max.iter) & ((stopping_prev - stopping) >= max(-self$tol, -1e-5))) {
        print(stopping)
        if (stopping_prev < stopping) onestepfit$epsilon.step <- -onestepfit$epsilon.step
        self$MOSS_A1$onestep_curve_update_mat()
        self$MOSS_A0$onestep_curve_update_mat()
        self$MOSS_A1$compute_EIC()
        self$MOSS_A0$compute_EIC()

        self$D1.t <- self$MOSS_A1$D1.t - self$MOSS_A0$D1.t
        self$Pn.D1.t <- colMeans(self$D1.t)
        self$MOSS_A1$D1.t <- self$D1.t
        self$MOSS_A0$D1.t <- self$D1.t
        self$MOSS_A1$Pn.D1.t <- self$Pn.D1.t
        self$MOSS_A0$Pn.D1.t <- self$Pn.D1.t

        iter_count <- iter_count + 1
        self$stopping_history[iter_count] <- stopping
        stopping_prev <- self$stopping_history[iter_count]
        stopping <- self$MOSS_A1$compute_stopping()
        if(iter_count>2){if(stopping>=stopping_prev) {break}}

        # if (iter_count %% 10 == 0) onestepfit$plot_onestep_curve(add = TRUE)
        # if (iter_count %% 10 == 0) plot(onestepfit$Pn.D1.t); abline(h = 0)
      }

      if (iter_count == self$max.iter) {
        warning('Max Iter count reached, stop iteration.')
      }

      self$MOSS_A1$compute_Psi()
      self$MOSS_A0$compute_Psi()
      self$Psi.hat <- self$MOSS_A1$Psi.hat - self$MOSS_A0$Psi.hat

      self$sd_EIC <- rep(NA, self$MOSS_A1$T.max)
      self$sd_EIC[self$MOSS_A1$T.uniq] <- apply(self$D1.t, 2, sd)
      self$sd_EIC <- zoo::na.locf(self$sd_EIC)
      self$upper_CI <- self$Psi.hat + 1.96 * self$sd_EIC/sqrt(self$MOSS_A1$n_sample)
      self$lower_CI <- self$Psi.hat - 1.96 * self$sd_EIC/sqrt(self$MOSS_A1$n_sample)
    },
    print = function(){
      out <- data.frame(1:self$MOSS_A1$T.max, self$Psi.hat, self$sd_EIC, self$upper_CI, self$lower_CI)
      colnames(out) <- c('Time', 'survival curve', 'std_err', 'upper_CI', 'lower_CI')
      return(out)
    },
    plot_onestep_curve = function(...){
      step_curve <- stepfun(x = 1:self$MOSS_A1$T.max, y = c(self$Psi.hat, self$Psi.hat[length(self$Psi.hat)]))
      # can `add`, `col`
      curve(step_curve, from = 0, to = self$MOSS_A1$T.max, ...)
    },
    plot_CI_pointwise = function(...){
      self$plot_onestep_curve(...)
      polygon(c(1:self$MOSS_A1$T.max, rev(1:self$MOSS_A1$T.max)), c(c(self$upper_CI), rev(c(self$lower_CI))),
                          col = rgb(0.7,0.7,0.7,0.4),
                          border = NA,
                          ...)
    },
    compute_CI_simultaneous = function(){
      Sigma_hat_EIC <- cor(self$D1.t)
      Sigma_hat_EIC[is.na(Sigma_hat_EIC)] <- 1e-10 * rnorm(n = sum(is.na(Sigma_hat_EIC))) # fill in where var is 0
      q_95_simCI <- simCI_quant(Sigma_hat_EIC, B = 500)
      CI_mat <- matrix(NA, nrow = self$MOSS_A1$T.max, ncol = 2)
      for (i in 1:self$MOSS_A1$T.max) {
        CI_mat[i,] <- ConfInt(est = self$Psi.hat[i], q = q_95_simCI, sd_est = self$sd_EIC[i]/sqrt(self$MOSS_A1$n_sample))
      }
      simul_CI <- data.frame(CI_mat)
      colnames(simul_CI) <- c('lower_CI', 'upper_CI')
      self$simul_CI <- simul_CI
    },
    plot_CI_simultaneous = function(...){
      self$plot_onestep_curve(...)
      polygon(c(1:self$MOSS_A1$T.max, rev(1:self$MOSS_A1$T.max)), c(c(self$simul_CI[,'upper_CI']), rev(c(self$simul_CI[,'lower_CI']))),
              col = rgb(0.7,0.7,0.7,0.4),
              border = NA,
              ...)
    }

  )
)
