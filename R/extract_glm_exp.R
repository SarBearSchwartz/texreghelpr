extract_glm_exp <- function (model,
                             include.any = FALSE,
                             include.aic = FALSE,
                             include.bic = FALSE,
                             include.loglik = FALSE,
                             include.deviance = FALSE,
                             include.nobs = FALSE,
                             ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$coef)
  coefficients <- s$coef[, 1]
  co_exp = exp(coefficients)

  ci.l <- exp(confint(model)[,1])
  ci.u <- exp(confint(model)[,2])

  standard.errors <- s$coef[, 2]
  significance <- s$coef[, 4]

  aic <- AIC(model)
  bic <- BIC(model)
  lik <- logLik(model)[1]
  dev <- deviance(model)
  n <- nobs(model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()


  if(include.any == FALSE){
    include.aic = FALSE
    include.bic = FALSE
    include.loglik = FALSE
    include.deviance = FALSE
    include.nobs = FALSE
  }
  if (include.aic == TRUE) {
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(coef.names = coefficient.names,
                     coef = co_exp,
                     ci.low = ci.l,
                     ci.up = ci.u,
                     pvalues = significance,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}
