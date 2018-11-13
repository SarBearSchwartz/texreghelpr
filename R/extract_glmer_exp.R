extract_glmer_exp <- function (model,
                               method = c("naive", "profile", "boot", "Wald"),
                               level = 0.95,
                               nsim = 1000,
                               include.aic = TRUE,
                               include.bic = TRUE,
                               include.dic = FALSE,
                               include.deviance = FALSE,
                               include.loglik = TRUE,
                               include.nobs = TRUE,
                               include.groups = TRUE,
                               include.variance = TRUE,
                               ...) {
  if (packageVersion("lme4") < 1) {
    message("Please update to a newer 'lme4' version for full compatibility.")
  }
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.aic == TRUE) {
    aic <- AIC(model)
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- BIC(model)
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.dic == TRUE) {
    is_REML <- lme4::isREML(model)
    llik <- logLik(model, REML = is_REML)
    dev <- deviance(lme4::refitML(model))
    n <- lme4::getME(model, "devcomp")$dims["n"]
    Dhat <- -2 * (llik)
    pD <- dev - Dhat
    DIC <- dev + pD[[1]]
    gof <- c(gof, DIC)
    gof.names <- c(gof.names, "DIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.deviance == TRUE) {
    dev <- deviance(lme4::refitML(model))
    gof <- c(gof, dev)
    gof.names <- c(gof.names, "Deviance")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.loglik == TRUE) {
    lik <- logLik(model)[1]
    gof <- c(gof, lik)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    n <- dim(model.frame(model))[1]
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (include.groups == TRUE) {
    grps <- sapply(model@flist, function(x) length(levels(x)))
    grp.names <- names(grps)
    grp.names <- paste("Num. groups:", grp.names)
    gof <- c(gof, grps)
    gof.names <- c(gof.names, grp.names)
    gof.decimal <- c(gof.decimal, rep(FALSE, length(grps)))
  }
  if (include.variance == TRUE) {
    vc <- as.data.frame(lme4::VarCorr(model))
    for (i in 1:nrow(vc)) {
      if (is.na(vc[i, 2]) && is.na(vc[i, 3])) {
        gof.names <- c(gof.names, "Var: Residual")
      }
      else if (is.na(vc[i, 3])) {
        gof.names <- c(gof.names, paste("Var:", vc[i,
                                                   1], vc[i, 2]))
      }
      else {
        gof.names <- c(gof.names, paste("Cov:", vc[i,
                                                   1], vc[i, 2], vc[i, 3]))
      }
      gof <- c(gof, vc[i, 4])
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }

  betas <- lme4::fixef(model, ...)


  if ("confint.merMod" %in% methods("confint") && method[1] !=
      "naive") {
    ci <- tryCatch({
      ci <- confint(model, method = method[1], level = level,
                    nsim = nsim, ...)
    }, error = function(err) {
      method <- "naive"
      message(paste("Confidence intervals not available for",
                    "this model. Using naive p values instead."))
    })
    if (is.null(ci)) {
      method <- "naive"
    }
    else {
      last <- nrow(ci)
      number <- length(betas)
      first <- last - number + 1
      ci <- ci[first:last, ]
      if (class(ci) == "matrix") {
        ci.l <- ci[, 1]
        ci.u <- ci[, 2]
      }
      else {
        ci.l <- ci[1]
        ci.u <- ci[2]
      }
    }
  }
  else if (method[1] != "naive") {
    method[1] <- "naive"
    message(paste("confint.merMod method not found. Using naive p values",
                  "instead."))
  }


  if (method[1] == "naive") {
    Vcov <- tryCatch({
      Vcov <- vcov(model, useScale = FALSE, ...)
    }, error = function(err) {
      stop(paste("Please load the Matrix package or update to the latest",
                 "development version of lme4 and run this command again."))
    })
    Vcov <- as.matrix(Vcov)
    se <- sqrt(diag(Vcov))
    zval <- betas/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    tr <- createTexreg(coef.names = names(betas),
                       coef = exp(betas),
                       ci.low = exp(betas - 1.98*se),
                       ci.up  = exp(betas + 1.98*se),
                       pvalues = pval,
                       gof.names = gof.names,
                       gof = gof,
                       gof.decimal = gof.decimal)
  }
  else {
    tr <- createTexreg(coef.names = names(betas),
                       coef   = exp(betas),
                       ci.low = exp(ci.l),
                       ci.up  = exp(ci.u),
                       gof.names = gof.names,
                       gof = gof,
                       gof.decimal = gof.decimal)
  }
  return(tr)
}
