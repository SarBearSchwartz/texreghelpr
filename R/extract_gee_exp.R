extract_gee_exp <- function (model,
                             robust = TRUE,
                             include.dispersion = TRUE,
                             include.nobs = TRUE,
                             ...){
  s <- summary(model, ...)
  names <- rownames(coef(s))
  co <- coef(s)[, 1]
  co_exp <- exp(co)

  if (robust == TRUE) {
    se <- coef(s)[, 4]
    zval <- coef(s)[, 5]
  }
  else {
    se <- coef(s)[, 2]
    zval <- coef(s)[, 3]
  }

  ci.l <- exp(co - 1.96*se)
  ci.u <- exp(co + 1.96*se)

  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
  n <- nobs(model)
  disp <- s$scale
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()

  if (include.dispersion == TRUE) {
    gof <- c(gof, disp)
    gof.names <- c(gof.names, "Dispersion")
    gof.decimal <- c(gof.decimal, TRUE)
  }

  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num. obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }

  tr <- createTexreg(coef.names = names,
                     coef = co_exp,
                     ci.low = ci.l,
                     ci.up = ci.u,
                     pvalues = pval,
                     gof.names = gof.names,
                     gof = gof,
                     gof.decimal = gof.decimal)
  return(tr)
}
