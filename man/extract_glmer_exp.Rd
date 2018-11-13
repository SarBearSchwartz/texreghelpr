\name{extract_glmer_exp}
\alias{extract_glmer_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract EXPONENTIATED coefficients and GOF measures from a `glmer` object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_glmer_exp(model,
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
                   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  the name of a model fit via the `lme4::glmer()` command
  }
  \item{method}{
  The method used to compute confidence intervals or p values. In lme4 models, the default value "naive" computes naive p values while the other methods compute confidence intervals using the confint function.
  }
  \item{include.aic}{
   should Akaike's information criterion (AIC) be reported?
  }
  \item{include.bic}{
   should Bayesian's information criterion (BIC) be reported?
  }
  \item{include.dic}{
   should the deviance information criterion (DIC) be reported?
  }
  \item{include.loglik}{
   should the log-likelihood be reported?
  }
  \item{include.deviance}{
   should the deviance be reported?
  }
  \item{include.nobs}{
   should the number of observations be reported?
  }
  \item{include.groups}{
   should the number of groups be reported?
  }
  \item{\dots}{
  Custom parameters which are handed over to subroutines.
  }
}
\details{
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
www.sarahschwartzstats.com
}
\author{
Sarah Schwartz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
