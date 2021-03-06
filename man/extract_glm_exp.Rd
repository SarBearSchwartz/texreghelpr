\name{extract_glm_exp}
\alias{extract_glm_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract EXPONENTIATED coefficients and GOF measures from a `glm` object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_glm_exp(model,
               include.aic = TRUE,
               include.bic = TRUE,
               include.loglik = TRUE,
               include.deviance = TRUE,
               include.nobs = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  the name of a model fit via the `stats::glm()` command
  }
  \item{include.aic}{
   should Akaike's information criterion (AIC) be reported?
  }
  \item{include.bic}{
   should Bayesian's information criterion (BIC) be reported?
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
  \item{\dots}{
  Custom parameters which are handed over to subroutines.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\examples{library(HSAUR)        # package with the dataset
library(gee)          # package to fit the GEE model
library(texreg)       # package to create the table
library(texreghelpr)  # package to exp the model parameters

## Read in the dataset
data("plasma", package = "HSAUR")

## Fit a GLM models
model_1 <- stats::glm(ESR ~ fibrinogen ,
                      data = plasma,
                      family = binomial(link = "logit"))

model_2 <- stats::glm(ESR ~ fibrinogen + globulin,
                    data = plasma,
                    family = binomial(link = "logit"))

## Ugly summary of the model
summary(model_1)

## Nice table for the model parameters
texreg::screenreg(list(texreghelpr::extract_glm_exp(model_1),
                       texreghelpr::extract_glm_exp(model_2)),
                  single.row = TRUE,
                  ci.test = 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
