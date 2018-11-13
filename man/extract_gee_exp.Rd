\name{extract_gee_exp}
\alias{extract_gee_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract EXPONENTIATED coefficients and GOF measures from a `gee` object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract_gee_exp(model,
                robust = TRUE,
                include.dispersion = TRUE,
                include.nobs = TRUE,
                ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  Rhe name of a model fit via the `gee::gee()` function
}
  \item{robust}{
  Report robust (sandwich) instead of naive standard errors?
}
  \item{include.dispersion}{
  Should the dispersion or scale parameter be reported?
}
  \item{include.nobs}{
  Should the number of observations be reported?
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
\examples{
library(HSAUR)        # package with the dataset
library(gee)          # package to fit the GEE model
library(texreg)       # package to create the table
library(texreghelpr)  # package to exp the model parameters

## Read in the dataset
data("epilepsy", package = "HSAUR")

## Fit a GEE model
model <- gee::gee(seizure.rate ~ base + age + treatment,
                  data = epilepsy,
                  family = poisson(link = "log"),
                  id = subject,
                  corstr = "exchangeable")

## Ugly summary of the model
summary(model)

## Nice table for the model parameters
texreg::screenreg(texreghelpr::extract_gee_exp(model),
                  single.row = TRUE,
                  ci.test = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
