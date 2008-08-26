\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint}
\alias{confint,ANY,missing-method}
\alias{confint,Estimate,missing-method}
\alias{confint,mle,missing-method}
\alias{confint,profile.mle,missing-method}
\title{ Methods for function confint in Package `distrMod' }
\description{
  Methods for function \code{confint} in package \pkg{distrMod};
  by default uses \code{\link[stats]{confint}} and its corresponding
  S3-methods, but also computes
  (asymptotic) confidence intervals for objects of class \code{Estimate}.
  Computes confidence intervals for one or more parameters in a fitted
  model.
}
\usage{
confint(object, method, \dots)
\S4method{confint}{ANY,missing}(object, method, parm, level = 0.95, \dots)
\S4method{confint}{Estimate,missing}(object, method, level = 0.95)
\S4method{confint}{mle,missing}(object, method, parm, level = 0.95, \dots)
\S4method{confint}{profile.mle,missing}(object, method, parm, level = 0.95, \dots)
}
\arguments{
  \item{object}{in default / signature \code{ANY} case: a fitted model object,
  in signature \code{Estimate} case, an object of class \code{Estimate}}
  \item{parm}{only used in default / signature \code{ANY} case: 
    a specification of which parameters are to be given
    confidence intervals, either a vector of numbers or a vector of
    names.  If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
  \item{method}{not yet used (only as \code{missing}; later to allow
                for various methods}
  \item{\dots}{additional argument(s) for methods.}
}
\value{
  \itemize{
  \item{signature \code{ANY,missing}:}{
  A matrix (or vector) with columns giving lower and upper confidence
  limits for each parameter. These will be labelled as (1-level)/2 and
  1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).}
  \item{signature \code{Estimate,missing}:}{An object of class \code{Confint}}
  }
}
\details{
  \code{confint} is a generic function.  Its behavior differs according to
  its arguments.
  \itemize{
   \item{signature \code{ANY,missing}:}{the default method;
   uses the S3 generic of package \pkg{stats},
    see \code{\link[stats]{confint}}.}
    \item{signature \code{Estimate,missing}:}{ 
  will return a corresponding confidence interval assuming
  asymptotic normality, and hence needs suitably filled slot \code{asvar}
  in argument \code{object}. Besides the actual bounds, organized in 
  an array just as in the S3 generic, the return value also captures 
  the name of the estimator for which it is produced, as well as the corresponding
  call producing the estimator, and the corresponding 
  \code{trafo} and \code{nuisance} slots/parts.}
  }
}
\seealso{
  \code{\link[stats]{confint}},
  \code{\link[MASS]{confint.glm}} and \code{\link[MASS]{confint.nls}} in
  package \pkg{MASS}, \code{\link{Confint-class}}.
}
\examples{
## for signature ANY examples confer stats::confint
## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Maximum likelihood estimator
res <- MLEstimator(x = x, ParamFamily = G)
confint(res)

### for comparison:
require(MASS)
(res1 <- fitdistr(x, "gamma"))
## add a convenient (albeit wrong)
## S3-method for vcov:
## --- wrong as in general cov-matrix
##     will not be diagonal
## but for conf-interval this does
## not matter...
vcov.fitdistr <- function(object, ...){
     v<-diag(object$sd^2)
     rownames(v) <- colnames(v) <- names(object$estimate) 
     v}

## explicitely transforming to
## MASS parametrization:
mtrafo <- function(x){
     nms0 <- names(c(main(param(G)),nuisance(param(G))))
     nms <- c("shape","rate")
     fval0 <- c(x[2], 1/x[1])
     names(fval0) <- nms
     mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
                     dimnames = list(nms,nms0))                          
     list(fval = fval0, mat = mat0)}

G2 <- G
trafo(G2) <- mtrafo
res2 <- MLEstimator(x = x, ParamFamily = G2)

old<-getdistrModOption("show.details")
distrModoptions("show.details" = "minimal")
res
res1
res2
confint(res)
confint(res1)
confint(res2)
confint(res,level=0.99)
distrModoptions("show.details" = old)
 
}
\keyword{models}
