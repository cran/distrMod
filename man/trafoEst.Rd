\name{trafoEst}
\alias{trafoEst}

\title{ Function trafoEst in Package `distrMod' }
\description{
  \code{trafoEst} takes a \eqn{\tau}{tau} like function (compare 
  \code{\link{trafo-methods}}) and transforms an existing estimator by means 
  of this transformation.
 }
\usage{
trafoEst(fct, estimator)
}
\arguments{
  \item{fct}{a \eqn{\tau}{tau} like function, i.e., a function
  in the main part \eqn{\theta}{theta} of the parameter returning a list \code{list(fval, mat)}
  where  \code{fval} is the function value \eqn{\tau(\theta)}{tau(theta)}
  of the transformation, and \code{mat}, its derivative matrix at
  \eqn{\theta}{theta}.}
  \item{estimator}{an object of class \code{Estimator}.}
}
\value{
  exactly the argument \code{estimator}, but with modified slots
  \code{estimate},   \code{asvar}, and   \code{trafo}.
  }  
\details{
The disadvantage of this proceeding is that the transformation is
not accounted for in determining the estimate (e.g. in a corresponding
optimality); it simply transforms an existing estimator, without
reapplying it to data. This becomes important in optimally robust
estimation.
}

\examples{
## Gaussian location and scale
NS <- NormLocationScaleFamily(mean=2, sd=3)
## generate data out of this situation
x <- r(distribution(NS))(30)

## want to estimate mu/sigma, sigma^2
## -> without new trafo slot:
mtrafo <- function(param){
  mu <- param["mean"]
  sd <- param["sd"]
  fval <- c(mu/sd, sd^2)
  nfval <- c("mu/sig", "sig^2")
  names(fval) <- nfval
  mat <- matrix(c(1/sd,0,-mu/sd^2,2*sd),2,2)
  dimnames(mat) <- list(nfval,c("mean","sd"))
  return(list(fval=fval, mat=mat))
}

## Maximum likelihood estimator in the original problem
res0 <- MLEstimator(x = x, ParamFamily = NS)
## transformation
res <- trafoEst(mtrafo, res0)
## confidence interval
 confint(res)
}
\keyword{models}
