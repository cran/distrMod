\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ShowDetails-method}

\title{Common `print' Methods for S4 classes in Package `distrMod' }
\description{
  Methods for \code{print} to the S4 classes in package \pkg{distrMod};
 }
\usage{
\S4method{print}{ShowDetails}(x, digits = getOption("digits"),
                show.details = c("maximal", "minimal", "medium"))
}
\arguments{
  \item{x}{object of class \code{ShowDetails}, a class union of
  classes \code{OptionalNumeric}, \code{OptionalMatrix},
 \code{MatrixorFunction}, \code{Estimate}, \code{MCEstimate}.}
  \item{digits}{unchanged w.r.t. default method of package \code{base}:
   a non-null value for 'digits' specifies the minimum number of
          significant digits to be printed in values.  The default,
          'NULL', uses 'getOption(digits)'.  (For the interpretation
          for complex numbers see 'signif'.) Non-integer values will be
          rounded down, and only values greater than or equal to 1 and
          no greater than 22 are accepted.}
  \item{show.details}{a character, controlling the degree of detailedness
  of the output; currently the following values are permitted:
  \code{"maximal"}, \code{"minimal"}, \code{"medium"}; for the meaning
  for the actual class, confer to the corresponding class help file.}
%  \item{\dots}{further arguments to be passed to or from other methods.}
}
\details{
      This method provides sort of a ''\code{show}
       with extra arguments'', in form of a common \code{print} method
       for the mentioned S4 classes. Essentially this \code{print} method
       just temporarily sets the global options according to the optional arguments
       \code{digits} and \code{show.details}, calls \code{show} and then re-sets
       the options to their global settings.}

\examples{
## set options to maximal detailedness
show.old <- getdistrModOption("show.details")
distrModoptions("show.details" = "maximal")
## define a model
NS <- NormLocationScaleFamily(mean=2, sd=3)
## generate data out of this situation
x <- r(distribution(NS))(30)

## want to estimate mu/sigma, sigma^2
## -> new trafo slot:
trafo(NS) <- function(param){
  mu <- param["mean"]
  sd <- param["sd"]
  fval <- c(mu/sd, sd^2)
  nfval <- c("mu/sig", "sig^2")
  names(fval) <- nfval
  mat <- matrix(c(1/sd,0,-mu/sd^2,2*sd),2,2)
  dimnames(mat) <- list(nfval,c("mean","sd"))
  return(list(fval=fval, mat=mat))
}
print(param(NS))
print(param(NS), show.details = "minimal")
print(param(NS), show.details = "medium")
## Maximum likelihood estimator
res <- MLEstimator(x = x, ParamFamily = NS)
print(res) #equivalent to 'show(res)' or 'res'
print(res, digits = 4)
print(res, show.details = "minimal")
print(res, show.details = "medium")
distrModoptions("show.details" = show.old)
}
\keyword{models}
