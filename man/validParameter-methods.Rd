\name{validParameter-methods}
\docType{methods}
\alias{validParameter-methods}
\alias{validParameter}
\alias{validParameter,ParamFamily-method}
\alias{validParameter,L2ScaleFamily-method}
\alias{validParameter,L2LocationFamily-method}
\alias{validParameter,L2LocationScaleFamily-method}
\alias{validParameter,BinomFamily-method}
\alias{validParameter,PoisFamily-method}
\alias{validParameter,GammaFamily-method}


\title{ Methods for function validParameter in Package `distrMod' }
\description{
  Methods for function \code{validParameter} in package \pkg{distrMod}
  to check whether a new parameter (e.g. "proposed" by an optimization)
  is valid.
}

\usage{
validParameter(object, \dots)
\S4method{validParameter}{ParamFamily}(object, param)
\S4method{validParameter}{L2ScaleFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{L2LocationFamily}(object, param)
\S4method{validParameter}{L2LocationScaleFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{BinomFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{PoisFamily}(object, param, tol=.Machine$double.eps)
\S4method{validParameter}{GammaFamily}(object, param, tol=.Machine$double.eps)
}

\details{
method for signature
\describe{
  \item{\code{ParamFamily}}{checks if all parameters are finite by \code{is.finite}
  if their length is between 1 and the joint length of main and nuisance parameter
  of \code{object}, and finally, if a call to \code{modifyParam(object)} with
  argument \code{param} would throw an error.
  }
  \item{\code{L2ScaleFamily}}{checks if the parameter is finite by \code{is.finite},
  if its length is 1, and if it is strictly larger than 0 (upto argument \code{tol})}
  \item{\code{L2LocationFamily}}{checks if the parameter is finite by \code{is.finite},
  if its length is 1}
  \item{\code{L2LocationScaleFamily}}{checks if the scale parameter is finite by \code{is.finite},
  if its length is 1 or 2 (e.g. if one features as nuisance parameter), and if it is 
  strictly larger than 0 (upto argument \code{tol})}
  \item{\code{BinomFamily}}{checks if the parameter is finite by \code{is.finite},
  if its length is 1, and if it is strictly larger than 0 and strictly smaller than 1 
  (upto argument \code{tol})}
  \item{\code{PoisFamily}}{checks if the parameter is finite by \code{is.finite},
  if its length is 1, and if it is strictly larger than 0 (upto argument \code{tol})}
  \item{\code{GammaFamily}}{checks if both parameters are finite by \code{is.finite},
  if their length is 1 or 2 (e.g. if one features as nuisance parameter), and if 
  both are strictly larger than 0 (upto argument \code{tol})}
}
}
\arguments{
  \item{object}{an object of class \code{ParamFamily}}
  \item{param}{either a numeric vector or an object of class 
  \code{ParamFamParameter}}
  \item{tol}{accuracy upto which the conditions have to be fulfilled}
  \item{\dots}{additional argument(s) for methods.}
}
\value{\code{logical} of length 1 --- valid or not}
\examples{
 NS <- NormLocationScaleFamily()
 validParameter(NS, c(scale=0.1, loc=2))
 validParameter(NS, c(scale=-0.1, loc=2))
 validParameter(NS, c(scale=0, loc=2))
 validParameter(NS, c(mean=2, sd=2))
}
\keyword{models}
