\name{trafo-methods}
\docType{methods}
\alias{trafo-methods}
\alias{trafo}
\alias{trafo,Estimate,missing-method}
\alias{trafo,Estimate,ParamFamParameter-method}
\alias{trafo,ParamFamParameter,missing-method}
\alias{trafo,ParamFamily,missing-method}
\alias{trafo,ParamFamily,ParamFamParameter-method}
\alias{trafo.fct}
\alias{trafo.fct-methods}
\alias{trafo.fct,ParamFamily-method}
\alias{trafo<-}
\alias{trafo<-,ParamFamParameter-method}
\alias{trafo<-,ParamFamily-method}
\title{ Methods for function trafo in Package `distrMod' }
\description{
  Methods for function \code{trafo} in package \pkg{distrMod};
  there are accessor (\code{trafo}) and replacement (\code{trafo<-})
  versions.
 }
\usage{
trafo(object, param, ...)
\S4method{trafo}{Estimate,missing}(object,param)
\S4method{trafo}{ParamFamParameter,missing}(object,param)
\S4method{trafo}{ParamFamily,missing}(object,param)
\S4method{trafo}{ParamFamily,ParamFamParameter}(object,param)
\S4method{trafo}{Estimate,ParamFamParameter}(object,param)
trafo.fct(object)
trafo(object) <- value
}
\arguments{
  \item{object}{an object of either class \code{Estimate},
  \code{ParamFamParameter}, \code{ParamFamily}}
  \item{param}{an object of class \code{ParamFamParameter}; the parameter
  value at which to evaluate the transformation}
  \item{value}{a matrix or a function; if it is a matrix, dimensions must
  be consistent to the parametric setting; if it is function, it should
  take one argument \code{param} of class \code{ParamFamParameter} and
  return a list of length two with named components \code{fval}
  (the function value, see below)
  and \code{mat} (a matrix --- with the same dimensions consistency
  conditions as above).}
  \item{\dots}{additional argument(s) for methods; not used so far.}
}
\value{
  The return value depends on the signature.
  For \code{trafo.fct}, we return the corresponding function 
  \eqn{\tau()}{tau()} (see below).
  For \code{trafo}, we have:
  \item{signature \code{Estimate,missing}:}{a list of length two with components
  \code{fct} and \code{mat} (see below)}
  \item{signature \code{Estimate,ParamFamParameter}:}{a list of length two with components
  \code{fct} and \code{mat} (see below)}
  \item{signature \code{ParamFamParameter,missing}:}{a matrix (see below)}
  \item{signature \code{ParamFamily,missing}:}{a matrix (see below)}
  \item{signature \code{ParamFamily,ParamFamParameter}:}{a list of length two
  with components \code{fct} and \code{mat} (see below)}
  }
\details{
  \code{trafo} is a slot of class \code{ParamFamParameter}, which
  in turn is a slot of class \code{ParamFamily}. It also sort of
  arises in class \code{Estimate}, i.e., all slots can be identified
  by the information contained in an instance thereof.

  \code{trafo} realizes partial influence curves; i.e.; we are only
  interested is some possibly lower dimensional smooth (not necessarily
  linear or even coordinate-wise) aspect/transformation \eqn{\tau}{tau}
  of the parameter \eqn{\theta}{theta}.

  To be coherent with the corresponding \emph{nuisance}
  implementation, we make the following convention:

  The full parameter \eqn{\theta}{theta} is split up coordinate-wise
  in a main parameter \eqn{\theta'}{theta'} and a nuisance parameter
  \eqn{\theta''}{theta''} (which is unknown, too, hence has to be
  estimated, but only is of secondary interest) and a fixed,
  known part \eqn{\theta'''}{theta'''}.
  
  Without loss of generality, we restrict ourselves to the case that
  transformation \eqn{\tau}{tau} only acts on the main parameter
  \eqn{\theta'}{theta'} --- if we want to transform the whole
  parameter, we only have to assume that both nuisance parameter
  \eqn{\theta''}{theta''} and fixed, known part of the parameter
  \eqn{\theta'''}{theta'''} have length 0.
  
  To the implementation:
  
  Slot \code{trafo} can either contain a (constant) matrix
  \eqn{D_\theta}{D_theta} or a function
  \deqn{\tau\colon \Theta' \to \tilde \Theta,\qquad \theta \mapsto \tau(\theta)}{tau:  Theta' -> TTheta,  theta |-> tau(theta)} 
  mapping main parameter
  \eqn{\theta'}{theta'} to some range \eqn{\tilde \Theta}{TTheta}.

  If \emph{slot value} \code{trafo} is a function, besides \eqn{\tau(\theta)}{tau(theta)}, 
  it will also return the corresponding derivative matrix
  \eqn{\frac{\partial}{\partial \theta}\tau(\theta)}{(d/d theta) (tau(theta))}.
  More specifically, the return value of this function \code{theta} is a
  list with entries \code{fval}, the function value \eqn{\tau(\theta)}{tau(theta)},
  and \code{mat}, the derivative matrix.

  In case \code{trafo} is a matrix \eqn{D}, we interpret it as such a derivative
  matrix \eqn{\frac{\partial}{\partial \theta}\tau(\theta)}{(d/d theta) (tau(theta))},
  and, correspondingly, \eqn{\tau(\theta)}{tau(theta)} as the linear mapping
  \eqn{\tau(\theta)=D\,\theta}{tau(theta)=D * theta}.\cr

  According to the signature, \emph{method} \code{trafo} will return different
  return value types. For signature
  \describe{
    \item{\code{Estimate,missing}:}{it will return a list with entries
    \code{fct}, the function  \eqn{\tau}{tau}, and \code{mat}, the matrix
    \eqn{\frac{\partial}{\partial \theta}\tau(\theta)}{(d/d theta) (tau(theta))}.
    function \eqn{\tau}{tau} will then return the list \code{list(fval, mat)}
    mentioned above.
    }
    \item{\code{Estimate,ParamFamParameter}:}{as signature
    \code{Estimate,missing}.}
    \item{\code{ParamFamParameter,missing}:}{it will just return the
     corresponding matrix.}
    \item{\code{ParamFamily,missing}:}{is just wrapper to signature
    \code{ParamFamParameter,missing}.}
    \item{\code{ParamFamily,ParamFamParameter}:}{as signature
    \code{Estimate,missing}.}
   }
  }

\examples{
## Gaussian location and scale
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

## Maximum likelihood estimator
(res <- MLEstimator(x = x, ParamFamily = NS))
## confidence interval
 confint(res)

}
\keyword{models}
