\name{mceCalc-methods}
\docType{methods}
\alias{mceCalc-methods}
\alias{mleCalc-methods}
\alias{mceCalc}
\alias{mleCalc}
\alias{mceCalc,numeric,ParamFamily-method}
\alias{mleCalc,numeric,ParamFamily-method}
\alias{mleCalc,numeric,BinomFamily-method}
\alias{mleCalc,numeric,PoisFamily-method}
\alias{mleCalc,numeric,NormLocationFamily-method}
\alias{mleCalc,numeric,NormScaleFamily-method}
\alias{mleCalc,numeric,NormLocationScaleFamily-method}
\title{ Methods for functions mceCalc and mleCalc in Package `distrMod' }
\description{
  Methods for functions \code{mceCalc} and \code{mleCalc} in package \pkg{distrMod};
 }
\usage{
mceCalc(x, PFam, ...)
mleCalc(x, PFam, ...)
\S4method{mceCalc}{numeric,ParamFamily}(x, PFam, criterion, 
                   startPar = NULL, penalty = 1e20, crit.name,
                   Infos = NULL, withthetaPar = FALSE,\dots)
\S4method{mleCalc}{numeric,ParamFamily}(x, PFam, startPar = NULL, 
                   penalty = 1e20, Infos = NULL, \dots)
\S4method{mleCalc}{numeric,BinomFamily}(x, PFam, \dots)
\S4method{mleCalc}{numeric,PoisFamily}(x, PFam, \dots)
\S4method{mleCalc}{numeric,NormLocationFamily}(x, PFam, \dots)
\S4method{mleCalc}{numeric,NormScaleFamily}(x, PFam, \dots)
\S4method{mleCalc}{numeric,NormLocationScaleFamily}(x, PFam, \dots)
}
\arguments{
  \item{x}{numeric; data at which to evaluate the estimator}
  \item{PFam}{an object of class \code{ParamFamily}; the parametric family at
              which to evaluate the estimator}
  \item{criterion}{a function measuring the ``goodness of fit''}
  \item{startPar}{in case \code{optim} is used: a starting value for the parameter fit;
  in case \code{optimize} is used: a vector containing a search interval for the
  (one-dim) parameter}
  \item{penalty}{numeric; penalizes non-permitted parameter values}
  \item{crit.name}{character; the name of the criterion; may be missing}
  \item{withthetaPar}{logical; shall Parameter theta be transmitted?}
  \item{Infos}{matrix; info slot to be filled in object of class \code{MCEstimate}; 
               may be missing}
  \item{\dots}{additional argument(s) for \code{optim} / \code{optimize}}
}
\value{
  a list with components
  \item{\code{estimate}}{ --- the estimate as a named vector of \code{numeric}}
  \item{\code{criterion}}{ --- the criterion value (i.e.; a \code{numeric} of length 1); 
                           e.g. the neg. log likelihood}
  \item{\code{est.name}}{ --- the name of the estimator}
  \item{\code{param}}{ --- estimate coerced to class \code{ParamFamParameter}}
  \item{\code{crit.fct}}{ --- a function with the named components of \code{theta}
                         as arguments returning the criterion value; used for
                         profiling / coercing to class \code{mle}}
  \item{\code{method}}{ --- a character reporting how the estimate was obtained, i.e.,
                       by  \code{optim}, by \code{optimize} or by explicit calculations}
  \item{crit.name}{character; the name of the criterion; may be \code{""}}
  \item{Infos}{matrix; info slot to be filled in object of class \code{MCEstimate}; 
               may be \code{NULL}}
  \item{samplesize}{numeric; sample size of \code{x}}
  }

\details{
  \code{mceCalc} is used internally by function \code{MCEstimator} to allow
   for method dispatch according to argument \code{PFam};
   similarly, and for the same purpose  \code{mleCalc} is used internally by function 
   \code{MLEstimator}.
  This way we / or any other developper can write particular methods
  for special cases where we may avoid using numerical optimization
  without interfering with existing code. 
For programming one's own \code{mleCalc} / \code{mceCalc} methods, there is
the helper function \code{\link{meRes}} to produce consistent return values. 

  }

\keyword{models}
