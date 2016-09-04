\name{internal_mlehelpers_for_distrMod}
\alias{internal_mlehelpers_for_distrMod}
\alias{.negLoglikelihood}
\alias{.get.criterion.fct}
\alias{.process.meCalcRes}
\alias{.samplesize}
\alias{.callParamFamParameter}

\title{Internal helper functions for treating MCEstimators in package distrMod}

\description{
These functions are used internally by functions \code{MCEstimator} and
\code{MLEstimator} in package ``distrMod''.}

\usage{
.negLoglikelihood(x, Distribution, ..., dropZeroDensity = TRUE)
.process.meCalcRes(res, PFam, trafo, res.name, call, asvar.fct, check.validity,
                   ..., .withEvalAsVar = TRUE)
.callParamFamParameter(PFam, theta, idx, nuis, fixed)
}


\arguments{
  \item{x}{numeric; the data at which to evaluate the negative log likelihood}
  \item{Distribution}{an object of class \code{Distribution}; 
                     the distribution at which to evaluate the negative log likelihood}
  \item{\dots}{in case of \code{.negLoglikelihood} simply to be able to use
               the general calling for the criterion function;
               in case of \code{.process.meCalcRes} further arguments to be passed 
               to function \code{asvar.fct}} 
  \item{res}{a list; the return value of \code{mleCalc} or \code{mceCalc}}
  \item{PFam}{an object of class \code{ParamFamily}; 
              the parametric family at which to evaluate the MCE}
  \item{trafo}{the trafo matrix or a function returning it or \code{NULL}} 
  \item{res.name}{the name of the corresponding estimator}
  \item{call}{the original call with which the estimate was produced}
  \item{asvar.fct}{a function producing an asymptotic covariance matrix or \code{NULL}.}
  \item{validity.check}{logical: shall return parameter value be checked for
    validity?}
  \item{theta}{parameter vector}
  \item{idx}{indices of the main part of the component}
  \item{nuis}{nuisance part of the parameter}
  \item{fixed}{fixed part of the parameter}
  \item{.withEvalAsVar}{logical: shall slot \code{asVar} be evaluated
                   (if \code{asvar.fct} is given) or
                   just the call be returned?}
  \item{dropZeroDensity}{logical of length 1; shall observations with
     density zero be dropped? Optimizers like \code{optim} require finite
     values, so get problems when negative loglikelihood is evaluated. }
}

\details{
\code{.negLoglikelihood} uses the \code{log} -argument of the corresponding \code{d}-slot
 of the distribution if available; else produces \code{log(d(Distribution)(x))}.
 
\code{.get.criterion.fct} produces a function \code{criterion.fct} 
      to fill slot \code{minuslogl} when an object of class \code{MCEstimate}
      is coerced to class \code{mle} (from package \pkg{stats4});
      this way we may use profiling methods introduced there also for objects
      of our classes. More specifically, we produce a function where all 
      coordinates/components of \code{theta} appear as separate named 
      arguments, which then calls \code{fun} with these separate arguments 
      again stacked to one (named) vector argument;

\code{.process.meCalcRes} processes the resulting return value list of methods
\code{mceCalc} and \code{mleCalc} to give a corresponding object of
class \code{MCEstimate}. 

\code{.callParamFamParameter} determines by means of the family-slot
parameter whether this is of a subclass of \code{ParamFamParameter},
and if so manipulates the call to generating function \code{ParamFamParameter}
accordingly (such that the result has the convenient type and the
convenient extra slots).
}

\value{
\item{negLoglikelihood}{the (numeric) value of the negative log likelihood
 evaluated at the data and the given distribution}
\item{.get.criterion.fct}{a function as described in section ``details''}
\item{.process.meCalcRes}{an object of class \code{MCEstimate}}
\item{.callParamFamParameter}{an object of class \code{ParamFamParameter}
or if necessary a corresponding subclass.}
}

\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
}

\seealso{
\code{\link{MCEstimate-class}},
\code{\link[stats4]{mle-class}},
}
\keyword{internal}
\concept{utilities}
