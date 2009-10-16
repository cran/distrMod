\name{meRes}
\alias{meRes}
\alias{get.criterion.fct}
\alias{samplesize,numeric-method}
\title{helper functions for mceCalc and mleCalc}

\description{helper functions to produce consistent lists to be digested 
in functions \code{\link{mceCalc}} and \code{\link{mleCalc}}}

\usage{
meRes(x, estimate, criterion.value, param, crit.fct, method = "explicit solution",
      crit.name = "Maximum Likelihood", Infos, warns = "", startPar = NULL)
get.criterion.fct(theta, Data, ParamFam, criterion.ff, fun, ...)
\S4method{samplesize}{numeric}(object)
}
\arguments{
  \item{x}{numeric; the data at which to evaluate the estimate}
  \item{estimate}{numeric; the estimate} 
  \item{criterion.value}{numeric; the value of the criterion}
  \item{param}{object of class \code{ParamFamParameter}; the parameter value} 
  \item{crit.fct}{a function to fill slot \code{minuslogl} when an object of 
      class \code{MCEstimate} is coerced to class \code{mle} 
      (from package \pkg{stats4}); to this end function \code{get.criterion.fct}
      (also see details below) is helpful (at least if the dimension of the 
      estimator is larger than 1).} 
  \item{method}{character; describes how the estimate was obtained} 
  \item{crit.name}{character; name of the criterion} 
  \item{Infos}{optional matrix of characters in two columns; 
               information to be attached to the estimate} 
  \item{warns}{collected warnings in optimization}
  \item{samplesize}{numeric; the sample size at which the estimator 
                    was evaluated}
  \item{theta}{the parameter value as named numeric vector}
  \item{Data}{numeric; the data at which to evaluate the MCE}
  \item{ParamFam}{an object of class \code{ParamFamily}; 
              the parametric family at which to evaluate the MCE}
  \item{criterion.ff}{the criterion function used in the MCE}
  \item{fun}{wrapper to the criterion function used in the MCE
   (with certain checking whether parameter value is permitted and possibly
    penalizing if not; see code to , for example.)}
  \item{startPar}{value of argument \code{StartPar} --- starting parameter used.}
  \item{\dots}{further arguments to be passed 
               to \code{optim} / \code{optimize}} 
  \item{object}{numeric; the data at which to evaluate the estimate}
  }

\details{ 
\code{get.criterion.fct} produces a function \code{criterion.fct} 
      to fill slot \code{minuslogl} when an object of class \code{MCEstimate}
      is coerced to class \code{mle} (from package \pkg{stats4});
      this way we may use profiling methods introduced there also for objects
      of our classes. More specifically, we produce a function where all 
      coordinates/components of \code{theta} appear as separate named 
      arguments, which then calls \code{fun} with these separate arguments 
      again stacked to one (named) vector argument;

\code{samplesize} determines the samplesize of argument \code{object},i.e.;
 if \code{object} has an attribute \code{dim}, it returns \code{dim(object)[2]},
 else \code{length(object)}.
}

\value{
\item{meRes}{a list of prescribed structure to be digested in functions 
\code{\link{mceCalc}} and \code{\link{mleCalc}} by the internal helper function 
\code{\link{.process.meCalcRes}}.}
\item{get.criterion.fct}{a function; see details below;}
\item{samplesize}{numeric}
}

\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\keyword{models}

