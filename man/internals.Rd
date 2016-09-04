\name{internals_for_distrMod}
\alias{internals_for_distrMod}
\alias{.inArgs}
\alias{.isUnitMatrix}
\alias{.validTrafo}
\alias{.CvMMDCovariance}
\alias{.show.with.sd}
\alias{.getLogDeriv}
\alias{.deleteDim}

\title{Internal functions of package distrMod}

\description{
These functions are used internally by package ``distrMod''.}

\usage{
.inArgs(arg, fct)
.isUnitMatrix(m)
.csimpsum(fx)
.validTrafo(trafo, dimension, dimensionwithN)
.CvMMDCovariance(L2Fam, param, mu = distribution(L2Fam),  
                 withplot = FALSE, withpreIC = FALSE,
                 N = getdistrOption("DefaultNrGridPoints")+1,
                 rel.tol=.Machine$double.eps^0.3, 
                 TruncQuantile = getdistrOption("TruncQuantile"), 
                 IQR.fac = 15, ...)
.show.with.sd(est, s)
.getLogDeriv(distr, 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
                         IQR.fac = getdistrExOption("IQR.fac"))
.deleteDim(x)
}


\arguments{
  \item{arg}{a formal argument as character}
  \item{fct}{a function}
  \item{m}{a matrix} 
  \item{est}{an estimator; usually a vector}
  \item{s}{a standard deviation} 
  \item{trafo}{an object of class \code{MatrixorFunction}}
  \item{dimension}{a numeric --- length of main part of the parameter}
  \item{dimensionwithN}{a numeric --- length of main and nuisance part of the parameter}

  \item{L2Fam}{an object of class \code{L2ParamFamily} --- for
  which we want to determine the IC resp. the as. [co]variance of the corresponding
  Minimum CvM estimator}
  \item{mu}{an object of class \code{UnivariateDistribution}: integration 
            measure (resp. distribution) for CvM distance}
  \item{withplot}{logical; defaults to \code{FALSE}; if \code{TRUE} for diagnostic  
                  purposes plots the influence function of the CvM-MDE}
  \item{withpreIC}{logical: shall we return a list with components \code{preIC} 
       and \code{var} or just \code{var}; here \code{var} is the corresponding
       asymptotic variance and \code{preIC} the corresponding
   \code{EuclRandVarList} featuring as argument \code{Curve} in \code{IC}s of
   package \pkg{RobAStBase}}
  \item{\dots}{currently not used} 
  \item{N}{integer; the number of grid points at which to evaluate the antiderivative
          in case of an absolutely continuous distribution; more precisely, internally
		  this becomes \eqn{2N+1}}
  \item{rel.tol}{relative tolerance for \code{distrExIntegrate}.} 
  \item{TruncQuantile}{numeric in (0,1); in case of an unbounded support of the 
                       distribution the quantile level at which to cut the distribution}
  \item{IQR.fac}{a positive numeric; a factor by which to multiply the IQR of the distribution
                 to obtain a sensible cut of point for the integration bounds}
  \item{lowerTruncQuantile}{lower quantile for quantile based integration range.}
  \item{upperTruncQuantile}{upper quantile for quantile based integration range.}
 \item{fx}{a vector of function evaluations multiplied by the gridwidth}
 \item{distr}{an object of class \code{AbscontDistribution}}
  \item{\dots}{further argument to be passed through --- so 
              \code{.CvMMDCovariance} can digest more arguments}
  \item{x}{a possibly named vector, which may have a \code{dim} attribute}
}

\details{
\code{.inArgs} (borrowed from package \pkg{distr})
      checks whether an argument \code{arg} is a formal argument of 
     \code{fct} --- not vectorized.

\code{.csimpsum} (borrowed from package \pkg{distr}) 
produces a primitimive function out of function evaluations by means 
of vectorized Simpson quadrature method, returning already the function values
of the prime function on a grid; it is to mimick the behaviour
of \code{cumsum}. 

\code{.isUnitMatrix} checks whether the argument is a unit matrix. 

\code{.validTrafo} checks whether the argument is a valid transformation.

\code{.CvMMDCovariance} computes the asymptotic covariance of the CvM-MDE 
according to H. Rieder (1994) "Robust Asymptotic Statistics".

\code{.show.with.sd} is code borrowed from \code{print.fitdistr} in  
package \pkg{MASS} by B.D. Ripley. It pretty-prints estimates with corresponding
sd's below.

\code{.getLogDeriv} determines numerically the negative logarithmic derivative of the 
density of distribution \code{distr}; to this end uses \code{D1ss}, 
\code{D2ss} from Martin Maechler's package \pkg{sfsmisc}. 

\code{.deleteDim} deletes a possible \code{dim} argument (sets it to \code{NULL})
but retains all other possible attributes, in particular a \code{name} attribute.
}

\value{
\item{.getLogderiv}{a function in one argument \code{x} --- the negative logarithmic
derivative of the density}
\item{.inArgs}{\code{logical} (length 1)}
\item{.csimpsum}{\code{numeric} (of length half the input length)}
\item{.isUnitMatrix}{\code{logical} (length 1)}
\item{.validTrafo}{\code{logical} (length 1)}
\item{.CvMMDCovariance}{corresponding as. [co]variance of
   the corresponding Minimum CvM estimator or list withcomponents
   \code{preIC} and \code{var} ---see above}
\item{.show.with.sd}{\code{invisible()}}
\item{.deleteDim}{vector \code{x} without \code{dim} attribute}
\item{.CvMMDCovariance}{if argument \code{withpreIC} is \code{TRUE}, then
a list with elements the IC of the CvM-MDE and its covariance is returned, otherwise
just the covariance}
}

\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}
  Matthias Kohl \email{Matthias.Kohl@stamats.de}
}

\seealso{
\code{\link{MLEstimator}},
\code{\link{Estimate-class}},
\code{\link{MCEstimate-class}},
\code{\link{Confint-class}},
\code{\link{ParamFamParameter-class}}
}
\keyword{internal}
\concept{utilities}
