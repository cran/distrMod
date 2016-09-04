\name{Estimator}
\alias{Estimator}

\title{ Function to compute estimates }
\description{
  The function \code{Estimator} provides a general way to compute 
  estimates.
}
\usage{
Estimator(x, estimator, name, Infos, asvar = NULL, nuis.idx,
          trafo = NULL, fixed = NULL, asvar.fct, na.rm = TRUE, ...,
          ParamFamily = NULL, .withEvalAsVar = TRUE)
}
\arguments{
  \item{x}{ (empirical) data }
  \item{estimator}{ function: estimator to be evaluated on \code{x}. }
  \item{name}{ optional name for estimator. }
  \item{Infos}{ character: optional informations about estimator }
  \item{asvar}{ optionally the asymptotic (co)variance of the estimator }
  \item{nuis.idx}{ optionally the indices of the estimate belonging 
                  to nuisance parameter}
  \item{fixed}{ optionally (numeric) the fixed part of the parameter}
  \item{trafo}{ an object of class \code{MatrixorFunction} -- a transformation
  for the main parameter}
  \item{asvar.fct}{optionally: a function to determine the corresponding
    asymptotic variance; if given, \code{asvar.fct} takes arguments
    \code{L2Fam}(the parametric model as object of class \code{L2ParamFamily})
    and \code{param} (the parameter value as object of class 
    \code{ParamFamParameter}); arguments are called by name; \code{asvar.fct}
     may also process further arguments passed through the \code{\dots} argument. }
  \item{na.rm}{logical: if  \code{TRUE}, the estimator is evaluated at \code{complete.cases(x)}.}
  \item{\dots}{ further arguments to \code{estimator}.}
  \item{ParamFamily}{an optional object of class \code{ParamFamily}. Passed on
     to \code{asvar.fct} to compute asymptotic variances.}
  \item{.withEvalAsVar}{logical: shall slot \code{asVar} be evaluated
                   (if \code{asvar.fct} is given) or
                   just the call be returned?}

}
\details{
  The argument \code{criterion} has to be a function with arguments the 
  empirical data as well as an object of class \code{"Distribution"} 
  and possibly \code{\dots}.
}
\value{
  An object of S4-class \code{"Estimate"}.
}
%\references{  }
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{Estimate-class}} }
\examples{
x <- rnorm(100)
Estimator(x, estimator = mean, name = "mean")

X <- matrix(rnorm(1000), nrow = 10)
Estimator(X, estimator = rowMeans, name = "mean")
}
\keyword{univar}
