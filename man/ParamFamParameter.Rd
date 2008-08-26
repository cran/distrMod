\name{ParamFamParameter}
\alias{ParamFamParameter}

\title{Generating function for ParamFamParameter-class}
\description{
  Generates an object of class \code{"ParamFamParameter"}.
}
\usage{
ParamFamParameter(name, main = numeric(0), nuisance, fixed, trafo)
}
\arguments{
  \item{name}{ (optional) character string: name of parameter }
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ (optional) numeric vector: nuisance paramter }
  \item{fixed}{ (optional) numeric vector: fixed part of the paramter }
  \item{trafo}{(optional) MatrixorFunction: transformation of the parameter }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{"parameter of a parametric family of probability measures"}
  is used. If \code{nuisance} is missing, the nuisance parameter is 
  set to \code{NULL}. The number of columns of \code{trafo} have
  to be equal and the number of rows have to be not larger than 
  the sum of the lengths of \code{main} and \code{nuisance}. 
  If \code{trafo} is missing, no transformation to the parameter
  is applied; i.e., \code{trafo} is set to an identity matrix.
}
\value{Object of class \code{"ParamFamParameter"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{ParamFamParameter-class}}}
\examples{
ParamFamParameter(main = 0, nuisance = 1, fixed = 2,
                  trafo = function(x) list(fval = sin(x), 
                                            mat = matrix(cos(x),1,1))
                  )                          

}

\concept{parameter}
\keyword{models}
