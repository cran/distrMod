\name{NonSymmetric-class}
\docType{class}
\alias{NonSymmetric-class}

\title{Class for Non-symmetric Functions}
\description{Class for non-symmetric functions.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("NonSymmetric")}.
  More frequently they are created via the generating function 
  \code{NonSymmetric}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      contains \dQuote{non-symmetric function} }
    \item{\code{SymmCenter}:}{Object of class \code{"NULL"} }
  }
}
\section{Extends}{
Class \code{"FunctionSymmetry"}, directly.\cr
Class \code{"Symmetry"}, by class \code{"FunctionSymmetry"}.
}
%\section{Methods}{}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{NonSymmetric}}}
\examples{
new("NonSymmetric")
}
\concept{symmetry}
\keyword{classes}
