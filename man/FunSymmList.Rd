\name{FunSymmList}
\alias{FunSymmList}

\title{Generating function for FunSymmList-class}
\description{
  Generates an object of class \code{"FunSymmList"}.
}
\usage{
FunSymmList(...)
}
\arguments{
  \item{\dots}{ Objects of class \code{"FunctionSymmetry"} which
    shall form the list of symmetry types. }
}
%\details{}
\value{ Object of class \code{"FunSymmList"} }
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{FunSymmList-class}}}
\examples{
FunSymmList(NonSymmetric(), EvenSymmetric(SymmCenter = 1), 
            OddSymmetric(SymmCenter = 2))

## The function is currently defined as
function (...){
    new("FunSymmList", list(...))
}
}
\concept{symmetry}
\keyword{math}
