\name{OddSymmetric}
\alias{OddSymmetric}

\title{Generating function for OddSymmetric-class}
\description{Generates an object of class \code{"OddSymmetric"}.}
\usage{OddSymmetric(SymmCenter = 0)}
\arguments{
  \item{SymmCenter}{ numeric: center of symmetry }
}
%\details{}
\value{Object of class \code{"OddSymmetric"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de}}
%\note{}
\seealso{\code{\link{OddSymmetric-class}}, \code{\link{FunctionSymmetry-class}}}
\examples{
OddSymmetric()

## The function is currently defined as
function(SymmCenter = 0){ 
    new("OddSymmetric", SymmCenter = SymmCenter) 
}
}
\concept{odd function}
\concept{symmetry}
\keyword{math}
