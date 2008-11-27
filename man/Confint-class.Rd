\name{Confint-class}
\docType{class}
\alias{Confint-class}
\alias{type,Confint-method}
\alias{call.estimate}
\alias{call.estimate,Confint-method}
\alias{confint,Confint,missing-method}
\alias{name.estimate}
\alias{name.estimate,Confint-method}
\alias{trafo.estimate}
\alias{trafo.estimate,Confint-method}
\alias{samplesize.estimate}
\alias{samplesize.estimate,Confint-method}
\alias{nuisance.estimate}
\alias{nuisance.estimate,Confint-method}
\alias{fixed.estimate}
\alias{fixed.estimate,Confint-method}
\alias{show,Confint-method}
\alias{print,Confint-method}

\title{Confint-class}
\description{Return value S4 classes for method ``confint''.}
\section{Objects from the Class}{
  Objects could in principle be created by calls of the 
  form \code{new("Confint", ...)}.
  The preferred form is to have them created via a call to 
  \code{\link{confint}}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}:
      type of the confidence interval (asymptotic, bootstrap,\dots). 
      Can be of length \code{>2}. Then in printing, the first element
      is printed in the gap '[...]' in 'an [...] confidence interval',
      while the other elements are printed below.}
    \item{\code{confint}:}{Object of class \code{"array"}:
      the confidence interval(s).}
    \item{\code{call.estimate}:}{Object of class \code{"call"}:
      the estimate(s) for which the confidence intervals are produced.}
    \item{\code{name.estimate}:}{Object of class \code{"character"}:
      the name of the estimate(s) for which the confidence intervals are produced.}
    \item{\code{samplesize.estimate}:}{Object of class \code{"numeric"}:
      the sample size of the estimate(s) for which the confidence intervals are produced.}
    \item{\code{trafo.estimate}:}{Object of class \code{"matrix"}:
      the trafo/derivative matrix of the estimate(s) for which 
      the confidence intervals are produced.}
    \item{\code{nuisance.estimate}:}{Object of class \code{"OptionalNumeric"}:
      the nuisance parameter (if any) at which the confidence 
      intervals are produced.}
    \item{\code{fixed.estimate}:}{Object of class \code{"OptionalNumeric"}:
      the fixed part of the parameter (if any) at which the confidence 
      intervals are produced.}
  }
}
\section{Methods}{
  \describe{
    \item{type}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{type}. }

    \item{confint}{\code{signature(object = "Confint", method = "missing")}: 
      accessor function for slot \code{type}. }
    
    \item{call.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{call.estimate}. }

    \item{name.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{name.estimate}. }

    \item{trafo.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{trafo.estimate}. }

    \item{samplesize.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{samplesize.estimate}. }

    \item{nuisance.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{nuisance.estimate}. }

    \item{fixed.estimate}{\code{signature(object = "Confint")}: 
      accessor function for slot \code{fixed.estimate}. }

    \item{show}{\code{signature(object = "Confint")}: shows a detailed view
     of the object; slots \code{nuisance.estimate} and
     \code{fixed.estimate} are only shown if non-null,
     and slot \code{trafo.estimate} only if different from a unit matrix.}
    \item{print}{\code{signature(object = "Confint")}: just as \code{show}, 
     but with additional arguments \code{digits}.}
  }
}
\section{Details for methods \code{show}, \code{print}}{
Detailedness of output by methods \code{show}, \code{print} is controlled
by the global option \code{show.details} to be set by 
\code{\link{distrModoptions}}.

As method \code{show} is used when inspecting an object by typing the object's 
name into the console, \code{show} comes without extra arguments and hence 
detailedness must be controlled  by global options. 

Method \code{print} may be called with a (partially matched) argument
\code{show.details}, and then the global option is temporarily set to this
value.

More specifically, when \code{show.detail} is matched to \code{"minimal"}
you will be shown only the type of the confidence interval(s) and its/their 
values. When \code{show.detail} is matched to \code{"medium"}, you will in
addition see the type of the estimator(s) for which it is produced,
the corresponding call of the estimater, its sample size, and, if present, the 
value of the corresponding nuisance parameter. 
Finally, when \code{show.detail} is matched to \code{"maximal"}, additionally 
you will be shown the fixed part of the parameter (if present) and 
the transformation of the estimator (if non-trivial, i.e. the 
identity) in form of its function code respectively of its derivative matrix.
            }
\note{The pretty-printing code for methods \code{show} and \code{print} 
has been borrowed from \code{confint.default} in package \pkg{stats}.}

%\references{}
\author{Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{Estimator}}, \code{\link{confint}}, \code{\link{Estimate-class}},
\code{\link{trafo-methods}}}
\examples{
## some transformation
mtrafo <- function(x){
     nms0 <- c("scale","shape")
     nms <- c("shape","rate")
     fval0 <- c(x[2], 1/x[1])
     names(fval0) <- nms
     mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
                     dimnames = list(nms,nms0))                          
     list(fval = fval0, mat = mat0)}

x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2, trafo = mtrafo)
## MLE
res <- MLEstimator(x = x, ParamFamily = G)
ci <- confint(res)
print(ci, digits = 4, show.details="maximal")
print(ci, digits = 4, show.details="medium")
print(ci, digits = 4, show.details="minimal")
}
\concept{confidence interval}
\keyword{classes}
