\name{MCEstimate-class}
\docType{class}
\alias{MCEstimate-class}
\alias{criterion}
\alias{criterion,MCEstimate-method}
\alias{criterion.fct}
\alias{criterion.fct,MCEstimate-method}
\alias{startPar,MCEstimate-method}
\alias{method}
\alias{method,MCEstimate-method}
\alias{optimwarn}
\alias{optimwarn,MCEstimate-method}
\alias{criterion<-}
\alias{criterion<-,MCEstimate-method}
\alias{coerce,MCEstimate,mle-method}
\alias{show,MCEstimate-method}
\alias{profile,MCEstimate-method}

\title{MCEstimate-class.}
\description{Class of minimum criterion estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MCEstimate", ...)}.
  More frequently they are created via the generating functions
  \code{MCEstimator}, \code{MDEstimator} or \code{MLEstimator}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate.}
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced.}
    \item{\code{criterion}}{Object of class \code{"numeric"}:
      minimum value of the considered criterion.}
    \item{\code{criterion.fct}}{Object of class \code{"function"}:
      the considered criterion function; used for compatibility with class
      \code{"mle"} from package \pkg{stats4}; should be a function
      returning the criterion; i.e. a numeric of length 1 and should have
      as arguments all named components of argument 
      \code{untransformed.estimate}}
    \item{\code{method}}{Object of class \code{"character"}:
      the method by which the estimate was calculated, i.e.;
      \code{"optim"}, \code{"optimize"}, or \code{"explicit calculation"}; 
      used for compatibility with class \code{"mle"} from package 
      \pkg{stats4}, could be any character value.}
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
    \item{\code{optimwarn}}{ object of class \code{"character"}
      warnings issued during optimization. }
    \item{\code{startPar}}{ --- object of class \code{"ANY"}; filled either
     with \code{NULL} (no starting value used) or with \code{"numeric"}
      --- the value of the starting parameter. }
    \item{\code{asvar}}{ object of class \code{"OptionalMatrix"}
      which may contain the asymptotic (co)variance of the estimator. }
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize at which the estimate was evaluated. }
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}: 
     indices of \code{estimate} belonging to the nuisance part}
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}:
     the fixed and known part of the parameter. }
    \item{\code{trafo}}{ object of class \code{"list"}:
     a list with components \code{fct} and \code{mat} (see below). }
    \item{\code{untransformed.estimate}}{Object of class \code{"ANY"}:
      untransformed estimate.}
    \item{\code{untransformed.asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed 
      estimator. }
    \item{\code{completecases}}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. }
    \item{\code{startPar}}{ object of class \code{"ANY"}; usually filled with
    argument \code{startPar} of generating function \code{MCEstimator},
    \code{MLEstimator}, \code{MDEstimator}.}
  }
}
\section{Extends}{
Class \code{"Estimate"}, directly.
}
\section{Methods}{
  \describe{
    \item{criterion}{\code{signature(object = "MCEstimate")}: 
      accessor function for slot \code{criterion}. }

    \item{criterion<-}{\code{signature(object = "MCEstimate")}: 
      replacement function for slot \code{criterion}. }

    \item{optimwarn}{\code{signature(object = "MCEstimate")}:
      accessor function for slot \code{optimwarn}. }

    \item{startPar}{\code{signature(object = "MCEstimate")}:
      accessor function for slot \code{startPar}. }

    \item{criterion.fct}{\code{signature(object = "MCEstimate")}:
      accessor function for slot \code{criterion.fct}. }

    \item{show}{\code{signature(object = "Estimate")}}

    \item{coerce}{\code{signature(from = "MCEstimate", to = "mle")}:
      create a \code{"mle"} object from a \code{"MCEstimate"} object}

    \item{profile}{\code{signature(fitted = "MCEstimate")}:
      coerces \code{fitted} to class \code{"mle"} and then calls
      the corresponding \code{\link[stats4:profile-methods]{profile}}-method
      from package \pkg{stats4}; for details we confer to the corresponding
      man page.}
  }
}

%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
%\note{}
\seealso{\code{\link{Estimate-class}}, \code{\link{MCEstimator}}, 
         \code{\link{MDEstimator}}, \code{\link{MLEstimator}}}
\examples{
## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

MDEstimator(x, G)
(m <- MLEstimator(x, G))
m.mle <- as(m,"mle")
par(mfrow=c(1,2))
plot(profile(m))
}
\concept{estimate}
\keyword{classes}
