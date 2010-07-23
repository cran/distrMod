\name{Estimate-class}
\docType{class}
\alias{Estimate-class}
\alias{name,Estimate-method}
\alias{name<-,Estimate-method}
\alias{estimate}
\alias{estimate,Estimate-method}
\alias{estimate.call}
\alias{estimate.call,Estimate-method}
\alias{Infos}
\alias{Infos,Estimate-method}
\alias{samplesize}
\alias{samplesize,Estimate-method}
\alias{completecases}
\alias{completecases,Estimate-method}
\alias{asvar}
\alias{asvar,Estimate-method}
\alias{fixed,Estimate-method}
\alias{asvar<-}
\alias{asvar<-,Estimate-method}
\alias{nuisance,Estimate-method}
\alias{main,Estimate-method}
\alias{Infos<-}
\alias{Infos<-,Estimate-method}
\alias{addInfo<-}
\alias{addInfo<-,Estimate-method}
\alias{show,Estimate-method}
\alias{print,Estimate-method}
\alias{untransformed.estimate}
\alias{untransformed.estimate,Estimate-method}
\alias{untransformed.asvar}
\alias{untransformed.asvar,Estimate-method}

\title{Estimate-class.}
\description{Class of estimates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Estimate", ...)}.
  More frequently they are created via the generating function 
  \code{Estimator}. 
}
\section{Slots}{
  \describe{
    \item{\code{name}}{Object of class \code{"character"}:
      name of the estimator. }
    \item{\code{estimate}}{Object of class \code{"ANY"}:
      estimate.}
    \item{\code{estimate.call}}{Object of class \code{"call"}:
      call by which estimate was produced.}
    \item{\code{Infos}}{ object of class \code{"matrix"}
      with two columns named \code{method} and \code{message}:
      additional informations. }
    \item{\code{asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the estimator. }
    \item{\code{samplesize}}{ object of class \code{"numeric"} ---
      the samplesize (only complete cases are counted)
      at which the estimate was evaluated. }
    \item{\code{completecases}}{ object of class \code{"logical"} ---
         complete cases at which the estimate was evaluated. }
    \item{\code{nuis.idx}}{ object of class \code{"OptionalNumeric"}:
     indices of \code{estimate} belonging to the nuisance part. }
    \item{\code{fixed}}{ object of class \code{"OptionalNumeric"}: 
     the fixed and known part of the parameter. }
    \item{\code{trafo}}{ object of class \code{"list"}: 
     a list with components \code{fct} and \code{mat} (see below). }
    \item{\code{untransformed.estimate}}{Object of class \code{"ANY"}:
      untransformed estimate. }
    \item{\code{untransformed.asvar}}{ object of class \code{"OptionalNumericOrMatrix"}
      which may contain the asymptotic (co)variance of the untransformed 
      estimator. }
  }
}
\section{Methods}{
  \describe{
    \item{name}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{name}. }

    \item{name<-}{\code{signature(object = "Estimate")}: 
      replacement function for slot \code{name}. }

    \item{estimate}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{estimate}. }

    \item{untransformed.estimate}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{untransformed.estimate}. }

    \item{estimate.call}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{estimate.call}. }

    \item{samplesize}{\code{signature(object = "Estimate")}:
      (with additional argument \code{onlycompletecases}
      defaulting to \code{TRUE} returns the sample size;
      in case there are any incomplete cases and argument
      \code{onlycompletecases} is \code{FALSE}, the number of
      these is added to slot \code{samplesize}. }

    \item{completecases}{\code{signature(object = "Estimate")}:
      accessor function for slot \code{completecases}. }

    \item{asvar}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{asvar}. }

    \item{asvar<-}{\code{signature(object = "Estimate")}: 
      replacement function for slot \code{asvar}. }

    \item{untransformed.asvar}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{untransformed.asvar}. }

    \item{nuisance}{\code{signature(object = "Estimate")}: 
      accessor function for \code{nuisance} part of slot \code{estimate}. }
      
    \item{main}{\code{signature(object = "Estimate")}: 
      accessor function for \code{main} part of slot \code{estimate}. }

    \item{fixed}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{fixed}. }

    \item{Infos}{\code{signature(object = "Estimate")}: 
      accessor function for slot \code{Infos}. }

    \item{Infos<-}{\code{signature(object = "Estimate")}: 
      replacement function for slot \code{Infos}. }

    \item{addInfo<-}{\code{signature(object = "Estimate")}: 
      function to add an information to slot \code{Infos}. }

    \item{show}{\code{signature(object = "Estimate")}}

    \item{print}{\code{signature(object = "Estimate")}: just as \code{show}, 
     but with additional arguments \code{digits}.}
  }
}
%\references{}
\section{Details for methods 'show', 'print'}{
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
  you will be shown only the name/type of the estimator, the value of its main
  part, and, if present, the corresponding standard errors, as well as, 
  also if present, the value of the nuisance part. 
  When \code{show.detail} is matched to \code{"medium"}, you will in
  addition see the class of the estimator, its call and its sample-size
  and, if present, the fixed part of the parameter and 
  the asymptotic covariance matrix. Also the information
  gathered in the \code{Infos} slot is shown.
  Finally, when \code{show.detail} is matched to \code{"maximal"}, 
  and if, in addition, you estimate non-trivial (i.e. not the identity)
  transformation of the parameter of the parametric family, you will also be 
  shown this transformation in form of its function and its derivative matrix
  at the estimated parameter value, as well as the estimator (with standard errors,
  if present) and (again, if present) the corresponding asymptotic covariance
  of the untransformed, total (i.e. main and nuisance part) parameter.

  \code{trafo} realizes partial influence curves; i.e.; we are only
  interested is some possibly lower dimensional smooth (not necessarily
  linear or even coordinate-wise) aspect/transformation \eqn{\tau}{tau}
  of the parameter \eqn{\theta}{theta}.

  To be coherent with the corresponding \emph{nuisance}
  implementation, we make the following convention:

  The full parameter \eqn{\theta}{theta} is split up coordinate-wise
  in a main parameter \eqn{\theta'}{theta'} and a nuisance parameter
  \eqn{\theta''}{theta''} (which is unknown, too, hence has to be
  estimated, but only is of secondary interest) and a fixed,
  known part \eqn{\theta'''}{theta'''}.
  
  Without loss of generality, we restrict ourselves to the case that
  transformation \eqn{\tau}{tau} only acts on the main parameter
  \eqn{\theta'}{theta'} --- if we want to transform the whole
  parameter, we only have to assume that both nuisance parameter
  \eqn{\theta''}{theta''} and fixed, known part of the parameter
  \eqn{\theta'''}{theta'''} have length 0.
  
  To the implementation:
  
  Slot \code{trafo} can either contain a (constant) matrix
  \eqn{D_\theta}{D_theta} or a function
  \deqn{\tau\colon \Theta' \to \tilde \Theta,\qquad \theta \mapsto \tau(\theta)}{tau:  Theta' -> TTheta,  theta |-> tau(theta)} 
  mapping main parameter
  \eqn{\theta'}{theta'} to some range \eqn{\tilde \Theta}{TTheta}.

  If \emph{slot value} \code{trafo} is a function, besides \eqn{\tau(\theta)}{tau(theta)}, 
  it will also return the corresponding derivative matrix
  \eqn{\frac{\partial}{\partial \theta}\tau(\theta)}{(d/d theta) (tau(theta))}.
  More specifically, the return value of this function \code{theta} is a
  list with entries \code{fval}, the function value \eqn{\tau(\theta)}{tau(theta)},
  and \code{mat}, the derivative matrix.

  In case \code{trafo} is a matrix \eqn{D}, we interpret it as such a derivative
  matrix \eqn{\frac{\partial}{\partial \theta}\tau(\theta)}{(d/d theta) (tau(theta))},
  and, correspondingly, \eqn{\tau(\theta)}{tau(theta)} as the linear mapping
  \eqn{\tau(\theta)=D\,\theta}{tau(theta)=D * theta}.\cr
}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{Peter.Ruckdeschel@itwm.fraunhofer.de}}
\note{
The pretty-printing code for methods \code{show} and \code{print} 
has been borrowed from \code{print.fitdistr} in package \pkg{MASS} by B.D. Ripley.
}
\seealso{\code{\link{Estimator}}}
\examples{
x <- rnorm(100)
Estimator(x, estimator = mean, name = "mean")

x1 <- x; x1[sample(1:100,10)] <- NA
myEst1 <- Estimator(x1, estimator = mean, name = "mean")
samplesize(myEst1)
samplesize(myEst1, onlycomplete = FALSE)
}
\concept{estimate}
\keyword{classes}
