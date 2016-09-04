\name{ParamFamily}
\alias{ParamFamily}

\title{Generating function for ParamFamily-class}
\description{
  Generates an object of class \code{"ParamFamily"}.
}
\usage{
ParamFamily(name, distribution = Norm(), distrSymm, modifyParam,
            main = main(param), nuisance = nuisance(param),
            fixed = fixed(param), trafo = trafo(param),
            param = ParamFamParameter(name = paste("Parameter of", 
                          name),  main = main, nuisance = nuisance, 
                                  fixed = fixed, trafo = trafo),
            props = character(0),
            startPar = NULL, makeOKPar = NULL)
}
\arguments{
  \item{name}{ character string: name of family }
  \item{distribution}{ object of class \code{"Distribution"}: 
    member of the family }
  \item{distrSymm}{ object of class \code{"DistributionSymmetry"}: 
      symmetry of \code{distribution}. }
  \item{startPar}{\code{startPar} is a function in the observations \code{x}
    returning initial information for \code{MCEstimator} used 
    by \code{optimize} resp. \code{optim}; i.e; if (total) parameter is of 
    length 1, \code{startPar} returns a search interval, else it returns an initial 
    parameter value.}
  \item{makeOKPar}{\code{makeOKPar} is a function in the (total)
    parameter \code{param}; used if \code{optim} resp. \code{optimize}---
           try to use ``illegal'' parameter values; then \code{makeOKPar} makes
           a valid parameter value out of the illegal one; if \code{NULL}
    slot \code{makeOKPar} of \code{ParamFamily} is used to produce it.}
  \item{main}{ numeric vector: main parameter }
  \item{nuisance}{ numeric vector: nuisance parameter }
  \item{fixed}{ numeric vector: fixed part of the parameter }
  \item{trafo}{ function in \code{param} or matrix: transformation of the parameter }
  \item{param}{ object of class \code{"ParamFamParameter"}:
    parameter of the family }
  \item{modifyParam}{ function: mapping from the parameter space 
  (represented by \code{"param"}) to the distribution space 
  (represented by \code{"distribution"}). }
  \item{props}{ character vector: properties of the family }
}
\details{
  If \code{name} is missing, the default 
  \dQuote{"parametric family of probability measures"} is used. 
  In case \code{distrSymm} is missing it is set to \code{NoSymmetry()}.
  If \code{param} is missing, the parameter is created via
  \code{main}, \code{nuisance} and \code{trafo} as described
  in \code{\link{ParamFamParameter}}.
  One has to specify a function which represents a mapping
  from the parameter space to the corresponding distribution
  space; e.g., in case of normal location a simple version of such
  a function would be \code{function(theta){ Norm(mean = theta) }}.
}
\value{Object of class \code{"ParamFamily"}}
%\references{}
\author{Matthias Kohl \email{Matthias.Kohl@stamats.de},\cr
Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de}}
%\note{}
\seealso{\code{\link{ParamFamily-class}}}
\examples{

## "default" (normal location)
F1 <- ParamFamily(modifyParam = function(theta){ Norm(mean = theta) })
plot(F1)

################################
## Some examples:
################################
## 1. Normal location family
theta <- 0
names(theta) <- "mean"
NL <- ParamFamily(name = "Normal location family",
          param = ParamFamParameter(name = "location parameter", main = theta),
          distribution = Norm(mean = 0, sd = 1), ## sd known!
          startPar = function(x,...) c(min(x),max(x)),
          distrSymm <- SphericalSymmetry(SymmCenter = 0),
          modifyParam = function(theta){ Norm(mean = theta, sd = 1) },
          props = paste(c("The normal location family is invariant under",
                    "the group of transformations 'g(x) = x + mean'",
                    "with location parameter 'mean'"), collapse = " "))
NL

## 2. Normal scale family
theta <- 1
names(theta) <- "sd"
NS <- ParamFamily(name = "Normal scale family",
          param = ParamFamParameter(name = "scale parameter", main = theta,
          .returnClsName = "ParamWithScaleFamParameter"),
          distribution = Norm(mean = 0, sd = 1), ## mean known!
          startPar = function(x,...) c(0,-min(x)+max(x)),
          distrSymm <- SphericalSymmetry(SymmCenter = 0),
          modifyParam = function(theta){ Norm(mean = 0, sd = theta) },
          props = paste(c("The normal scale family is invariant under",
                    "the group of transformations 'g(y) = sd*y'",
                    "with scale parameter 'sd'"), collapse = " "))
NS

## 3. Normal location and scale family
theta <- c(0, 1)
names(theta) <- c("mean", "sd")
NLS <- ParamFamily(name = "Normal location and scale family",
          param = ParamFamParameter(name = "location and scale parameter",
                                    main = theta,
                                 .returnClsName = "ParamWithScaleFamParameter"),
          distribution = Norm(mean = 0, sd = 1),
          startPar = function(x,...) c(median(x),mad(x)),
          makeOKPar = function(param) {param[2]<-abs(param[2]); return(param)},
          distrSymm <- SphericalSymmetry(SymmCenter = 0),
          modifyParam = function(theta){
                            Norm(mean = theta[1], sd = theta[2])
                        },
          props = paste(c("The normal location and scale family is",
                    "invariant under the group of transformations",
                    "'g(x) = sd*x + mean' with location parameter",
                    "'mean' and scale parameter 'sd'"),
                    collapse = " "))
NLS

## 4. Binomial family
theta <- 0.3
names(theta) <- "prob"
B <- ParamFamily(name = "Binomial family",
         param = ParamFamParameter(name = "probability of success", 
                                   main = theta),
         startPar = function(x,...) c(0,1),
         distribution = Binom(size = 15, prob = 0.3), ## size known!
         modifyParam = function(theta){ Binom(size = 15, prob = theta) },
         props = paste(c("The Binomial family is symmetric with respect",
                   "to prob = 0.5; i.e.,",
                   "d(Binom(size, prob))(k)=d(Binom(size,1-prob))(size-k)"),
                   collapse = " "))
B

## 5. Poisson family
theta <- 7
names(theta) <- "lambda"
P <- ParamFamily(name = "Poisson family",
          param = ParamFamParameter(name = "positive mean", main = theta),
          startPar = function(x,...) c(0,max(x)),
          distribution = Pois(lambda = 7),
          modifyParam = function(theta){ Pois(lambda = theta) })
P


## 6. Exponential scale family
theta <- 2
names(theta) <- "scale"
ES <- ParamFamily(name = "Exponential scale family",
          param = ParamFamParameter(name = "scale parameter", main = theta,
                           .returnClsName = "ParamWithScaleFamParameter"),
          startPar = function(x,...) c(0,max(x)-min(x)),
          distribution = Exp(rate = 1/2),
          modifyParam = function(theta){ Exp(rate = 1/theta) },
          props = paste(c("The Exponential scale family is invariant under",
                    "the group of transformations 'g(y) = scale*y'",
                    "with scale parameter 'scale = 1/rate'"),
                    collapse = " " ))
ES

## 7. Lognormal scale family
theta <- 2
names(theta) <- "scale"
LS <- ParamFamily(name = "Lognormal scale family",
          param = ParamFamParameter(name = "scale parameter", main = theta,
                           .returnClsName = "ParamWithScaleFamParameter"),
          startPar = function(x,...) c(0,max(x)-min(x)),
          distribution = Lnorm(meanlog = log(2), sdlog = 2),## sdlog known!
          modifyParam = function(theta){ 
                            Lnorm(meanlog = log(theta), sdlog = 2) 
                        },
          props = paste(c("The Lognormal scale family is invariant under",
                    "the group of transformations 'g(y) = scale*y'",
                    "with scale parameter 'scale = exp(meanlog)'"),
                    collapse = " "))
LS

## 8. Gamma family
theta <- c(1, 2)
names(theta) <- c("scale", "shape")
G <- ParamFamily(name = "Gamma family",
        param = ParamFamParameter(name = "scale and shape", main = theta,
                           withPosRestr = TRUE,
                           .returnClsName = "ParamWithScaleAndShapeFamParameter"),
        startPar = function(x,...) {E <- mean(x); V <- var(X); c(V/E,E^2/V)},
        makeOKPar = function(param) abs(param),
        distribution = Gammad(scale = 1, shape = 2),
        modifyParam = function(theta){ 
                          Gammad(scale = theta[1], shape = theta[2]) 
                      },
        props = paste(c("The Gamma family is scale invariant via the",
                  "parametrization '(nu,shape)=(log(scale),shape)'"),
                  collapse = " "))
G
}
\concept{parametric family}
\keyword{distribution}
\keyword{models}
