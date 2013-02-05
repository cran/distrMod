pkgname <- "distrMod"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('distrMod')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("BetaFamily")
### * BetaFamily

flush(stderr()); flush(stdout())

### Name: BetaFamily
### Title: Generating function for Beta families
### Aliases: BetaFamily
### Keywords: models

### ** Examples

(B1 <- BetaFamily())
FisherInfo(B1)
checkL2deriv(B1)



cleanEx()
nameEx("BiasType-class")
### * BiasType-class

flush(stderr()); flush(stdout())

### Name: BiasType-class
### Title: Bias Type
### Aliases: BiasType-class name,BiasType-method name<-,BiasType-method
### Keywords: classes

### ** Examples

aB <- positiveBias()
name(aB)



cleanEx()
nameEx("BinomFamily")
### * BinomFamily

flush(stderr()); flush(stdout())

### Name: BinomFamily
### Title: Generating function for Binomial families
### Aliases: BinomFamily
### Keywords: models

### ** Examples

(B1 <- BinomFamily(size = 25, prob = 0.25))
plot(B1)
FisherInfo(B1)
checkL2deriv(B1)



cleanEx()
nameEx("CauchyLocationScaleFamily")
### * CauchyLocationScaleFamily

flush(stderr()); flush(stdout())

### Name: CauchyLocationScaleFamily
### Title: Generating function for Cauchy location and scale families
### Aliases: CauchyLocationScaleFamily
### Keywords: models

### ** Examples

(C1 <- CauchyLocationScaleFamily())
plot(C1)
FisherInfo(C1)
### need smaller integration range:
distrExoptions("ElowerTruncQuantile"=1e-4,"EupperTruncQuantile"=1e-4)
checkL2deriv(C1)
distrExoptions("ElowerTruncQuantile"=1e-7,"EupperTruncQuantile"=1e-7)



cleanEx()
nameEx("Confint-class")
### * Confint-class

flush(stderr()); flush(stdout())

### Name: Confint-class
### Title: Confint-class
### Aliases: Confint-class type,Confint-method call.estimate
###   call.estimate,Confint-method confint,Confint,missing-method
###   name.estimate name.estimate,Confint-method trafo.estimate
###   trafo.estimate,Confint-method samplesize.estimate
###   samplesize.estimate,Confint-method completecases.estimate
###   completecases.estimate,Confint-method nuisance.estimate
###   nuisance.estimate,Confint-method fixed.estimate
###   fixed.estimate,Confint-method show,Confint-method
###   print,Confint-method
### Keywords: classes

### ** Examples

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



cleanEx()
nameEx("Estimate-class")
### * Estimate-class

flush(stderr()); flush(stdout())

### Name: Estimate-class
### Title: Estimate-class.
### Aliases: Estimate-class name,Estimate-method name<-,Estimate-method
###   estimate estimate,Estimate-method estimate.call
###   estimate.call,Estimate-method Infos Infos,Estimate-method samplesize
###   samplesize,Estimate-method completecases
###   completecases,Estimate-method asvar asvar,Estimate-method
###   fixed,Estimate-method asvar<- asvar<-,Estimate-method
###   nuisance,Estimate-method main,Estimate-method Infos<-
###   Infos<-,Estimate-method addInfo<- addInfo<-,Estimate-method
###   show,Estimate-method print,Estimate-method untransformed.estimate
###   untransformed.estimate,Estimate-method untransformed.asvar
###   untransformed.asvar,Estimate-method
### Keywords: classes

### ** Examples

x <- rnorm(100)
Estimator(x, estimator = mean, name = "mean")

x1 <- x; x1[sample(1:100,10)] <- NA
myEst1 <- Estimator(x1, estimator = mean, name = "mean")
samplesize(myEst1)
samplesize(myEst1, onlycomplete = FALSE)



cleanEx()
nameEx("Estimator")
### * Estimator

flush(stderr()); flush(stdout())

### Name: Estimator
### Title: Function to compute estimates
### Aliases: Estimator
### Keywords: univar

### ** Examples

x <- rnorm(100)
Estimator(x, estimator = mean, name = "mean")

X <- matrix(rnorm(1000), nrow = 10)
Estimator(X, estimator = rowMeans, name = "mean")



cleanEx()
nameEx("EvenSymmetric-class")
### * EvenSymmetric-class

flush(stderr()); flush(stdout())

### Name: EvenSymmetric-class
### Title: Class for Even Functions
### Aliases: EvenSymmetric-class
### Keywords: classes

### ** Examples

new("EvenSymmetric")



cleanEx()
nameEx("EvenSymmetric")
### * EvenSymmetric

flush(stderr()); flush(stdout())

### Name: EvenSymmetric
### Title: Generating function for EvenSymmetric-class
### Aliases: EvenSymmetric
### Keywords: math

### ** Examples

EvenSymmetric()

## The function is currently defined as
function(SymmCenter = 0){ 
    new("EvenSymmetric", SymmCenter = SymmCenter) 
}



cleanEx()
nameEx("ExpScaleFamily")
### * ExpScaleFamily

flush(stderr()); flush(stdout())

### Name: ExpScaleFamily
### Title: Generating function for exponential scale families
### Aliases: ExpScaleFamily
### Keywords: models

### ** Examples

(E1 <- ExpScaleFamily())
plot(E1)
Map(L2deriv(E1)[[1]])
checkL2deriv(E1)



cleanEx()
nameEx("FunSymmList-class")
### * FunSymmList-class

flush(stderr()); flush(stdout())

### Name: FunSymmList-class
### Title: List of Symmetries for a List of Functions
### Aliases: FunSymmList-class
### Keywords: classes

### ** Examples

new("FunSymmList", list(NonSymmetric(), EvenSymmetric(SymmCenter = 1), 
                        OddSymmetric(SymmCenter = 2)))



cleanEx()
nameEx("FunSymmList")
### * FunSymmList

flush(stderr()); flush(stdout())

### Name: FunSymmList
### Title: Generating function for FunSymmList-class
### Aliases: FunSymmList
### Keywords: math

### ** Examples

FunSymmList(NonSymmetric(), EvenSymmetric(SymmCenter = 1), 
            OddSymmetric(SymmCenter = 2))

## The function is currently defined as
function (...){
    new("FunSymmList", list(...))
}



cleanEx()
nameEx("GammaFamily")
### * GammaFamily

flush(stderr()); flush(stdout())

### Name: GammaFamily
### Title: Generating function for Gamma families
### Aliases: GammaFamily
### Keywords: models

### ** Examples

(G1 <- GammaFamily())
FisherInfo(G1)
checkL2deriv(G1)



cleanEx()
nameEx("InfoNorm")
### * InfoNorm

flush(stderr()); flush(stdout())

### Name: InfoNorm
### Title: Generating function for InfoNorm-class
### Aliases: InfoNorm
### Keywords: robust

### ** Examples

InfoNorm()

## The function is currently defined as
function(){ new("InfoNorm") }



cleanEx()
nameEx("L2GroupFamily-class")
### * L2GroupFamily-class

flush(stderr()); flush(stdout())

### Name: L2GroupParamFamily-class
### Title: L2 differentiable parametric group family
### Aliases: L2GroupParamFamily-class LogDeriv
###   LogDeriv,L2GroupParamFamily-method LogDeriv<-
###   LogDeriv<-,L2GroupParamFamily-method
### Keywords: classes models

### ** Examples

F1 <- new("L2GroupParamFamily")
plot(F1)



cleanEx()
nameEx("L2LocationFamily-class")
### * L2LocationFamily-class

flush(stderr()); flush(stdout())

### Name: L2LocationFamily-class
### Title: L2 differentiable parametric group family
### Aliases: L2LocationFamily-class
### Keywords: classes models

### ** Examples

F1 <- new("L2LocationFamily")
plot(F1)



cleanEx()
nameEx("L2LocationFamily")
### * L2LocationFamily

flush(stderr()); flush(stdout())

### Name: L2LocationFamily
### Title: Generating function for L2LocationFamily-class
### Aliases: L2LocationFamily
### Keywords: models

### ** Examples

F1 <- L2LocationFamily()
plot(F1)



cleanEx()
nameEx("L2LocationScaleFamily-class")
### * L2LocationScaleFamily-class

flush(stderr()); flush(stdout())

### Name: L2LocationScaleFamily-class
### Title: L2 differentiable parametric group family
### Aliases: L2LocationScaleFamily-class
### Keywords: classes models

### ** Examples

F1 <- new("L2LocationScaleFamily")
plot(F1)



cleanEx()
nameEx("L2LocationScaleFamily")
### * L2LocationScaleFamily

flush(stderr()); flush(stdout())

### Name: L2LocationScaleFamily
### Title: Generating function for L2LocationScaleFamily-class
### Aliases: L2LocationScaleFamily
### Keywords: models

### ** Examples

F1 <- L2LocationScaleFamily()
plot(F1)



cleanEx()
nameEx("L2LocationUnknownScaleFamily")
### * L2LocationUnknownScaleFamily

flush(stderr()); flush(stdout())

### Name: L2LocationUnknownScaleFamily
### Title: Generating function for L2LocationScaleFamily-class in nuisance
###   situation
### Aliases: L2LocationUnknownScaleFamily
### Keywords: models

### ** Examples

F1 <- L2LocationUnknownScaleFamily()
plot(F1)



cleanEx()
nameEx("L2ParamFamily-class")
### * L2ParamFamily-class

flush(stderr()); flush(stdout())

### Name: L2ParamFamily-class
### Title: L2 differentiable parametric family
### Aliases: plot plot-methods L2ParamFamily-class FisherInfo
###   FisherInfo,L2ParamFamily,missing-method
###   FisherInfo,L2ParamFamily,ParamFamParameter-method L2deriv
###   L2deriv,L2ParamFamily,missing-method
###   L2deriv,L2ParamFamily,ParamFamParameter-method L2derivSymm
###   L2derivSymm,L2ParamFamily-method L2derivDistr
###   L2derivDistr,L2ParamFamily-method L2derivDistrSymm
###   L2derivDistrSymm,L2ParamFamily-method
###   checkL2deriv,L2ParamFamily-method
###   E,L2ParamFamily,EuclRandVariable,missing-method
###   E,L2ParamFamily,EuclRandMatrix,missing-method
###   E,L2ParamFamily,EuclRandVarList,missing-method
###   plot,L2ParamFamily,missing-method
### Keywords: classes models

### ** Examples

F1 <- new("L2ParamFamily")
plot(F1)

## selection of subpanels for plotting
F2 <- L2LocationScaleFamily()
layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
plot(F2,mfColRow = FALSE,
     to.draw.arg=c("p","q","loc"))
plot(F2,mfColRow = FALSE, inner=list("empirical cdf","pseudo-inverse",
     "L2-deriv, loc.part"), to.draw.arg=c("p","q","loc"))



cleanEx()
nameEx("L2ParamFamily")
### * L2ParamFamily

flush(stderr()); flush(stdout())

### Name: L2ParamFamily
### Title: Generating function for L2ParamFamily-class
### Aliases: L2ParamFamily
### Keywords: models

### ** Examples

F1 <- L2ParamFamily()
plot(F1)



cleanEx()
nameEx("L2ScaleFamily-class")
### * L2ScaleFamily-class

flush(stderr()); flush(stdout())

### Name: L2ScaleFamily-class
### Title: L2 differentiable parametric group family
### Aliases: L2ScaleFamily-class
### Keywords: classes models

### ** Examples

F1 <- new("L2ScaleFamily")
plot(F1)



cleanEx()
nameEx("L2ScaleFamily")
### * L2ScaleFamily

flush(stderr()); flush(stdout())

### Name: L2ScaleFamily
### Title: Generating function for L2ScaleFamily-class
### Aliases: L2ScaleFamily
### Keywords: models

### ** Examples

F1 <- L2ScaleFamily()
plot(F1)



cleanEx()
nameEx("L2ScaleUnknownLocationFamily")
### * L2ScaleUnknownLocationFamily

flush(stderr()); flush(stdout())

### Name: L2ScaleUnknownLocationFamily
### Title: Generating function for L2LocationScaleFamily-class in nuisance
###   situation
### Aliases: L2ScaleUnknownLocationFamily
### Keywords: models

### ** Examples

F1 <- L2ScaleUnknownLocationFamily()
plot(F1)



cleanEx()
nameEx("LnormScaleFamily")
### * LnormScaleFamily

flush(stderr()); flush(stdout())

### Name: LnormScaleFamily
### Title: Generating function for lognormal scale families
### Aliases: LnormScaleFamily
### Keywords: models

### ** Examples

(L1 <- LnormScaleFamily())
plot(L1)
Map(L2deriv(L1)[[1]])
checkL2deriv(L1)



cleanEx()
nameEx("MCEstimate-class")
### * MCEstimate-class

flush(stderr()); flush(stdout())

### Name: MCEstimate-class
### Title: MCEstimate-class.
### Aliases: MCEstimate-class criterion criterion,MCEstimate-method
###   criterion.fct criterion.fct,MCEstimate-method
###   startPar,MCEstimate-method method method,MCEstimate-method optimwarn
###   optimwarn,MCEstimate-method criterion<- criterion<-,MCEstimate-method
###   coerce,MCEstimate,mle-method show,MCEstimate-method
###   profile,MCEstimate-method
### Keywords: classes

### ** Examples

## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

MDEstimator(x, G)
(m <- MLEstimator(x, G))
m.mle <- as(m,"mle")
par(mfrow=c(1,2))
profileM <- profile(m)
## plot-profile throws an error



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("MCEstimator")
### * MCEstimator

flush(stderr()); flush(stdout())

### Name: MCEstimator
### Title: Function to compute minimum criterion estimates
### Aliases: MCEstimator
### Keywords: univar

### ** Examples

## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Maximum Likelihood estimator
## Note: you can directly use function MLEstimator!
negLoglikelihood <- function(x, Distribution){
    res <- -sum(log(Distribution@d(x)))
    names(res) <- "Negative Log-Likelihood"
    return(res)
}
MCEstimator(x = x, ParamFamily = G, criterion = negLoglikelihood)

## Kolmogorov(-Smirnov) minimum distance estimator
## Note: you can also use function MDEstimator!
MCEstimator(x = x, ParamFamily = G, criterion = KolmogorovDist, 
            crit.name = "Kolmogorov distance")

## Total variation minimum distance estimator
## Note: you can also use function MDEstimator!
## discretize Gamma distribution
MCEstimator(x = x, ParamFamily = G, criterion = TotalVarDist, 
            crit.name = "Total variation distance")

## or smooth empirical distribution (takes some time!)
#MCEstimator(x = x, ParamFamily = G, criterion = TotalVarDist, 
#            asis.smooth.discretize = "smooth", crit.name = "Total variation distance")

## Hellinger minimum distance estimator
## Note: you can also use function MDEstimator!
## discretize Gamma distribution
distroptions(DistrResolution = 1e-8)
MCEstimator(x = x, ParamFamily = G, criterion = HellingerDist, 
            crit.name = "Hellinger Distance", startPar = c(1,2))
distroptions(DistrResolution = 1e-6)

## or smooth empirical distribution (takes some time!)
#MCEstimator(x = x, ParamFamily = G, criterion = HellingerDist, 
#            asis.smooth.discretize = "smooth", crit.name = "Hellinger distance")



cleanEx()
nameEx("MDEstimator")
### * MDEstimator

flush(stderr()); flush(stdout())

### Name: MDEstimator
### Title: Function to compute minimum distance estimates
### Aliases: MDEstimator
### Keywords: univar robust

### ** Examples

## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Kolmogorov(-Smirnov) minimum distance estimator
MDEstimator(x = x, ParamFamily = G, distance = KolmogorovDist)

## von Mises minimum distance estimator with default mu
MDEstimator(x = x, ParamFamily = G, distance = CvMDist)

## don't run to reduce check time on CRAN
## Not run: 
##D ## von Mises minimum distance estimator with default mu
##D MDEstimator(x = x, ParamFamily = G, distance = CvMDist,
##D             asvar.fct = distrMod:::.CvMMDCovariance)
##D #*** variance routine is still in testing phase so not yet
##D #*** exported to namespace
##D ## von Mises minimum distance estimator with mu = N(0,1)
##D MDEstimator(x = x, ParamFamily = G, distance = CvMDist, mu = Norm())
##D 
##D ## Total variation minimum distance estimator
##D ## gamma distributions are discretized
##D MDEstimator(x = x, ParamFamily = G, distance = TotalVarDist)
##D ## or smoothing of emprical distribution (takes some time!)
##D #MDEstimator(x = x, ParamFamily = G, distance = TotalVarDist, asis.smooth.discretize = "smooth")
##D 
##D ## Hellinger minimum distance estimator
##D ## gamma distributions are discretized
##D distroptions(DistrResolution = 1e-10)
##D MDEstimator(x = x, ParamFamily = G, distance = HellingerDist, startPar = c(1,2))
##D distroptions(DistrResolution = 1e-6) # default
##D ## or smoothing of emprical distribution (takes some time!)
##D #MDEstimator(x = x, ParamFamily = G, distance = HellingerDist, asis.smooth.discretize = "smooth")
## End(Not run)



cleanEx()
nameEx("MLEstimator")
### * MLEstimator

flush(stderr()); flush(stdout())

### Name: MLEstimator
### Title: Function to compute maximum likelihood estimates
### Aliases: MLEstimator
### Keywords: univar

### ** Examples

#############################
## 1. Binomial data
#############################
## (empirical) data
x <- rbinom(100, size=25, prob=.25)

## ML-estimate
MLEstimator(x, BinomFamily(size = 25))


#############################
## 2. Poisson data
#############################
## Example: Rutherford-Geiger (1910); cf. Feller~(1968), Section VI.7 (a)
x <- c(rep(0, 57), rep(1, 203), rep(2, 383), rep(3, 525), rep(4, 532), 
       rep(5, 408), rep(6, 273), rep(7, 139), rep(8, 45), rep(9, 27), 
       rep(10, 10), rep(11, 4), rep(12, 0), rep(13, 1), rep(14, 1))

## ML-estimate
MLEstimator(x, PoisFamily())


#############################
## 3. Normal (Gaussian) location and scale
#############################
## (empirical) data
x <- rnorm(100)

## ML-estimate
MLEstimator(x, NormLocationScaleFamily())
## compare:
c(mean(x),sd(x))


#############################
## 4. Gamma model
#############################
## (empirical) data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Maximum likelihood estimator
(res <- MLEstimator(x = x, ParamFamily = G))

## Asymptotic (CLT-based) confidence interval
confint(res)

## some profiling
par(mfrow=c(1,2))
plot(profile(res))
par(mfrow=c(1,1))

## implementation of ML-estimator of package MASS
require(MASS)
(res1 <- fitdistr(x, "gamma"))

## comparison
## shape
estimate(res)[2]
## rate
1/estimate(res)[1]

## minor differences due to the fact that by default, fitdistr uses
## BFGS, while we use Nelder-Mead instead

## log-likelihood
res1$loglik
## negative log-likelihood
criterion(res)


## don't run to reduce check time on CRAN
## Not run: 
##D ## explicitely transforming to
##D ## MASS parametrization:
##D mtrafo <- function(x){
##D      nms0 <- names(c(main(param(G)),nuisance(param(G))))
##D      nms <- c("shape","rate")
##D      fval0 <- c(x[2], 1/x[1])
##D      names(fval0) <- nms
##D      mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
##D                      dimnames = list(nms,nms0))                          
##D      list(fval = fval0, mat = mat0)}
##D 
##D G2 <- G
##D trafo(G2) <- mtrafo
##D res2 <- MLEstimator(x = x, ParamFamily = G2)
##D 
##D old <- getdistrModOption("show.details")
##D distrModoptions("show.details" = "minimal")
##D res1
##D res2
##D 
##D ## some profiling
##D par(mfrow=c(1,2))
##D plot(profile(res2))
##D par(mfrow=c(1,1))
## End(Not run)

#############################
## 5. Cauchy Location Scale model
#############################
(C <- CauchyLocationScaleFamily())
loc.true <- 1
scl.true <- 2

## (empirical) data
x <- rcauchy(50, location = loc.true, scale = scl.true)

## Maximum likelihood estimator
(res <- MLEstimator(x = x, ParamFamily = C))
## Asymptotic (CLT-based) confidence interval
confint(res)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("NBinomFamily")
### * NBinomFamily

flush(stderr()); flush(stdout())

### Name: NbinomFamily
### Title: Generating function for Nbinomial families
### Aliases: NbinomFamily NbinomwithSizeFamily NbinomMeanSizeFamily
### Keywords: models

### ** Examples

(N1 <- NbinomFamily(size = 25, prob = 0.25))
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)
(N1.w <- NbinomwithSizeFamily(size = 25, prob = 0.25))
plot(N1.w)
FisherInfo(N1.w)
checkL2deriv(N1.w)
(N2.w <- NbinomMeanSizeFamily(size = 25, mean = 75))
plot(N2.w)
FisherInfo(N2.w)
checkL2deriv(N2.w)




cleanEx()
nameEx("NonSymmetric-class")
### * NonSymmetric-class

flush(stderr()); flush(stdout())

### Name: NonSymmetric-class
### Title: Class for Non-symmetric Functions
### Aliases: NonSymmetric-class
### Keywords: classes

### ** Examples

new("NonSymmetric")



cleanEx()
nameEx("NonSymmetric")
### * NonSymmetric

flush(stderr()); flush(stdout())

### Name: NonSymmetric
### Title: Generating function for NonSymmetric-class
### Aliases: NonSymmetric
### Keywords: math

### ** Examples

NonSymmetric()

## The function is currently defined as
function(){ new("NonSymmetric") }



cleanEx()
nameEx("NormLocationFamily")
### * NormLocationFamily

flush(stderr()); flush(stdout())

### Name: NormLocationFamily
### Title: Generating function for normal location families
### Aliases: NormLocationFamily
### Keywords: models

### ** Examples

(N1 <- NormLocationFamily())
plot(N1)
L2derivDistr(N1)



cleanEx()
nameEx("NormLocationScaleFamily")
### * NormLocationScaleFamily

flush(stderr()); flush(stdout())

### Name: NormLocationScaleFamily
### Title: Generating function for normal location and scale families
### Aliases: NormLocationScaleFamily
### Keywords: models

### ** Examples

(N1 <- NormLocationScaleFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)



cleanEx()
nameEx("NormLocationUnknownScaleFamily")
### * NormLocationUnknownScaleFamily

flush(stderr()); flush(stdout())

### Name: NormLocationUnknownScaleFamily
### Title: Generating function for normal location families with unknown
###   scale as nuisance
### Aliases: NormLocationUnknownScaleFamily
### Keywords: models

### ** Examples

(N1 <- NormLocationUnknownScaleFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)



cleanEx()
nameEx("NormScaleFamily")
### * NormScaleFamily

flush(stderr()); flush(stdout())

### Name: NormScaleFamily
### Title: Generating function for normal scale families
### Aliases: NormScaleFamily
### Keywords: models

### ** Examples

(N1 <- NormScaleFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)



cleanEx()
nameEx("NormScaleUnknownLocationFamily")
### * NormScaleUnknownLocationFamily

flush(stderr()); flush(stdout())

### Name: NormScaleUnknownLocationFamily
### Title: Generating function for normal scale families with unknown
###   location as nuisance
### Aliases: NormScaleUnknownLocationFamily
### Keywords: models

### ** Examples

(N1 <- NormScaleUnknownLocationFamily())
plot(N1)
FisherInfo(N1)
checkL2deriv(N1)



cleanEx()
nameEx("NormType-class")
### * NormType-class

flush(stderr()); flush(stdout())

### Name: NormType-class
### Title: Norm Type
### Aliases: NormType-class name,NormType-method name<-,NormType-method fct
###   fct<- fct,NormType-method fct<-,NormType-method
### Keywords: classes

### ** Examples

EuclNorm <- NormType("EuclideanNorm",EuclideanNorm)
fct(EuclNorm)
name(EuclNorm)



cleanEx()
nameEx("NormType")
### * NormType

flush(stderr()); flush(stdout())

### Name: NormType
### Title: Generating function for NormType-class
### Aliases: NormType
### Keywords: math

### ** Examples

NormType()



cleanEx()
nameEx("OddSymmetric-class")
### * OddSymmetric-class

flush(stderr()); flush(stdout())

### Name: OddSymmetric-class
### Title: Class for Odd Functions
### Aliases: OddSymmetric-class
### Keywords: classes

### ** Examples

new("OddSymmetric")



cleanEx()
nameEx("OddSymmetric")
### * OddSymmetric

flush(stderr()); flush(stdout())

### Name: OddSymmetric
### Title: Generating function for OddSymmetric-class
### Aliases: OddSymmetric
### Keywords: math

### ** Examples

OddSymmetric()

## The function is currently defined as
function(SymmCenter = 0){ 
    new("OddSymmetric", SymmCenter = SymmCenter) 
}



cleanEx()
nameEx("ParamFamParameter-class")
### * ParamFamParameter-class

flush(stderr()); flush(stdout())

### Name: ParamFamParameter-class
### Title: Parameter of a parametric family of probability measures
### Aliases: ParamFamParameter-class ParamWithScaleFamParameter-class
###   ParamWithScaleAndShapeFamParameter-class
###   ParamWithShapeFamParameter-class length,ParamFamParameter-method
###   dimension,ParamFamParameter-method main main,ParamFamParameter-method
###   main,ParamWithScaleAndShapeFamParameter-method main<-
###   main<-,ParamFamParameter-method nuisance
###   nuisance,ParamFamParameter-method
###   nuisance,ParamWithScaleAndShapeFamParameter-method nuisance<-
###   nuisance<-,ParamFamParameter-method fixed
###   fixed,ParamFamParameter-method
###   fixed,ParamWithScaleAndShapeFamParameter-method fixed<-
###   fixed<-,ParamFamParameter-method withPosRestr
###   withPosRestr,ParamWithShapeFamParameter-method withPosRestr<-
###   withPosRestr<-,ParamWithShapeFamParameter-method
###   show,ParamFamParameter-method show,ParamWithShapeFamParameter-method
###   show,ParamWithScaleAndShapeFamParameter-method
### Keywords: classes

### ** Examples

new("ParamFamParameter")



cleanEx()
nameEx("ParamFamParameter")
### * ParamFamParameter

flush(stderr()); flush(stdout())

### Name: ParamFamParameter
### Title: Generating function for ParamFamParameter-class
### Aliases: ParamFamParameter
### Keywords: models

### ** Examples

ParamFamParameter(main = 0, nuisance = 1, fixed = 2,
                  trafo = function(x) list(fval = sin(x), 
                                            mat = matrix(cos(x),1,1))
                  )                          




cleanEx()
nameEx("ParamFamily-class")
### * ParamFamily-class

flush(stderr()); flush(stdout())

### Name: ParamFamily-class
### Title: Parametric family of probability measures.
### Aliases: ParamFamily-class main,ParamFamily-method
###   nuisance,ParamFamily-method fixed,ParamFamily-method
###   param,ParamFamily-method modifyParam modifyParam,ParamFamily-method
###   fam.call fam.call,ParamFamily-method startPar
###   startPar,ParamFamily-method makeOKPar makeOKPar,ParamFamily-method
###   plot,ParamFamily,missing-method show,ParamFamily-method
### Keywords: classes models

### ** Examples

F1 <- new("ParamFamily") # prototype
plot(F1)



cleanEx()
nameEx("ParamFamily")
### * ParamFamily

flush(stderr()); flush(stdout())

### Name: ParamFamily
### Title: Generating function for ParamFamily-class
### Aliases: ParamFamily
### Keywords: distribution models

### ** Examples


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



cleanEx()
nameEx("PoisFamily")
### * PoisFamily

flush(stderr()); flush(stdout())

### Name: PoisFamily
### Title: Generating function for Poisson families
### Aliases: PoisFamily
### Keywords: models

### ** Examples

(P1 <- PoisFamily(lambda = 4.5))
plot(P1)
FisherInfo(P1)
checkL2deriv(P1)



cleanEx()
nameEx("QFNorm")
### * QFNorm

flush(stderr()); flush(stdout())

### Name: QFNorm
### Title: Generating function for QFNorm-class
### Aliases: QFNorm
### Keywords: math

### ** Examples

QFNorm()

## The function is currently defined as
function(){ new("QFNorm") }



cleanEx()
nameEx("SelfNorm")
### * SelfNorm

flush(stderr()); flush(stdout())

### Name: SelfNorm
### Title: Generating function for SelfNorm-class
### Aliases: SelfNorm
### Keywords: robust

### ** Examples

SelfNorm()

## The function is currently defined as
function(){ new("SelfNorm") }



cleanEx()
nameEx("addAlphTrsp2col")
### * addAlphTrsp2col

flush(stderr()); flush(stdout())

### Name: addAlphTrsp2col
### Title: "addAlphTrsp2col"
### Aliases: addAlphTrsp2col
### Keywords: distribution

### ** Examples

  addAlphTrsp2col(rgb(1,0.3,0.03), 25)
  addAlphTrsp2col("darkblue", 25)
  addAlphTrsp2col("#AAAAAAAA",25)
  palette(rainbow(6))
  addAlphTrsp2col(2, 25)



cleanEx()
nameEx("asBias-class")
### * asBias-class

flush(stderr()); flush(stdout())

### Name: asBias-class
### Title: Standardized Asymptotic Bias
### Aliases: asBias-class
### Keywords: classes

### ** Examples

new("asBias")



cleanEx()
nameEx("asBias")
### * asBias

flush(stderr()); flush(stdout())

### Name: asBias
### Title: Generating function for asBias-class
### Aliases: asBias
### Keywords: robust

### ** Examples

asBias()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
     new("asBias",biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("asCov-class")
### * asCov-class

flush(stderr()); flush(stdout())

### Name: asCov-class
### Title: Asymptotic covariance
### Aliases: asCov-class
### Keywords: classes

### ** Examples

new("asCov")



cleanEx()
nameEx("asCov")
### * asCov

flush(stderr()); flush(stdout())

### Name: asCov
### Title: Generating function for asCov-class
### Aliases: asCov
### Keywords: robust

### ** Examples

asCov()

## The function is currently defined as
function(){ new("asCov") }



cleanEx()
nameEx("asHampel-class")
### * asHampel-class

flush(stderr()); flush(stdout())

### Name: asHampel-class
### Title: Asymptotic Hampel risk
### Aliases: asHampel-class bound bound,asHampel-method
###   show,asHampel-method
### Keywords: classes

### ** Examples

new("asHampel")



cleanEx()
nameEx("asHampel")
### * asHampel

flush(stderr()); flush(stdout())

### Name: asHampel
### Title: Generating function for asHampel-class
### Aliases: asHampel
### Keywords: robust

### ** Examples

asHampel()

## The function is currently defined as
function(bound = Inf, biastype = symmetricBias(), normtype = NormType()){ 
    new("asHampel", bound = bound, biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("asMSE-class")
### * asMSE-class

flush(stderr()); flush(stdout())

### Name: asMSE-class
### Title: Asymptotic mean square error
### Aliases: asMSE-class
### Keywords: classes

### ** Examples

new("asMSE")



cleanEx()
nameEx("asMSE")
### * asMSE

flush(stderr()); flush(stdout())

### Name: asMSE
### Title: Generating function for asMSE-class
### Aliases: asMSE
### Keywords: robust

### ** Examples

asMSE()

## The function is currently defined as
function(biastype = symmetricBias(), normtype = NormType()){ 
         new("asMSE", biastype = biastype, normtype = normtype) }



cleanEx()
nameEx("asSemivar-class")
### * asSemivar-class

flush(stderr()); flush(stdout())

### Name: asSemivar-class
### Title: Semivariance Risk Type
### Aliases: asSemivar-class sign,asSemivar-method sign<-,asSemivar-method
### Keywords: classes

### ** Examples

asSemivar()



cleanEx()
nameEx("asSemivar")
### * asSemivar

flush(stderr()); flush(stdout())

### Name: asSemivar
### Title: Generating function for asSemivar-class
### Aliases: asSemivar
### Keywords: robust

### ** Examples

asSemivar()




cleanEx()
nameEx("asUnOvShoot-class")
### * asUnOvShoot-class

flush(stderr()); flush(stdout())

### Name: asUnOvShoot-class
### Title: Asymptotic under-/overshoot probability
### Aliases: asUnOvShoot-class width width,asUnOvShoot-method
###   show,asUnOvShoot-method
### Keywords: classes

### ** Examples

new("asUnOvShoot")



cleanEx()
nameEx("asUnOvShoot")
### * asUnOvShoot

flush(stderr()); flush(stdout())

### Name: asUnOvShoot
### Title: Generating function for asUnOvShoot-class
### Aliases: asUnOvShoot
### Keywords: robust

### ** Examples

asUnOvShoot()

## The function is currently defined as
function(width = 1.960, biastype = symmetricBias()){ 
     new("asUnOvShoot", width = width, biastype = biastype) }



cleanEx()
nameEx("asymmetricBias-class")
### * asymmetricBias-class

flush(stderr()); flush(stdout())

### Name: asymmetricBias-class
### Title: asymmetric Bias Type
### Aliases: asymmetricBias-class nu,asymmetricBias-method
###   nu<-,asymmetricBias-method nu nu<-
### Keywords: classes

### ** Examples

asymmetricBias()
## The function is currently defined as
function(){ new("asymmetricBias", name = "asymmetric Bias", nu = c(1,1)) }

aB <- asymmetricBias()
nu(aB)
try(nu(aB) <- -2) ## error
nu(aB) <- c(0.3,1)



cleanEx()
nameEx("asymmetricBias")
### * asymmetricBias

flush(stderr()); flush(stdout())

### Name: asymmetricBias
### Title: Generating function for asymmetricBias-class
### Aliases: asymmetricBias
### Keywords: robust

### ** Examples

asymmetricBias()

## The function is currently defined as
function(){ new("asymmetricBias", name = "asymmetric Bias", nu = c(1,1)) }



cleanEx()
nameEx("checkL2deriv")
### * checkL2deriv

flush(stderr()); flush(stdout())

### Name: checkL2deriv
### Title: Generic function for checking L2-derivatives
### Aliases: checkL2deriv
### Keywords: models

### ** Examples

F1 <- new("L2ParamFamily")
checkL2deriv(F1)



cleanEx()
nameEx("confint-methods")
### * confint-methods

flush(stderr()); flush(stdout())

### Name: confint-methods
### Title: Methods for function confint in Package 'distrMod'
### Aliases: confint-methods confint confint,ANY,missing-method
###   confint,Estimate,missing-method confint,mle,missing-method
###   confint,profile.mle,missing-method
### Keywords: models

### ** Examples

## for signature ANY examples confer stats::confint
## (empirical) Data
x <- rgamma(50, scale = 0.5, shape = 3)

## parametric family of probability measures
G <- GammaFamily(scale = 1, shape = 2)

## Maximum likelihood estimator
res <- MLEstimator(x = x, ParamFamily = G)
confint(res)

### for comparison:
require(MASS)
(res1 <- fitdistr(x, "gamma"))
## add a convenient (albeit wrong)
## S3-method for vcov:
## --- wrong as in general cov-matrix
##     will not be diagonal
## but for conf-interval this does
## not matter...
vcov.fitdistr <- function(object, ...){
     v<-diag(object$sd^2)
     rownames(v) <- colnames(v) <- names(object$estimate) 
     v}

## explicitely transforming to
## MASS parametrization:
mtrafo <- function(x){
     nms0 <- names(c(main(param(G)),nuisance(param(G))))
     nms <- c("shape","rate")
     fval0 <- c(x[2], 1/x[1])
     names(fval0) <- nms
     mat0 <- matrix( c(0, -1/x[1]^2, 1, 0), nrow = 2, ncol = 2,
                     dimnames = list(nms,nms0))                          
     list(fval = fval0, mat = mat0)}

G2 <- G
trafo(G2) <- mtrafo
res2 <- MLEstimator(x = x, ParamFamily = G2)

old<-getdistrModOption("show.details")
distrModoptions("show.details" = "minimal")
res
res1
res2
confint(res)
confint(res1)
confint(res2)
confint(res,level=0.99)
distrModoptions("show.details" = old)
 



cleanEx()
nameEx("distrModMASK")
### * distrModMASK

flush(stderr()); flush(stdout())

### Name: distrModMASK
### Title: Masking of/by other functions in package "distrMod"
### Aliases: distrModMASK MASKING
### Keywords: programming distribution documentation

### ** Examples

distrModMASK()



cleanEx()
nameEx("distrModOptions")
### * distrModOptions

flush(stderr()); flush(stdout())

### Name: distrModOptions
### Title: Function to change the global variables of the package
###   'distrMod'
### Aliases: distrModOptions getdistrModOption distrModoptions show.details
### Keywords: misc distribution

### ** Examples

distrModoptions()
distrModoptions("show.details")
distrModoptions("show.details" = "maximal")
distrModOptions("show.details" = "minimal")
# or
getdistrModOption("show.details")



cleanEx()
nameEx("fiBias-class")
### * fiBias-class

flush(stderr()); flush(stdout())

### Name: fiBias-class
### Title: Finite-sample Bias
### Aliases: fiBias-class
### Keywords: classes

### ** Examples

new("fiBias")



cleanEx()
nameEx("fiBias")
### * fiBias

flush(stderr()); flush(stdout())

### Name: fiBias
### Title: Generating function for fiBias-class
### Aliases: fiBias
### Keywords: robust

### ** Examples

fiBias()

## The function is currently defined as
function(){ new("fiBias") }



cleanEx()
nameEx("fiCov-class")
### * fiCov-class

flush(stderr()); flush(stdout())

### Name: fiCov-class
### Title: Finite-sample covariance
### Aliases: fiCov-class
### Keywords: classes

### ** Examples

new("fiCov")



cleanEx()
nameEx("fiCov")
### * fiCov

flush(stderr()); flush(stdout())

### Name: fiCov
### Title: Generating function for fiCov-class
### Aliases: fiCov
### Keywords: robust

### ** Examples

fiCov()

## The function is currently defined as
function(){ new("fiCov") }



cleanEx()
nameEx("fiHampel-class")
### * fiHampel-class

flush(stderr()); flush(stdout())

### Name: fiHampel-class
### Title: Finite-sample Hampel risk
### Aliases: fiHampel-class bound,fiHampel-method show,fiHampel-method
### Keywords: classes

### ** Examples

new("fiHampel")



cleanEx()
nameEx("fiHampel")
### * fiHampel

flush(stderr()); flush(stdout())

### Name: fiHampel
### Title: Generating function for fiHampel-class
### Aliases: fiHampel
### Keywords: robust

### ** Examples

fiHampel()

## The function is currently defined as
function(bound = Inf){ new("fiHampel", bound = bound) }



cleanEx()
nameEx("fiMSE-class")
### * fiMSE-class

flush(stderr()); flush(stdout())

### Name: fiMSE-class
### Title: Finite-sample mean square error
### Aliases: fiMSE-class
### Keywords: classes

### ** Examples

new("fiMSE")



cleanEx()
nameEx("fiMSE")
### * fiMSE

flush(stderr()); flush(stdout())

### Name: fiMSE
### Title: Generating function for fiMSE-class
### Aliases: fiMSE
### Keywords: robust

### ** Examples

fiMSE()

## The function is currently defined as
function(){ new("fiMSE") }



cleanEx()
nameEx("fiUnOvShoot-class")
### * fiUnOvShoot-class

flush(stderr()); flush(stdout())

### Name: fiUnOvShoot-class
### Title: Finite-sample under-/overshoot probability
### Aliases: fiUnOvShoot-class width,fiUnOvShoot-method
###   show,fiUnOvShoot-method
### Keywords: classes

### ** Examples

new("fiUnOvShoot")



cleanEx()
nameEx("fiUnOvShoot")
### * fiUnOvShoot

flush(stderr()); flush(stdout())

### Name: fiUnOvShoot
### Title: Generating function for fiUnOvShoot-class
### Aliases: fiUnOvShoot
### Keywords: robust

### ** Examples

fiUnOvShoot()

## The function is currently defined as
function(width = 1.960){ new("fiUnOvShoot", width = width) }



cleanEx()
nameEx("isKerAinKerB")
### * isKerAinKerB

flush(stderr()); flush(stdout())

### Name: isKerAinKerB
### Title: isKerAinKerB
### Aliases: isKerAinKerB
### Keywords: algebra array

### ** Examples

ma <- cbind(1,1,c(1,1,7))
D <- t(ma %*% c(0,1,-1))
isKerAinKerB(D,ma)
isKerAinKerB(ma,D)



cleanEx()
nameEx("negativeBias")
### * negativeBias

flush(stderr()); flush(stdout())

### Name: negativeBias
### Title: Generating function for onesidedBias-class
### Aliases: negativeBias
### Keywords: robust

### ** Examples

negativeBias()

## The function is currently defined as
function(){ new("onesidedBias", name = "negative Bias", sign = -1) }



cleanEx()
nameEx("norms")
### * norms

flush(stderr()); flush(stdout())

### Name: norm
### Title: Norm functions
### Aliases: EuclideanNorm QuadFormNorm
### Keywords: robust

### ** Examples

mm <- matrix(rnorm(20),2,10)
EuclideanNorm(mm)
QuadFormNorm(mm, A = PosSemDefSymmMatrix(matrix(c(3,1,1,1),2,2)))



cleanEx()
nameEx("onesidedBias-class")
### * onesidedBias-class

flush(stderr()); flush(stdout())

### Name: onesidedBias-class
### Title: onesided Bias Type
### Aliases: onesidedBias-class sign sign<- sign,onesidedBias-method
###   sign<-,onesidedBias-method
### Keywords: classes

### ** Examples

positiveBias()
## The function is currently defined as
function(){ new("onesidedBias", name = "positive Bias", sign = 1) }

negativeBias()
## The function is currently defined as
function(){ new("onesidedBias", name = "negative Bias", sign = -1) }

pB <- positiveBias()
sign(pB)
try(sign(pB) <- -2) ## error
sign(pB) <- -1



cleanEx()
nameEx("positiveBias")
### * positiveBias

flush(stderr()); flush(stdout())

### Name: positiveBias
### Title: Generating function for onesidedBias-class
### Aliases: positiveBias
### Keywords: robust

### ** Examples

positiveBias()

## The function is currently defined as
function(){ new("onesidedBias", name = "positive Bias", sign = 1) }



cleanEx()
nameEx("print-methods")
### * print-methods

flush(stderr()); flush(stdout())

### Name: print-methods
### Title: Common 'print' Methods for S4 classes in Package 'distrMod'
### Aliases: print-methods print,ShowDetails-method
### Keywords: models

### ** Examples

## set options to maximal detailedness
show.old <- getdistrModOption("show.details")
distrModoptions("show.details" = "maximal")
## define a model
NS <- NormLocationScaleFamily(mean=2, sd=3)
## generate data out of this situation
x <- r(distribution(NS))(30)

## want to estimate mu/sigma, sigma^2
## -> new trafo slot:
trafo(NS) <- function(param){
  mu <- param["mean"]
  sd <- param["sd"]
  fval <- c(mu/sd, sd^2)
  nfval <- c("mu/sig", "sig^2")
  names(fval) <- nfval
  mat <- matrix(c(1/sd,0,-mu/sd^2,2*sd),2,2)
  dimnames(mat) <- list(nfval,c("mean","sd"))
  return(list(fval=fval, mat=mat))
}
print(param(NS))
print(param(NS), show.details = "minimal")
print(param(NS), show.details = "medium")
## Maximum likelihood estimator
res <- MLEstimator(x = x, ParamFamily = NS)
print(res) #equivalent to 'show(res)' or 'res'
print(res, digits = 4)
print(res, show.details = "minimal")
print(res, show.details = "medium")
distrModoptions("show.details" = show.old)



cleanEx()
nameEx("qqplot")
### * qqplot

flush(stderr()); flush(stdout())

### Name: qqplot
### Title: Methods for Function qqplot in Package 'distrMod'
### Aliases: qqplot qqplot-methods qqplot,ANY,ProbFamily-method
###   qqplot,ANY,UnivariateDistribution-method
### Keywords: hplot distribution

### ** Examples

qqplot(r(Norm(15,sqrt(30)))(40), Chisq(df=15))
qqplot(r(Norm(15,sqrt(30)))(40), NormLocationFamily())



cleanEx()
nameEx("symmetricBias-class")
### * symmetricBias-class

flush(stderr()); flush(stdout())

### Name: symmetricBias-class
### Title: symmetric Bias Type
### Aliases: symmetricBias-class
### Keywords: classes

### ** Examples

symmetricBias()
## The function is currently defined as
function(){ new("symmetricBias", name = "symmetric Bias") }



cleanEx()
nameEx("symmetricBias")
### * symmetricBias

flush(stderr()); flush(stdout())

### Name: symmetricBias
### Title: Generating function for symmetricBias-class
### Aliases: symmetricBias
### Keywords: robust

### ** Examples

symmetricBias()

## The function is currently defined as
function(){ new("symmetricBias", name = "symmetric Bias") }



cleanEx()
nameEx("trAsCov-class")
### * trAsCov-class

flush(stderr()); flush(stdout())

### Name: trAsCov-class
### Title: Trace of asymptotic covariance
### Aliases: trAsCov-class
### Keywords: classes

### ** Examples

new("trAsCov")



cleanEx()
nameEx("trAsCov")
### * trAsCov

flush(stderr()); flush(stdout())

### Name: trAsCov
### Title: Generating function for trAsCov-class
### Aliases: trAsCov
### Keywords: robust

### ** Examples

trAsCov()

## The function is currently defined as
function(){ new("trAsCov") }



cleanEx()
nameEx("trFiCov-class")
### * trFiCov-class

flush(stderr()); flush(stdout())

### Name: trFiCov-class
### Title: Trace of finite-sample covariance
### Aliases: trFiCov-class
### Keywords: classes

### ** Examples

new("trFiCov")



cleanEx()
nameEx("trFiCov")
### * trFiCov

flush(stderr()); flush(stdout())

### Name: trFiCov
### Title: Generating function for trFiCov-class
### Aliases: trFiCov
### Keywords: robust

### ** Examples

trFiCov()

## The function is currently defined as
function(){ new("trFiCov") }



cleanEx()
nameEx("trafo-methods")
### * trafo-methods

flush(stderr()); flush(stdout())

### Name: trafo-methods
### Title: Methods for function trafo in Package 'distrMod'
### Aliases: trafo-methods trafo trafo,Estimate,missing-method
###   trafo,Estimate,ParamFamParameter-method
###   trafo,ParamFamParameter,missing-method
###   trafo,ParamWithScaleAndShapeFamParameter,missing-method
###   trafo,ParamFamily,missing-method
###   trafo,ParamFamily,ParamFamParameter-method trafo.fct
###   trafo.fct-methods trafo.fct,ParamFamily-method trafo<-
###   trafo<-,ParamFamParameter-method trafo<-,ParamFamily-method
### Keywords: models

### ** Examples

## Gaussian location and scale
NS <- NormLocationScaleFamily(mean=2, sd=3)
## generate data out of this situation
x <- r(distribution(NS))(30)

## want to estimate mu/sigma, sigma^2
## -> new trafo slot:
trafo(NS) <- function(param){
  mu <- param["mean"]
  sd <- param["sd"]
  fval <- c(mu/sd, sd^2)
  nfval <- c("mu/sig", "sig^2")
  names(fval) <- nfval
  mat <- matrix(c(1/sd,0,-mu/sd^2,2*sd),2,2)
  dimnames(mat) <- list(nfval,c("mean","sd"))
  return(list(fval=fval, mat=mat))
}

## Maximum likelihood estimator
(res <- MLEstimator(x = x, ParamFamily = NS))
## confidence interval
 confint(res)




cleanEx()
nameEx("trafoEst")
### * trafoEst

flush(stderr()); flush(stdout())

### Name: trafoEst
### Title: Function trafoEst in Package 'distrMod'
### Aliases: trafoEst
### Keywords: models

### ** Examples

## Gaussian location and scale
NS <- NormLocationScaleFamily(mean=2, sd=3)
## generate data out of this situation
x <- r(distribution(NS))(30)

## want to estimate mu/sigma, sigma^2
## -> without new trafo slot:
mtrafo <- function(param){
  mu <- param["mean"]
  sd <- param["sd"]
  fval <- c(mu/sd, sd^2)
  nfval <- c("mu/sig", "sig^2")
  names(fval) <- nfval
  mat <- matrix(c(1/sd,0,-mu/sd^2,2*sd),2,2)
  dimnames(mat) <- list(nfval,c("mean","sd"))
  return(list(fval=fval, mat=mat))
}

## Maximum likelihood estimator in the original problem
res0 <- MLEstimator(x = x, ParamFamily = NS)
## transformation
res <- trafoEst(mtrafo, res0)
## confidence interval
 confint(res)



cleanEx()
nameEx("validParameter-methods")
### * validParameter-methods

flush(stderr()); flush(stdout())

### Name: validParameter-methods
### Title: Methods for function validParameter in Package 'distrMod'
### Aliases: validParameter-methods validParameter
###   validParameter,ParamFamily-method validParameter,L2ScaleUnion-method
###   validParameter,L2ScaleFamily-method
###   validParameter,L2LocationFamily-method
###   validParameter,L2LocationScaleFamily-method
###   validParameter,BinomFamily-method validParameter,PoisFamily-method
###   validParameter,L2ScaleShapeUnion-method
### Keywords: models

### ** Examples

 NS <- NormLocationScaleFamily()
 validParameter(NS, c(scale=0.1, loc=2))
 validParameter(NS, c(scale=-0.1, loc=2))
 validParameter(NS, c(scale=0, loc=2))
 validParameter(NS, c(mean=2, sd=2))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
