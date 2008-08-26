### some further examples:
require(distrMod)

### Poisson Family
P <- PoisFamily(3)
# generate data
x <- r(P)(40)
MLEstimator(x,P)
MDEstimator(x,P)
MDEstimator(x,P, distance = CvMDist, asvar.fct = distrMod:::.CvMMDCovariance)
MDEstimator(x,P, distance = CvMDist, mu = Norm())
MDEstimator(x,P, distance = TotalVarDist)


### Beta Family
B <- BetaFamily(2,4)
# generate data
x <- r(B)(40)
distroptions(DistrResolution = 1e-10)
MDEstimator(x, B, distance = TotalVarDist)
MDEstimator(x, B)
MDEstimator(x, B, distance = CvMDist, asvar.fct = distrMod:::.CvMMDCovariance)
(MLE<-MLEstimator(x, B))
confint(MLE)

### a new central distribution
my3d <- AbscontDistribution( d = function(x) exp(-abs(x)^3), withS = TRUE)
plot(my3d)
my3dF <- L2LocationScaleFamily(name = "my3dF",
                               centraldistribution = my3d)

plot(my3dF)

### generate some data out of the model
x <- r(my3dF)(40)*3+2
### evaluate the MLE:
MLEstimator(x,my3dF)
MDE = MDEstimator(x = x, ParamFamily = my3dF, distance = CvMDist)

MDE.asvar <- distrMod:::.CvMMDCovariance(my3dF, 
                 param = ParamFamParameter(main= estimate(MDE)),
                 expon = 2, withplot = TRUE)

