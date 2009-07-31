### some further examples:
require(distrMod)
options("newDevice"=TRUE)


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
