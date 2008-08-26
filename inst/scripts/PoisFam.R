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


