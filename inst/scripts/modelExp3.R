##########################################################
### a new central distribution  in a location scale model
##########################################################

require(distrMod)

my3d <- AbscontDistribution( d = function(x) exp(-abs(x)^3), withS = TRUE)
plot(my3d)
## in a location scale model
scl.true <- 2; loc.true <- 3
my3dF <- L2LocationScaleFamily(loc = loc.true,
                               scale = scl.true,
                               name = "my3dF",
                               centraldistribution = my3d)

plot(my3dF)

### generate some data out of the model
x <- r(my3dF)(40)

### evaluate the MLE:
mledistrMod <- MLEstimator(x,my3dF)

#some profiling
par(mfrow=c(1,2))
plot(profile(mledistrMod))

# a confidence interval
confint(mledistrMod)

(mde.kolm <- (x = x, ParamFamily = my3dF))
(mde.CvM <- MDEstimator(x = x, ParamFamily = my3dF, distance = CvMDist))
asvar(mde.CvM) <- distrMod:::.CvMMDCovariance(my3dF, 
                  param = ParamFamParameter(main= estimate(MDE)),
                  expon = 2, withplot = TRUE)
# a confidence interval
confint(mde.CvM)
