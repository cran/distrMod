require(distrMod)

## example to CvM MDE for Normal Location and Scale

x=rnorm(30)
NF=NormLocationScaleFamily()

system.time(print(MDEstimator(x,NF,CvMDist)))
#with useApply
system.time(print(MDEstimator(x,NF,CvMDist,useApply=TRUE)))

MDEstimator(rnorm(30),NF,CvMDist)
#another sample
MDEstimator(rnorm(30),NF,CvMDist)
# larger sample size
MDEstimator(rnorm(300),NF,CvMDist)

MDEstimator(rnorm(300,mean=2,sd=2),NF,CvMDist)
#another sample
MDEstimator(rnorm(300,mean=2,sd=2),NF,CvMDist)
