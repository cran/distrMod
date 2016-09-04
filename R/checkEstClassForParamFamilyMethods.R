setMethod(".checkEstClassForParamFamily",
              signature=signature(PFam="ANY",estimator="ANY"),
              function(PFam, estimator) estimator)
