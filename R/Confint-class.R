setMethod("confint", signature(object="Confint", method="missing"),
           function(object, method) object@confint)

setMethod("call.estimate", signature(object="Confint"),
           function(object) object@call.estimate)
setMethod("name.estimate", signature(object="Confint"),
           function(object) object@name.estimate)
setMethod("samplesize.estimate", signature(object="Confint"),
           function(object) object@samplesize.estimate)
setMethod("nuisance.estimate", signature(object="Confint"),
           function(object) object@nuisance.estimate)
setMethod("trafo.estimate", signature(object="Confint"),
           function(object) object@trafo.estimate)
setMethod("fixed.estimate", signature(object="Confint"),
           function(object) object@fixed.estimate)
setMethod("type", signature(object="Confint"),
           function(object) object@type)

