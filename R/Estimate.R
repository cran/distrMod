###############################################################################
## Functions and methods for "Estimate" classes
###############################################################################

setMethod("name", "Estimate", function(object) object@name)
setReplaceMethod("name", "Estimate", 
                  function(object, value) {object@name <- value; object})

setMethod("estimate", "Estimate", function(object) object@estimate)
setMethod("untransformed.estimate", "Estimate", 
           function(object) object@untransformed.estimate)
setMethod("estimate.call", "Estimate", function(object) object@estimate.call)

setMethod("trafo", signature(object = "Estimate", param = "missing"), 
           function(object, param) object@trafo)

setMethod("trafo", signature(object = "Estimate", param = "ParamFamParameter"), 
   function(object, param){
        if(is.function(trafo(param))) 
             return(list(fct = trafo(param), 
                         mat = (trafo(param)(object@untransformed.estimate))$mat))
        else return(list(fct = function(x) trafo(param)%*%x, 
                         mat = trafo(param)))
           
   })

setMethod("fixed", signature(object = "Estimate"), 
           function(object) object@fixed)

setMethod("Infos", "Estimate", function(object) object@Infos)
setReplaceMethod("Infos", "Estimate", 
    function(object, value){ 
        object@Infos <- value 
        if(!is.character(value))
            stop("'value' is no matrix of characters")
        if(ncol(value)!=2)
            stop("'value' has to be a matrix with two columns")
        object
    })

setMethod("addInfo<-", "Estimate", 
    function(object, value){ 
        object@Infos <- rbind(object@Infos, " " = value) 
        if(length(value)!=2)
            stop("length of 'value' is != 2")
        if(!is.character(value))
            stop("'value' is no vector of characters")
        object 
    })

setMethod("samplesize", "Estimate", function(object) object@samplesize)
setMethod("asvar", "Estimate", function(object) object@asvar)

setReplaceMethod("asvar", "Estimate", 
                  function(object, value){ 
          mat <- trafo(object)$mat
          if(.isUnitMatrix(mat)){
             object@asvar <- value
          }else{   
             object@untransformed.asvar <- value
             object@asvar <- mat%*%value%*%t(mat)
          }
          object})

setMethod("untransformed.asvar", "Estimate", function(object) 
           object@untransformed.asvar)

setMethod("criterion", "MCEstimate", function(object) object@criterion)
setMethod("criterion.fct", "MCEstimate", function(object) object@criterion.fct)
setMethod("method", "MCEstimate", function(object) object@method)

setReplaceMethod("criterion", "MCEstimate", 
                  function(object, value) {object@criterion <- value; object})


setMethod("nuisance", "Estimate", function(object) { 
      if(is.null(object@nuis.idx))
         return(NULL)
      else return (object@estimate[object@nuis.idx])    
      })
setMethod("main", "Estimate", function(object) { 
      if(is.null(object@nuis.idx))
         return(object@estimate)
      else return (object@estimate[-object@nuis.idx])    
      })


