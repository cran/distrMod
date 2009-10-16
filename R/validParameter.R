 ################## validParameter  Methods #############################
 
 setMethod("validParameter", signature(object = "ParamFamily"),
          function(object, param){
             if(is(param,"ParamFamParameter") && length(nuisance(object)))
                  theta <- c(main(param), nuisance(param))
             else {if (is(param,"ParamFamParameter"))  param <- main(param) 
                   theta <- param
                   }
             if(!all(is.finite(theta))) return(FALSE)
             if(length(param)<0 || length(theta) > length(param(object)))
                return(FALSE)
             if(is( try(dum <- object@modifyParam(theta), silent = TRUE),
                "try-error"))  return(FALSE)
             return(TRUE)})

 setMethod("validParameter", signature(object = "L2ScaleFamily"),
          function(object, param, tol=.Machine$double.eps){
             if(is(param,"ParamFamParameter"))
                param <- main(param)
             if(!all(is.finite(param))) return(FALSE)
             if(length(param)!=1) return(FALSE)
             return(param > tol)})

 setMethod("validParameter", signature(object = "L2LocationFamily"),
          function(object, param){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)!=1) return(FALSE)
          TRUE})

 setMethod("validParameter", signature(object = "L2LocationScaleFamily"),
          function(object, param, tol=.Machine$double.eps){
             if(is(param,"ParamFamParameter") && length(nuisance(object)))
                  theta <- c(main(param), nuisance(param))
             else {if (is(param,"ParamFamParameter"))  param <- main(param) 
                   theta <- param
                   }
          if(!all(is.finite(theta))) return(FALSE)
          if(length(theta)>2||length(theta)<1) return(FALSE)
          if(any(names(theta)%in% c("scale", "sd"))){
              return(theta[names(theta)%in% c("scale", "sd")]>tol)
              }
          return(TRUE)
          })


 setMethod("validParameter", signature(object = "BinomFamily"),
          function(object, param, tol=.Machine$double.eps){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)!=1) return(FALSE)
          if(param<= tol || param>= 1-tol)
             return(FALSE)
          return(TRUE)
          })
 setMethod("validParameter", signature(object = "PoisFamily"),
          function(object, param, tol=.Machine$double.eps){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)!=1) return(FALSE)
          if(param<= tol) return(FALSE)
          return(TRUE)
          })

 setMethod("validParameter", signature(object = "GammaFamily"),
          function(object, param, tol=.Machine$double.eps){
          if(is(param,"ParamFamParameter"))
                param <- main(param)
          if(!all(is.finite(param))) return(FALSE)
          if(length(param)>2||length(param)<1) return(FALSE)
          if(any(param<= tol)) return(FALSE)
          return(TRUE)
          })

