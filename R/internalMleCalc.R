#################### internal helper functions in connection with
#################### m[l,c]eCalc

### not exported:
.negLoglikelihood <- function(x, Distribution, ...){
          ### increase accuracy:
           if(Distribution@.withSim||!.inArgs("log",d(Distribution)))
               res <- -sum(log(Distribution@d(x, ...)))
           else
               res <- -sum(Distribution@d(x, log = TRUE, ...))
           return(res)
    }


#

##########################################################################
#internal helper
##########################################################################
.process.meCalcRes <- function(res, PFam, trafo, res.name, call,
                               asvar.fct, ...){
    lmx <- length(main(PFam))
    lnx <- length(nuisance(PFam))
    idx <- 1:lmx
    jdx <-      if(lnx) lmx + 1:lnx else idx
    nuis.idx <- if(lnx) jdx else NULL

    theta <- res$estimate
    crit <- res$criterion
    param <- res$param
    
    crit.name <- ""
    if(res$crit.name == ""){
       if(!is.null(names(res$criterion))) 
           if(names(res$criterion) != "") 
              crit.name <- names(res$criterion)
            
    }else crit.name <- res$crit.name
    names(crit) <- crit.name
    est.name <-  if(crit.name=="") "Minimum criterion estimate"  else
                    paste("Minimum", crit.name, "estimate", sep = " ") 

    if(is.null(res$Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                        dimnames=list(character(0), c("method", "message")))
    else{
        Infos <- matrix(c(rep("MCEstimator", length(Infos)), Infos), ncol = 2)
        colnames(Infos) <- c("method", "message")
    }



    traf1 <- PFam@param@trafo
    if(is.null(trafo))
         {if(is.matrix(traf1))
             traf0 <- list(fct = function(x)
                                 list(fval = traf1 %*% x, mat = traf1),
                           mat = traf1)
          else
             traf0 <- list(fct = traf1, mat = traf1(main(param))$mat)
         }
    else {if(is.matrix(trafo))
             traf0 <- list(fct = function(x)
                                 list(fval = trafo %*% x, mat = trafo),
                           mat = trafo)
          else
             traf0 <- list(fct = trafo, mat = trafo(main(param))$mat)
         }

    
    if(!validParameter(PFam,param))
          {warning("Optimization for MCE did not give a valid result. You could try to use argument 'penalty'.")
           theta <- as.numeric(rep(NA, lnx+lmx))
           res <- new("MCEstimate", name = est.name, estimate = theta,
                       criterion = crit, Infos = Infos, samplesize = samplesize,
                       nuis.idx = nuis.idx, estimate.call = call,
                       trafo = traf0)
           return(res)}
    
    estimate <- theta[idx]

    asvar <- NULL
    if(!missing(asvar.fct))
       asvar <- asvar.fct(PFam, param, ...)

    untransformed.estimate <- theta
    untransformed.asvar <- asvar

    if(!.isUnitMatrix(traf0$mat)){
       estimate <- traf0$fct(estimate)$fval
       trafm <- traf0$mat
       if(!is.null(asvar)){
           asvar <- trafm%*%asvar[idx,idx]%*%t(trafm)
           rownames(asvar) <- colnames(asvar) <- c(names(estimate))
          }
    }

    new("MCEstimate", name = est.name, estimate = estimate, criterion = crit,
         asvar = asvar, Infos = Infos, samplesize = res$samplesize,
         nuis.idx = nuis.idx, estimate.call = call, trafo = traf0,
         untransformed.estimate = untransformed.estimate,
         untransformed.asvar = untransformed.asvar,
         criterion.fct = res$crit.fct, method = res$method,
         fixed = fixed(param))
}

##########################################################################

## caching to speed up things:
.inArgs <- distr:::.inArgs

