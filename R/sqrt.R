setMethod("sqrt", signature(x = "PosSemDefSymmMatrix"), function(x){
            er <- eigen(x)
            d <- sqrt(er$values)
            return(er$vectors %*% diag(d) %*% t(er$vectors))
})

