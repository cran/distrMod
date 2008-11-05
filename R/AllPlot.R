### from Matthias' thesis / ROptEst
setMethod("plot", signature(x = "ParamFamily", y = "missing"),
    function(x, ...){ 
        e1 <- x@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")

        plot(e1) 
    })
setMethod("plot", signature(x = "L2ParamFamily", y = "missing"),
    function(x, withSweave = getdistrOption("withSweave"), 
             main = FALSE, inner = TRUE, sub = FALSE, 
             col.inner = par("col.main"), cex.inner = 0.8, 
             bmar = par("mar")[1], tmar = par("mar")[3], ...,
             mfColRow = TRUE){

        xc <- match.call(call = sys.call(sys.parent(1)))$x
        dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
        
        if(!is.logical(inner)){
           if(!is.list(inner)||length(inner) != 4)
               stop("Argument 'inner' must either be 'logical' or a 'list' vector of length 4")
           innerD <- inner[1:3]
           innerL <- inner[4] 
        }else{innerD <- innerL <- inner}
        
        if(!is.null(dots[["lty"]]))  dots["lty"] <- NULL
        if(!is.null(dots[["type"]])) dots["type"] <- NULL
        if(!is.null(dots[["xlab"]])) dots["xlab"] <- NULL
        if(!is.null(dots[["ylab"]])) dots["ylab"] <- NULL

        e1 <- x@distribution
        if(!is(e1, "UnivariateDistribution")) stop("not yet implemented")


        if(is(e1, "AbscontDistribution")){
            lower <- ifelse(is.finite(q(e1)(0)), q(e1)(0), q(e1)(getdistrOption("TruncQuantile")))
            upper <- ifelse(is.finite(q(e1)(1)), q(e1)(1), q(e1)(1 - getdistrOption("TruncQuantile")))
            h <- upper - lower
            x.vec <- seq(from = lower - 0.1*h, to = upper + 0.1*h, length = 1000)
            plty <- "l"
            lty <- "solid"
        }else{
            if(is(e1, "DiscreteDistribution")){
                x.vec <- support(e1)
                plty <- "p"
                lty <- "dotted"
            }else{
                x.vec <- r(e1)(1000)
                x.vec <- sort(unique(x.vec))
                plty <- "p"
                lty <- "dotted"
            }
        }

        dims <- length(x@param)
        L2deriv <- as(diag(dims) %*% x@L2deriv, "EuclRandVariable")

        mainL <- FALSE
        subL <- FALSE
        lineT <- NA

     .mpresubs <- function(inx)
                    distr:::.presubs(inx, c("%C", "%D", "%A"),
                          c(as.character(class(x)[1]),
                            as.character(date()),
                            as.character(deparse(xc))))

     if (hasArg(main)){
         mainL <- TRUE
         if (is.logical(main)){
             if (!main) mainL <-  FALSE
             else
                  main <- gettextf("Distribution Plot for model %%A") ###
                          ### double  %% as % is special for gettextf
             }
         main <- .mpresubs(main)
         if (mainL) {
             if(missing(tmar))
                tmar <- 5
             if(missing(cex.inner))
                cex.inner <- .65
             lineT <- 0.6
             }
     }
     if (hasArg(sub)){
         subL <- TRUE
         if (is.logical(sub)){
             if (!sub) subL <-  FALSE
             else       sub <- gettextf("generated %%D")
                          ### double  %% as % is special for gettextf
         }
         sub <- .mpresubs(sub)
         if (subL)
             if (missing(bmar)) bmar <- 6
     }

     if(is.logical(innerL)){
        innerT <- paste(gettextf("Component "), 1:dims,
                        gettextf(" of L_2 derivative\nof"),
                        name(x)[1],
                        gettextf("\nwith main parameter ("),
                        paste(round(x@param@main, 3), collapse = ", "),")")
        if(!is.null(x@param@nuisance))
            innerT <- paste(innerT,
                        gettextf("\nand nuisance parameter ("),
                        paste(round(x@param@nuisance, 3), collapse = ", "),
                        ")",
                        sep=""  )
        if(!is.null(x@param@fixed))
            innerT <- paste(innerT,
                        gettextf("\nand fixed known parameter ("),
                        paste(round(x@param@fixed, 3), collapse = ", "),
                        ")",
                        sep=""  )
     }else{
        innerT <- rep(sapply(inner, .mpresubs), length.out=dims)
     }

        dotsT <- dots
        dotsT["main"] <- NULL
        dotsT["cex.main"] <- NULL
        dotsT["col.main"] <- NULL
        dotsT["line"] <- NULL

     do.call(plot, c(list(e1,withSweave = withSweave, 
             main = main, inner = innerD, sub = sub, 
             col.inner = col.inner, cex.inner = 1.5*cex.inner),
             dots, mfColRow=mfColRow))       
     
      w0 <- options("warn")
        options(warn = -1)
        opar <- par()
        if (!withSweave)
             devNew()
        nrows <- trunc(sqrt(dims))
        ncols <- ceiling(dims/nrows)
        
        if(mfColRow)
           parArgs <- list(mfrow = c(nrows, ncols))

        omar <- par("mar")
        parArgs <- c(parArgs,list(mar = c(bmar,omar[2],tmar,omar[4])))
       
     do.call(par,args=parArgs)
        for(i in 1:dims){
            do.call(plot, args=c(list(x=x.vec, y=sapply(x.vec, L2deriv@Map[[i]]),
                                 type = plty, lty = lty,
                                 xlab = "x",
                                 ylab = expression(paste(L[2], " derivative"))),
                                 dots))
            if(is(e1, "DiscreteDistribution")){
                x.vec1 <- seq(from = min(x.vec), to = max(x.vec), length = 1000)
                do.call(lines, args=c(list(x.vec1, sapply(x.vec1, L2deriv@Map[[i]]),
                              lty = "dotted"),dots))
            }
            do.call(title, args = c(list(main = innerT[i]), dotsT, line = lineT,
                    cex.main = cex.inner, col.main = col.inner))
        }

        if(!hasArg(cex.main)) cex.main <- par("cex.main") else cex.main <- dots$"cex.main"
        if(!hasArg(col.main)) col.main <- par("col.main") else col.main <- dots$"col.main"
        if (mainL)
            mtext(text = main, side = 3, cex = cex.main, adj = .5,
                  outer = TRUE, padj = 1.4, col = col.main)

        if(!hasArg(cex.sub)) cex.sub <- par("cex.sub") else cex.sub <- dots$"cex.sub"
        if(!hasArg(col.sub)) col.sub <- par("col.sub") else col.sub <- dots$"col.sub"
        if (subL)
            mtext(text = sub, side = 1, cex = cex.sub, adj = .5,
                  outer = TRUE, line = -1.6, col = col.sub)

     par(opar)
     options(w0)
     invisible()
    })
