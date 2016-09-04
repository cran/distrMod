################################################################
# return level - Plot functions in package distrMod
################################################################


setMethod("returnlevelplot", signature(x = "ANY",
                              y = "UnivariateDistribution"),
    function(x,    ### observations
             y,    ### distribution
             n = length(x), ### number of points to be plotted
             withIdLine = TRUE, ### shall line y=x be plotted in
             withConf = TRUE,   ### shall confidence lines be plotted
             withConf.pw  = withConf,   ### shall pointwise confidence lines be plotted
             withConf.sim = withConf,   ### shall simultaneous confidence lines be plotted
             plot.it = TRUE,    ### shall be plotted at all (inherited from stats::qqplot)
             datax = FALSE,     ### as in qqnorm
             MaxOrPOT = c("Max","POT"), ### used for block maxima or points over threshold
             npy = 365, ### number of observations per year
             threshold = if(is(y,"GPareto")) NA else 0,
             xlab = deparse(substitute(x)), ## x-label
             ylab = deparse(substitute(y)), ## y-label
             main = "",
             ...,                 ## further parameters
             width = 10,          ## width (in inches) of the graphics device opened
             height = 5.5,        ## height (in inches) of the graphics device opened}
             withSweave = getdistrOption("withSweave"), ## logical: if \code{TRUE}
             ##               (for working with \command{Sweave}) no extra device is opened and height/width are not set
             mfColRow = TRUE,     ## shall we use panel partition mfrow=c(1,1)?
             n.CI = n,            ## number of points to be used for CI
             withLab = FALSE,     ## shall observation labels be plotted in
             lab.pts = NULL,      ## observation labels to be used
             which.lbs = NULL,    ## which observations shall be labelled
             which.Order = NULL,  ## which of the ordered (remaining) observations shall be labelled
             order.traf = NULL,   ## an optional trafo; by which the observations are ordered (as order(trafo(obs))
             col.IdL = "red",     ## color for the identity line
             lty.IdL = 2,         ## line type for the identity line
             lwd.IdL = 2,         ## line width for the identity line
             alpha.CI = .95,      ## confidence level
             exact.pCI = (n<100), ## shall pointwise CIs be determined with exact Binomial distribution?
             exact.sCI = (n<100), ## shall simultaneous CIs be determined with exact kolmogorov distribution?
             nosym.pCI = FALSE,   ## shall we use (shortest) asymmetric CIs?
             col.pCI = "orange",  ## color for the pointwise CI
             lty.pCI = 3,         ## line type for the pointwise CI
             lwd.pCI = 2,         ## line width for the pointwise CI
             pch.pCI = par("pch"),## symbol for points (for discrete mass points) in pointwise CI
             cex.pCI = par("cex"),## magnification factor for points (for discrete mass points) in pointwise CI
             col.sCI = "tomato2", ## color for the simultaneous CI
             lty.sCI = 4,         ## line type for the simultaneous CI
             lwd.sCI = 2,         ## line width for the simultaneous CI
             pch.sCI = par("pch"),## symbol for points (for discrete mass points) in simultaneous CI
             cex.sCI = par("cex"),## magnification factor for points (for discrete mass points) in simultaneous CI
             added.points.CI = TRUE, ## should the CIs be drawn through additional points?
             cex.pch = par("cex"),## magnification factor for the plotted symbols
             col.pch = par("col"),## color for the plotted symbols
             cex.lbl = par("cex"),## magnification factor for the plotted observation labels
             col.lbl = par("col"),## color for the plotted observation labels
             adj.lbl = NULL,      ## adj parameter for the plotted observation labels
             alpha.trsp = NA,     ## alpha transparency to be added afterwards
             jit.fac = 0,         ## jittering factor used for discrete distributions
             jit.tol = .Machine$double.eps, ## tolerance for jittering: if distance 
                                 #is smaller than jit.tol, points are considered replicates
             check.NotInSupport = TRUE, ## shall we check if all x lie in support(y)
             col.NotInSupport = "red", ## if preceding check TRUE color of x if not in support(y)
             with.legend = TRUE,  ## shall a legend be plotted
             legend.bg = "white", ## background for the legend
             legend.pos = "topleft", ## position for the legend
             legend.cex = 0.8,     ## magnification factor for the legend
             legend.pref = "",     ## prefix for legend  text
             legend.postf = "",    ## postfix for legend text
             legend.alpha = alpha.CI, ## nominal level of CI
             debug = FALSE, ## shall additional debug output be printed out?
             withSubst = TRUE
    ){ ## return value as in stats::qqplot

    MaxOrPOT <- match.arg(MaxOrPOT)
    mc <- match.call(call = sys.call(sys.parent(1)))
    xcc <- as.character(deparse(mc$x))

   .mpresubs <- if(withSubst){
                   function(inx) 
                    .presubs(inx, c("%C", "%A", "%D" ),
                          c(as.character(class(x)[1]), 
                            as.character(date()), 
                            xcc))
               }else function(inx)inx

    if(missing(xlab)) mc$xlab <- paste(gettext("Return level of"),
                                       as.character(deparse(mc$x)))
    if(missing(ylab)) mc$ylab <- gettext("Return period (years)")
    if(missing(main)) mc$main <- gettext("Return level plot")
    mcl <- as.list(mc)[-1]
    mcl$datax <- NULL
    mcl$MaxOrPOT <- NULL
    mcl$npy <- NULL
    mcl$withSweave <- NULL
    mcl$mfColRow <- NULL
    mcl$type <-NULL
    mcl$debug <- NULL
    mcl$added.points.CI <- NULL
    if(is.null(mcl$datax)) datax <- FALSE
    force(x)

    thresh0 <- threshold 
    if(is(y,"GPareto")){ 
       if(is.na(threshold)) thresh0 <- location(y)
       y <- y - thresh0
       x <- x + thresh0
    }              

    xj <- sort(x)
    if(any(.isReplicated(x, jit.tol))&&jit.fac>0)
       xj[.isReplicated(x, jit.tol)] <- jitter(x[.isReplicated(x, jit.tol)], factor=jit.fac)

    xj <- sort(xj)
    ord.x <- order(xj)

    p2rl <- function(pp){
               pp <- p(y)(pp)
               return(if(MaxOrPOT=="Max") -1/log(pp) else  1/(1-pp)/npy)
    }

    pp <- ppoints(length(xj))
    yc.o <- q(y)(pp)
    ycl <- p2rl(yc.o)

    ### extend range somewhat
#    pyn <- p(y)(10^(seq(-1, 3.75 + log10(npy), by = 0.1)))
    xyall <- force(sort(unique(c(yc.o,x,
                    q(y)(c(seq(0.01, 0.09, by = 0.01),(1:9)/10,
                         0.95, 0.99, 0.995, 0.999))
                         ))))
    rxyall  <- (max(xyall)-min(xyall))*0.6
    rxymean <- (max(xyall)+min(xyall))/2

    xyallc  <- seq(from=rxymean-rxyall,to=rxymean+rxyall, length.out=400)
#    print(xyallc)
    pxyall  <- p(y)(xyallc)
#    print(pxyall)

    pxyallc <- p2rl(xyallc)
     xyallc <-  xyallc[pxyall>0.00001 & pxyall<0.99999]
    pxyallc <- pxyallc[pxyall>0.00001 & pxyall<0.99999]

#    print(cbind(pxyallc,xyallc))

    if("support" %in% names(getSlots(class(y))))
       ycl <- sort(jitter(ycl, factor=jit.fac))

    alp.v <- .makeLenAndOrder(alpha.trsp,ord.x)
    alp.t <- function(x,a1) if(is.na(x)) x else addAlphTrsp2col(x,a1)
    alp.f <- if(length(alpha.trsp)==1L && is.na(alpha.trsp))
             function(x,a) x else function(x,a) mapply(x,alp.t,a1=a)
    cex.pch <- .makeLenAndOrder(cex.pch,ord.x)
    cex.lbl <- .makeLenAndOrder(cex.lbl,ord.x)
    col.pch <- alp.f(.makeLenAndOrder(col.pch,ord.x),alp.v)
    col.lbl <- alp.f(.makeLenAndOrder(col.lbl,ord.x),alp.v)

    if(withLab){
      if(is.null(lab.pts)) lab.pts <- paste(ord.x)
      else lab.pts <- .makeLenAndOrder(lab.pts,ord.x)
    }

    if(check.NotInSupport){
       xo <- x[ord.x]
       nInSupp <- which(xo < q(y)(0))

       nInSupp <- unique(sort(c(nInSupp,which( xo > q(y)(1)))))
       if("support" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( ! xo %in% support(y)))))
       if("gaps" %in% names(getSlots(class(y))))
          nInSupp <- unique(sort(c(nInSupp,which( .inGaps(xo,gaps(y))))))
       if(length(nInSupp)){
          col.pch[nInSupp] <- col.NotInSupport
          if(withLab)
#             col.lbl[ord.x[nInSupp]] <- col.NotInSupport
             col.lbl[nInSupp] <- col.NotInSupport
       }
    }


    if(n!=length(x)) withLab <- FALSE

    mcl <- .deleteItemsMCL(mcl)
    mcl$cex <- cex.pch
    mcl$col <- col.pch

    if (!withSweave){
           devNew(width = width, height = height)
    }
    opar <- par("mfrow", no.readonly = TRUE)
    if(mfColRow) on.exit(do.call(par, list(mfrow=opar, no.readonly = TRUE)))

    if(mfColRow) opar1 <- par(mfrow = c(1,1), no.readonly = TRUE)

    ret <- list(x=xj,y=ycl)

    if(plot.it){
       xallc1 <- sort(c(xj,xyallc))
       yallc1 <- sort(c(ycl,pxyallc))
       mcl$x <- mcl$y <- NULL
       if(datax){
          mcl$xlab <- xlab
          mcl$ylab <- ylab
          do.call(plot, c(list(x=xallc1, y=yallc1, log="y",type="n"),mcl))
          do.call(points, c(list(x=xj, y=ycl), mcl))
    #       ret <- do.call(stats::qqplot, args=mcl0, log="y", ylim = c(0.1,1000))
       }else{
          mcl$ylab <- xlab
          mcl$xlab <- ylab
          do.call(plot, c(list(x=yallc1, y=xallc1, log="x",type="n"),mcl))
          do.call(points, c(list(x=ycl, y=xj),mcl))
       }
    }

    if(withLab&& plot.it){
       lbprep <- .labelprep(xj,yc.o,lab.pts,
                            col.lbl,cex.lbl,which.lbs,which.Order,order.traf)
       lbprep$y0 <- p2rl(lbprep$y0)
       xlb0 <- if(datax) lbprep$x0 else lbprep$y0
       ylb0 <- if(datax) lbprep$y0 else lbprep$x0
       text(x = xlb0, y = ylb0, labels = lbprep$lab,
            cex = lbprep$cex, col = lbprep$col, adj = adj.lbl)
    }

    if(withIdLine){
       if(plot.it){
          if(datax){
             lines(xyallc,pxyallc,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
          }else{
             lines(pxyallc,xyallc,col=col.IdL,lty=lty.IdL,lwd=lwd.IdL)
          }
       }
       qqb <- NULL
       if(#is(y,"AbscontDistribution")&&
       withConf){

          if(added.points.CI){
             xy <- unique(sort(c(x,xj,xyallc,yc.o)))
          }else{
             xy <- unique(sort(c(x,xj,yc.o)))
          }
          xy <- xy[!.NotInSupport(xy,y)]
          lxy <- length(xy)
          if(is(y,"DiscreteDistribution")){
             n0 <- min(n.CI, length(support(y)))
             n1 <- max(n0-lxy,0)
             if (n1 >0 ){
                 notyetInXY <- setdiff(support(y), xy)
                 xy0 <- sample(notyetInXY, n1)
                 xy <- sort(unique(c(xy,xy0)))
             }
          }else{
             if(lxy < n.CI){
                n1 <- (n.CI-lxy)%/%3
                xy0 <- seq(min(xy),max(xy),length=n1)
                xy1 <- r(y)(n.CI-lxy-n1)
                xy <- sort(unique(c(xy,xy0,xy1)))
             }
          }

        #qqb <- qqbounds(sort(unique(xy)),y,alpha.CI,n,withConf.pw, withConf.sim,
        #                   exact.sCI,exact.pCI,nosym.pCI, debug = debug)
        #qqb$crit <- p2rl(qqb$crit)
        if(plot.it){
          qqb <- .confqq(xy, y, datax, withConf.pw, withConf.sim, alpha.CI,
                      col.pCI, lty.pCI, lwd.pCI, pch.pCI, cex.pCI,
                      col.sCI, lty.sCI, lwd.sCI, pch.sCI, cex.sCI,
                  n, exact.sCI = exact.sCI, exact.pCI = exact.pCI,
                  nosym.pCI = nosym.pCI, with.legend = with.legend,
                  legend.bg = legend.bg, legend.pos = legend.pos,
                  legend.cex = legend.cex, legend.pref = legend.pref,
                  legend.postf = legend.postf, legend.alpha = legend.alpha,
                  qqb0=NULL, transf0=p2rl, debug = debug)
       }
    }}
    return(invisible(c(ret,qqb)))
    })

## into distrMod
setMethod("returnlevelplot", signature(x = "ANY",
                              y = "ProbFamily"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- paste(gettext("Return Level of"), as.character(deparse(mc$x)))
    if(missing(ylab)) mc$ylab <- paste(gettext("Return Period at"), as.character(deparse(mc$y)))
    mcl <- as.list(mc)[-1]

    mcl$y <- yD <- y@distribution
    if(!is(yD,"UnivariateDistribution"))
       stop("Not yet implemented.")

    return(invisible(do.call(getMethod("returnlevelplot", signature(x="ANY", y="UnivariateDistribution")),
            args=mcl)))
    })

setMethod("returnlevelplot", signature(x = "ANY",
                              y = "Estimate"), function(x, y,
                              n = length(x), withIdLine = TRUE, withConf = TRUE,
    withConf.pw  = withConf,  withConf.sim = withConf,
    plot.it = TRUE, xlab = deparse(substitute(x)),
    ylab = deparse(substitute(y)), ...){

    mc <- match.call(call = sys.call(sys.parent(1)))
    if(missing(xlab)) mc$xlab <- paste(gettext("Return Level of"), as.character(deparse(mc$x)))
    mcl <- as.list(mc)[-1]

    param <- ParamFamParameter(main=untransformed.estimate(y), nuisance=nuisance(y),
                               fixed=fixed(y))

    es.call <- y@estimate.call
    nm.call <- names(es.call)
    PFam <- NULL
    if("ParamFamily" %in% nm.call)
       PFam <- eval(as.list(es.call)[["ParamFamily"]])
    if(is.null(PFam))
       stop("There is no object of class 'ProbFamily' in the call of 'x'")

    PFam0 <- modifyModel(PFam, param)
    mcl$y <- PFam0
    if(missing(ylab)) mc$ylab <- paste(gettext("Return Period at fitted"), name(PFam0))

    return(invisible(do.call(getMethod("returnlevelplot", signature(x="ANY", y="ProbFamily")),
            args=mcl)))
    })
