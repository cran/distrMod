\name{internals_for_qqplot}
\alias{internals_for_qqplot}
\alias{.labelprep}


\title{Internal functions for qqplot of package distrMod}

\description{
These functions are used internally by qqplot of package distrMod.}

\usage{
.labelprep(x,y,lab.pts,col.lbl,cex.lbl,adj.lbl,which.lbs,which.Order,order.traf)
}


\arguments{
\item{x}{a (numeric) vector}
\item{y}{a (numeric) vector of same length as \code{x}}
\item{cex.lbl}{magnification factor for the plotted observation labels}
\item{col.lbl}{color for the plotted observation labels}
\item{adj.lbl}{adj parameter for the plotted observation labels}
\item{lab.pts}{character or \code{NULL}; observation labels to be used}
\item{which.lbs}{integer or \code{NULL}; which observations shall be labelled}
\item{which.Order}{integer or \code{NULL}; which of the ordered (remaining) observations shall be labelled}
\item{order.traf}{function or \code{NULL}; an optional trafo by which the observations are ordered (as order(trafo(obs)).}
}

\details{
\code{.labelprep} prepares the labels for later output:
first some indices among \code{x} according
to \code{which.lbs} are (optionally) selected, 
then orders both \code{x} and \code{y} according to the
order of the remaining (optionally transformed by \code{order.traf}) values of
\code{x} (decreasing); finally only the order statistics (w.r.t. the preceding
ordering) as specified in \code{which.Order} are returned;
correspondingly, labels, colors, and magnification factors, are ordered
respectively thinned out.
}


\value{
\item{.labelprep}{a list with components
\code{x0} (the thinned out and ordered vector \code{x}),
\code{y0} (the thinned out and ordered vector \code{y}),
\code{lab} (the thinned out and ordered vector of labels \code{lab.pts}),
\code{col} (the thinned out and ordered vector of colors \code{col.lbs}),
\code{cex} (the thinned out and ordered vector of magnification factors \code{cex.lbs}).
\code{adj} (the thinned out and ordered vector of adjacencies \code{adj.lbs}).
}
}

\author{
  Peter Ruckdeschel \email{peter.ruckdeschel@uni-oldenburg.de},
}

\seealso{\code{\link[stats]{ks.test}}, \code{\link[distr]{qqplot}}
,\code{\link[distr:internals-qqplot]{internals_qqplot}}
, \code{\link[distrMod]{qqplot}}, \code{\link[RobAStBase]{qqplot}}
}
\keyword{internal}
\concept{utilities}
