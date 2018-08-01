\name{internals_for_qqplot}
\alias{internals_for_qqplot}
\alias{.labelprep}


\title{Internal functions for qqplot of package distrMod}

\description{
These functions are used internally by qqplot of package distrMod.}

\usage{
.labelprep(x,y,lab.pts,col.lbs,cex.lbs,adj.lbs,which.lbs,
           which.Order,order.traf, which.nonlbs)
}


\arguments{
\item{x}{a (numeric) vector}
\item{y}{a (numeric) vector of same length as \code{x}}
\item{cex.lbs}{magnification factor for the plotted observation labels}
\item{col.lbs}{color for the plotted observation labels}
\item{adj.lbs}{adjustment factor for the plotted observation labels}
\item{lab.pts}{character or \code{NULL}; observation labels to be used}
\item{which.lbs}{integer or \code{NULL}; which observations shall be labelled}
\item{which.Order}{integer or \code{NULL}; which of the ordered (remaining) observations shall be labelled}
\item{order.traf}{function or \code{NULL}; an optional trafo by which the observations are ordered (as order(trafo(obs)).}
\item{which.nonlbs}{indices of the observations which should be plotted but
       not labelled; either an integer vector with the indices of the observations
          to be plotted into graph or \code{NULL} --- then all non-labelled
          observations are plotted.}
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
\code{cex} (the thinned out and ordered vector of magnification factors \code{cex.lbs}),
\code{ord} (the (ordered according to y) indices of the observations selected for labelling),
\code{ns} (the (ordered) indices of the selected non-labelled observations).
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
