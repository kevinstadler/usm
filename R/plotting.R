prepareplots <- function(nplots) {
  par(mfrow=c(nplots%/%2, if (nplots>1) 2 else 1), mar=c(4.1, 3.1, 3.1, 2.5), pty="s", xaxs="i", yaxs="i", mgp=c(2, 0.7, 0))
}

newplot <- function(xlim=0:1, xlab="interactions/agent", ...)
  plot(xlim, 0:1, type="n", xlab=xlab, ylab="x", xaxs="i", yaxs="i", ...)

addmomentumaxis <- function(col="gray", ...) {
    axis(4, c(0, 1/4, 1/2, 3/4, 1), c(-1.0, -0.5, 0, 0.5, 1.0), col=col)
    mtext("m\'", side=4, 2, col=col)
}

addmomentumline <- function(col="gray", lty=2, ...)
  abline(h=0.5, col=col, lty=lty, ...)

addtrajectory <- function(xs, data, maxmomentum, dataindices=c("xrange", "x", "mrange", "m"), invert=FALSE, alpha=0.5, ...) {
  for (dataindex in dataindices) {
    if (dataindex == "xrange") {
      if (invert)
        polygon(c(xs, rev(xs)), c(1-data[,"xmin"], rev(1-data[,"xmax"])), col=gray(0.4, alpha), border=NA)
      else
        polygon(c(xs, rev(xs)), c(data[,"xmin"], rev(data[,"xmax"])), col=gray(0.4, alpha), border=NA)
    } else if (dataindex == "mrange") {
      polygon(c(xs, rev(xs)), c(0.5-sign(as.numeric(invert)-0.5)*data[,"mmin"]/(2*maxmomentum), rev(0.5-sign(as.numeric(invert)-0.5)*data[,"mmax"]/(2*maxmomentum))), col=gray(0.8, alpha), border=NA)
    } else if (substr(dataindex,1,1) == "x" || dataindex <= 4) {
      if (invert) {
        lines(xs, 1-data[,dataindex], ...)
      } else {
        lines(xs, data[,dataindex], ...)
      }
    } else if (substr(dataindex,1,1) == "m" || dataindex > 4) {
      lines(xs, 0.5-sign(as.numeric(invert)-0.5)*data[,dataindex]/(2*maxmomentum), col=gray(0.3, 0.5), ...)
    } else {
      warning("Unknown data plotting index: ", dataindex)
    }
  }
}

# turn a vector of parameter specifications into an output filename
resultfilename <- function(pars) {
  do.call(paste, c(as.list(mapply(c, parnames, pars)), sep=''))
}


plotmomentum <- function(data, xindex=1, xmin=2, xmax=3, mindex=4, mnormalisation=1, t=1:nrow(data), xlab="interactions/agent", alpha=0.3, ...) {
  plot(range(t), 0:1, type="n", xaxs="i", yaxs="i", xlab=xlab, ylab="x")
  lines((1+data[,mindex])/(2*mnormalisation), col="gray")
  lines(data[,xindex])
  axis(4, c(0, 1/4, 1/2, 3/4, 1), c(-1, -0.5, 0, 0.5, 1), col="gray")
  mtext("m\'", side=4, 2, col="gray")
}

plotgenerational <- function(data, alpha=0.3) {
  plot(c(1, nrow(data)), 0:1, type="n", xaxs="i", yaxs="i", xlab="time (years)", ylab="Proportion incoming variant")
  polygon(c(1:nrow(data),nrow(data):1), c(data$maxold, rev(data$minold)), col=rgb(0,0,1,alpha), border=FALSE)
  polygon(c(1:nrow(data),nrow(data):1), c(data$maxyoung, rev(data$minyoung)), col=rgb(1,0,0,alpha), border=FALSE)
  polygon(c(1:nrow(data),nrow(data):1), (1+c(data$maxmomentum, rev(data$minmomentum)))/2, col=rgb(0,1,0,alpha), border=FALSE)
  polygon(c(1:nrow(data),nrow(data):1), c(data$maxx, rev(data$minx)), col=gray(0.2,alpha), border=FALSE)
  lines(data$old, col="blue")
  lines(data$young, col="red")
  lines(data$x)
  lines(data$momentum, col="darkgreen")
  lines((1+data$young-data$old)/2, col="green")
}
# par(mfrow=c(1,2))
# plotgenerational(read.table("generational.csv", header=TRUE))
# plotgenerational(read.table("generational2.csv", header=TRUE))

plottransition <- function(transition, ...) {
  data <- read.table(paste(transition$run, resultfilename(transition[parnames]), sep="/"))
  plot(transition$start:transition$end, data[transition$start:transition$end,1], type="l", ylim=0:1, xlab="interactions/agent", ylab="x", xaxs="i", yaxs="i")
  lines(transition$start:transition$end, data[transition$start:transition$end,3], lty=2)
  lines(transition$start:transition$end, data[transition$start:transition$end,4], lty=2)
}
