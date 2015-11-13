source("~/usm/R/plotting.R")
setwd("data")
source("conditions.R")

# all parameter combinations
allconditions <- as.matrix(expand.grid(N, b, t, x0, i, a, m))
parnames <- c("N", "b", "T", "x0", "i", "a", "m")
colnames(allconditions) <- parnames
allconditions <- data.frame(allconditions)

# column names of the individual runs
colnames <- c("x", "xsd", "xmin", "xmax", "m", "msd", "mmin", "mmax")

readtrajectory <- function(run, condition)
  read.table(paste(run, resultfilename(condition), sep="/"), col.names=colnames)

# perform a simple segmentation based on crossing thresholds
# returns a dataframe with columns start, end, direction (1 for up, -1 for down), success (0 for interrupted changes)
segmentrun <- function (data, threshold) {
  nonlow <- data > threshold
  nonhigh <- data < (1-threshold)
  # T=1, F=0 -> -1=low, 0 = mid, 1=high
  data <- nonlow - nonhigh
  crossings <- abs(diff(data))
  run <- NULL
  t <- 2
  lastextreme <- data[t]
  while (t < length(data)) {
    nextchange <- match(1, crossings[t:length(data)])
    if (is.na(nextchange)) {
      nextchange <- length(data)
    } else {
      nextchange <- nextchange+t
      if (data[nextchange] != 0) {
        # something interesting happened
        run <- rbind(run, c(start=t, end=nextchange, direction=-lastextreme, success=data[nextchange]!=lastextreme))
        lastextreme <- data[nextchange]
      }
    }
    t <- nextchange
  }
  return(data.frame(run))
}

stepstomaxdifference <- function(alpha, gamma)
  log(alpha / gamma) / (alpha - gamma)

expdecay <- function(lambda, t)
  exp(lambda*-t)

maxdecaydifference <- function(alpha, gamma) {
  t  <- stepstomaxdifference(alpha, gamma)
  expdecay(alpha, t) - expdecay(gamma, t)
}

simplesummary <- function(conditions=allconditions, runs=1:24) {
  alltransitions <- data.frame()
  for (run in runs) {
    print(paste("Run", run))
    for (i in 1:nrow(conditions)) {
      data <- readtrajectory(run, conditions[i,])
      transitions <- segmentrun(data$x, 0.05)
      ntransitions <- nrow(transitions)
      #ncompleted <- 0
      if (ntransitions > 0) {
        #ncompleted <- nrow(subset(transitions, success==1))
        for (j in 1:nrow(transitions)) {
          f <- if (transitions[j,"direction"] == 1) max else min
          transitions[j,"t0"] <- if (transitions[j,"success"] == 0)  NA else match(transitions[j,"direction"]==1, data$x[transitions[j,"start"]:transitions[j,"end"]] >= 0.5)-1
          transitions[j,"mmax"] <- f(data$m[transitions[j,"start"]:transitions[j,"end"]])/maxdecaydifference(conditions[i,"a"], conditions[i,"a"]*conditions[i,"m"])
          transitions[j,"xmax"] <- f(data$x[transitions[j,"start"]:transitions[j,"end"]])
          alltransitions <- rbind(alltransitions, c(run=run, conditions[i,], transitions[j,]))
        }
      }
    } # incremental writeout
    write.table(alltransitions, "~/simpletransitions.csv", row.names=FALSE, append=TRUE)
  }
  return(alltransitions)
}

if (file.exists("~/simpletransitions.csv")) {
  alltransitions <- read.table("~/simpletransitions.csv", header=TRUE)
} else {
  library(doParallel)
  registerDoParallel(cores=8)
  system.time(alltransitions <- foreach(run=1:24, .combine=rbind) %dopar% simplesummary(allconditions, run))
  write.table(alltransitions, "~/simpletransitions.csv", row.names=FALSE)
}
alltransitions$xmax <- abs(1-alltransitions$poslow-alltransitions$xmax)
completed <- subset(alltransitions, success==1)
interrupted <- subset(alltransitions, success==0)
