###
# NADA for R by Lopaka Lee.
#
# Version 1.1-2
# Copyright (2004, 2005) Lopaka Lee
#
# A S-language software module based on 
# methodologies described by Dennis R. Helsel in his book 
# Nondetects and Data Analysis: Statistics for Censored Environmental Data.
#
# NADA is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
# 
# NADA is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.  You should have received a copy of the GNU General
# Public License along with NADA; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
###


###-->> BEGIN R-package initialization routines

.First.lib =
function(libname, pkgname)
{
    require(survival, quietly=TRUE)
}

###-->> END R-package initialization routines

###-->> BEGIN Generic methods for S/R objects

# Here we make some S3 generics that don't usually exist.
#
# This shuts up the warning we get for masking functions in the base package
.conflicts.OK = 1
#
# abline() as a S3 generic
#   Note that abline is in graphics:: in R > 1.9.x
abline.default =
    if (exists("abline", where = NULL, inherits = FALSE, mode = "function"))
        base::abline else graphics::abline

abline = 
function (a = NULL, b = NULL, h = NULL, v = NULL, reg = NULL,
          coef = NULL, ...) UseMethod("abline")

# mean() as a S3 generic
#   Note that mean is in stats:: in R > 1.9.x
mean.default = 
    if (exists("mean", where = NULL, inherits = FALSE, mode = "function"))
        base::mean.default else stats::mean

mean = function(x, ...) UseMethod("mean")

# sd() as a S3 generic
#   Note that sd is in stats:: in R > 1.9.x
sd.default = 
    if (exists("sd", where = NULL, inherits = FALSE, mode = "function"))
        base::sd else stats::sd 

sd = function(x, na.rm) UseMethod("sd")


# median() as a S3 generic
#   Note that median is in stats:: in R > 1.9.x
median.default = 
    if (exists("median", where = NULL, inherits = FALSE, mode = "function"))
        base::median else stats::median 

median = function(x, na.rm) UseMethod("median")
#
#
###-->> END Generic methods for S/R objects


###-->> BEGIN Regression on Order Statistics (ROS) functions

###-->> BEGIN General utility functions

##
# splitQual extracts qualifed and unqualifed vectors from a vector
# containing concatenated qualifiying characters and numeric values
# like "<0.5". See man page.
splitQual =
function(v, qual.symbol = "<")
{
    qual.index = grep(qual.symbol, x=as.character(v))  

    qual.chars = as.character(v[qual.index])
    qual = as.numeric(sub(qual.symbol, "", qual.chars))

    unqual.index = -1 * qual.index
    unqual = as.numeric(as.character(v[unqual.index]))
    
    return(list(qual         = qual, 
                unqual       = unqual, 
                qual.index   = qual.index,
                unqual.index = unqual.index))
}

## pct.censored  -- Simple function to save some typing
pct.censored =
function(obs, censored)
{
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    return(100*(length(obs[censored])/length(obs)))
}

###-->> END General utility functions

# ros -- Regression on Order Statistics (ROS).
# An implementation of ROS for left-censored data (less-thans) 
# containing one to multiple censoring thresholds. See man page.
ros = 
function(obs, 
         censored, 
         forwardT = "log", 
         reverseT = "exp")
{
    if (is.null(forwardT) || is.null(reverseT)) {
        forwardT = reverseT = ".trueT"
    }
    else if (!exists(forwardT)) {
         stop("Can not find Forward Transformation function: ", forwardT, "\n")
    }
    else if (!exists(reverseT)) {
        stop("Can not find Reverse Transformation function: ", reverseT, "\n")
    }

    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    if ( (length(obs[censored])/length(obs)) > 0.8 ) {
        warning("Input > 80% censored -- Results are tenuous.\n")
    }

    ix = order(obs)
    obs = obs[ix]
    censored = censored[ix]

    pp = hc.ppoints(obs, censored)

    pp.nq = qnorm(pp[!censored])
    obs.transformed = get(forwardT)(obs[!censored])
    hc = lm(obs.transformed ~ pp.nq)

    oldClass(hc) = c("ros", "lm")
    hc$obs      = obs
    hc$modeled  = obs 
    hc$pp       = pp
    hc$censored = censored 
    hc$reverseT = reverseT
    hc$forwardT = forwardT

    hc$modeled[censored] = predict.ros(hc, qnorm(pp[censored]))

    return(hc)
}

# ROS is technically "Linear Regression on Order Statistics" 
# or LROS.  Version 1.0-1 had funtions and objects as "lros".  
# As of Version 1.0-2, I've renamed things to "ros" to keep consistent
# with the literature on this method which calls it just ROS.
# This is provided for backward compatibility.
lros = ros

#  .trueT is provided so that ros() can be used with no transforms.
#  (It is a quick hack -- rather than recoding ros and family)
.trueT =
function(x)
{
    return(x)
}


print.ros =
function(x, ...)
{
    cat("\nMultiply-Censored ROS Model\n\n")

    uncen.n = length(x$modeled[!x$censored])
    cen.n   = length(x$modeled[x$censored])

    cat("           N:", (uncen.n + cen.n), "\n")
    cat("    Censored:", cen.n, "\n")
    cat("  % Censored:", format((cen.n/(cen.n+uncen.n))*100, digits=4), "\n")
    cat("\n");
    
    cat("        Mean:", format(mean.ros(x), digits=4), "\n")
    cat("      StdDev:", format(sd.ros(x), digits=4), "\n")
    cat("\n");
    cat("   Quantiles:\n")
    print(quantile(x))
    cat("\n");
    cat("Use summary() to view the linear regression model\n\n")
}

summary.ros =
function(object, plot=FALSE, ...)
{
    ret = summary.lm(object, ...)
    if (plot)
      {
        oldClass(object) = "lm"
        plot.lm(object, ...)
      }
    return(ret)
}

# as.data.frame conversion discards all linear model info
as.data.frame.ros =
function (x, row.names = NULL, optional = FALSE)
{
    x = list(obs=x$obs, censored=x$censored, pp=x$pp, modeled=x$modeled)
    as.data.frame(x, row.names, optional)
}

mean.ros =
function(x, ...)
{
    mean(x$modeled, ...)
}

sd.ros =
function(x, na.rm=FALSE)
{
    sqrt(var(x$modeled, na.rm = na.rm))
}


median.ros =
function(x, na.rm=FALSE)
{
    median(x$modeled, na.rm)
}

## Query and prediction functions for ROS objects

quantile.ros =
function(x, probs=c(0.05,0.10,0.25,0.50,0.75,0.90,0.95),...)
{
    quantile(x$modeled, probs, ...)
}

# predict.ros -- given new normal quantiles of plotting positions,
# returns the corresponding modeled values.
predict.ros =
function(object, newdata, ...)
{
    predicted = 
      as.vector(predict.lm(object, newdata=data.frame(pp.nq=newdata), ...))

    return(get(object$reverseT)(predicted))
}

## Generic plotting functions

# Adds the line from a ROS model to an existing ROS prob-plot.
abline.ros =
function (a = NULL, ...)


{
    # We can't delegate this to the default abline because of transformations.
    minmax.nq = qnorm(c(min(a$pp), max(a$pp)))
    lines(x=minmax.nq, y=predict.ros(a, minmax.nq), ...)
}

# Constructs a prob-plot representation of a ROS model
plot.ros =
function(x, 
         plot.censored = FALSE, 
         lm.line = TRUE, 
         grid    = TRUE, 
         ylab    = "Value", 
         pch     = 16,
         ... )
{
    ##
    # To do:
    #   Refactor this routine -- it is long and ugly.
    #   Constrain y ticks to factors of 10?
    #   Draw in grid lines at y ticks
    #
    uncen = x$modeled[!x$censored]
    cen   = x$modeled[x$censored]

    pp.uncen.nq = qnorm(x$pp[!x$censored])
    pp.cen.nq   = qnorm(x$pp[x$censored])

    ymin = min(c(uncen, cen))
    ymax = max(c(uncen, cen))
    xmin = min(c(pp.uncen.nq, pp.cen.nq))
    xmax = max(c(pp.uncen.nq, pp.cen.nq))

    if (x$forwardT == "log") {
        plot(y    = uncen, 
             x    = pp.uncen.nq, 
             ylim = c(ymin, ymax), 
             xlim = c(xmin, xmax), 
             ylab = ylab,
             xlab = "Normal Quantiles",
             log  = "y",
             yaxt = "n",
             pch  = 16,
             ...
        )
     }
     else {
        plot(y    = uncen, 
             x    = pp.uncen.nq, 
             ylim = c(ymin, ymax), 
             xlim = c(xmin, xmax), 
             ylab = ylab,
             xlab = "Normal Quantiles",
             yaxt = "n",
             pch  = 16,
             ...
        )
     }

    axis(2, axTicks(2))

    # Draw a line through the predicted xmin and xmax points 
    if (lm.line)
      {
        abline.ros(x) 
      }

    # Draw top "Chance of Exceedance" axis
    labels = c("5","10","25","50","75","90","95")
    atv = qnorm(c(.05, .1, .25, .50, .75, .90, .95))
    axis(3, at=atv, labels=rev(labels), las=2)
    mtext("Percent Chance of Exceedance", side=3, line=3)

    # Plot the synthetic censored points -- if requested
    if (plot.censored) points(y=cen, x=pp.cen.nq)

    # Draw in grid at major divisions  -- if requested
    if (grid)
      {
        abline.default(v=atv, lty="dotted")
        abline.default(h=axTicks(2), lty="dotted")
      }
}


## Routines for calculating Helsel-Cohn style plotting positions

# hc.ppoints calculates computes Wiebull-type plotting postions of data
# containing mixed uncensored and left-censored data. See man page.
# If there are no censored values, the plotting postitions are calculated
# using the standard ppoints() function.
hc.ppoints = 
function(obs, censored)
{    
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    pp = numeric(length(obs))

    if (!any(censored)) pp = ppoints(obs)
    else
      {
        uncen = obs[!censored]
        cen   = obs[censored]

        cn = .cohnN(uncen, cen)
        pp[!censored] = hc.ppoints.uncen(cn=cn)
        pp[censored]  = hc.ppoints.cen(cn=cn)
      }

    return(pp)
}

# .cohnN Calculates 'Cohn' numbers -- quantities described by 
# Helsel and Cohn's (1988) reformulation of a prob-plotting formula
# described by Hirsch and Stedinger (1987).
.cohnN =
function(uncen, cen)
{
    alldata = c(uncen, cen)
    A = B = C = P = numeric()

    limit = sort(unique(cen))

    a = length(uncen[uncen < limit[1]])
    if (a > 0) { limit = c(0, limit) }

    i = length(limit)

    A[i] = length(uncen[ uncen >= limit[i] ])
    B[i] = length(alldata[alldata <= limit[i]])
    C[i] = length(cen[ cen == limit[i] ])
    P[i] = A[i]/(A[i] + B[i])

    i = i - 1
    while (i > 0)
      {
        A[i] = length(uncen[ uncen >= limit[i] & uncen < limit[i + 1] ])
        B[i] = length(alldata[alldata <= limit[i]])
        C[i] = length(cen[ cen == limit[i] ])
        P[i] = P[i + 1] + ((A[i]/(A[i] + B[i])) * (1 - P[i + 1]))

        i = i - 1
      }
    return(list(A=A, B=B, C=C, P=P, limit=limit))
}

# hc.ppoints.uncen calculates plotting postions for uncensored data.
hc.ppoints.uncen =
function(obs, censored, cn)
{
    if (missing(cn)) { cn = .cohnN(obs, censored) }

    nonzero = (cn$A != 0)
    A     = cn$A[nonzero]
    B     = cn$B[nonzero]
    P     = cn$P[nonzero]
    limit = cn$limit[nonzero]

    pp = numeric()
    n = length(limit)
    for (i in 1:n)
      {
        R = 1:A[i] 

        k = P[i+1]
        if (is.na(k)) k = 0

        for (r in 1:length(R))
          {
            pp = c(pp, (1 - P[i]) + ((P[i] - k) * R[r])/(A[i] + 1))
          }
      }
    return(pp)
}

# hc.ppoints.cen calculates plotting postions for censored data.
hc.ppoints.cen =
function(obs, censored, cn)
{    
    if (missing(cn)) { cn = .cohnN(obs, censored) }

    C     = cn$C
    P     = cn$P
    limit = cn$limit

    if (P[1] == 1) 
      {
        C     = C[-1]
        P     = P[-1]
        limit = limit[-1]
      }
    
    pp = numeric()
    for (i in 1:length(limit)) 
      {
        c.i = C[i]
        for (r in 1:c.i) 
          {
            pp = c(pp, (1 - P[i]) * r/(c.i + 1))
          }
      }
    return(pp)
}

###-->> END Regression on Order Statistics (ROS) functions

###-->> BEGIN Survival Analysis-based functions -- 'cencode' part of NADA4R

# utility to create ECDFs that are useful to us
cenECDF = function(x, y) 
{
    y = y[order(x)]
    x = sort(x)
    return(stepfun(x[-1], y))
}

# Cen() is analgous to Surv() in survival package.
# The difference is that Cen() "flips" the our left-censored data
# into right-censored data so that the survival code can be used.
Cen =
function(obs, censored, type="left", ...) 
{
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    m = Surv(time=obs, time2=(!censored), type=type, ...)

    oldClass(m) = c("Cen", class(m))
    attr(m, 'flipFactor') = max(obs) 

    if (type == "left") m = flip(m)

    return(m)
}

# flip() a Cen object back into its original scale.
# All functions transparently do this for the user.
flip =
function(x) 
{
    if (!is(x, "Cen")) stop ("Input must be a Cen object");

    x[,1] = attr(x, 'flipFactor') - x[,1]

    if      (attr(x, 'type') == "right") attr(x, 'type') =  "left"
    else if (attr(x, 'type') == "left")  attr(x, 'type') =  "right"
    else stop ("Can only flip \"left\" or \"right\" censored Cen objects")

    return(x)
}

cenfit =
function (formula, data=NULL, ...) 
{
    ## This turns a Cen object in a proper model object.
    #  Typical S-Language ugliness follows ...
    call <- match.call()
    if ((mode(call[[2]]) == "call" && 
         call[[2]][[1]] == as.name("Cen")) ||
         inherits(formula, "Cen")) 
      {
        formula = eval(parse(text=paste(deparse(call[[2]]), 1, sep = "~")))
        environment(formula) <- parent.frame()
      }
    m = match.call(expand = FALSE)
    m$formula = terms(formula)
    m[[1]] = as.name("model.frame")

    if (is.null(data)) m = eval(m, parent.frame())
    else m = eval(m, data, parent.frame())
    ##

    # construct the survival curve
    sf = NULL
    if (is.null(data)) sf = survfit(formula, ...)
    else sf = survfit(formula, data=data, ...)

    # flip results back into the expected scale
    cenObj = model.response(m)
    sf$time = attr(cenObj, "flipFactor") - sf$time

    oldClass(sf) = c("cenfit", "survfit")

    return(sf)
}


# Indexing cenfit objects -- 
# When this happens, we throw away the strata information.
"[.cenfit" = 
function(x, i, j, drop=F) 
{
    if (is.null(x$strata)) stop("can't index; object contains only one ECDF")

    x$n = x$strata.all[i]

    class(x) = "survfit"
    x = NextMethod("[")
    class(x) = c("cenfit", "survfit")

    x$strata        = NULL
    x$strata.all    = NULL
    x$ntimes.strata = NULL

    return (x)
}

summary.cenfit =
function(object, times, censored=FALSE, scale=1, ...)
{
    class(object) = "survfit"
    x = summary(object, times, censored, scale, ...)
    # To do: modify the call object to reflect cenfit call
    # for now, just nullify it -- users typically remember the call.
    x$call = NULL
    # Reverse vectors -- our world sees things from the left ;-)
    x = lapply(x, rev)
    class(x) = c("summary.cenfit", "summary.survfit")
    return(x)
}

plot.cenfit =
function(x, 
         log  = 'x', 
         axLimitFactor = 0.8, 
         ylab = "Probability", 
         xlab = "Value", 
         lty  = seq(1,6),
         ...)
{
    firstx = (min(x$time)*axLimitFactor)

    oldClass(x) = "survfit"
    plot(x, log=log, firstx=firstx, ylab=ylab, xlab=xlab, lty=lty, ...)

    # Draw a vertical line at the minimum observation -- this is
    # the lower extent of the step curve or ECDF.
    abline(v=min(x$time))

    abline(h=1.0)
}

# NADA internal function -- used by predict.cenfit and quantile.cenfit
.predict.cenfit =
function(newdata, ecdF, ecdF.l, ecdF.u, do.conf, conf.int) 
{
    obj  = NULL
    pred = ecdF(newdata)

    if (!do.conf) obj = pred
    else
      {
        pred.l = ecdF.l(newdata)
        pred.u = ecdF.u(newdata)

        obj = data.frame(newdata, pred, pred.l, pred.u)

        l.lbl = paste(conf.int, "LCL", sep='')
        u.lbl = paste(conf.int, "UCL", sep='')

        names(obj) = c("newdata", "prediction", l.lbl, u.lbl)
      }

    return(obj)
}

predict.cenfit =
function(object, newdata, do.conf=FALSE, ...) 
{
    x = object
    obj = NULL

    if (is.null(x$strata))
      {
        ecdF   = cenECDF(x$time, x$surv)
        ecdF.l = cenECDF(x$time, x$lower)
        ecdF.u = cenECDF(x$time, x$upper)
        
        obj = 
        .predict.cenfit(newdata, ecdF, ecdF.l, ecdF.u, do.conf, x$conf.int)
      }
    else
      {
        nstrata = length(x$strata)
        for (i in 1:nstrata)
          {
            ix = (rep(1:nstrata, x$strata) == i)

            ecdF   = cenECDF(x$time[ix], x$surv[ix])
            ecdF.l = cenECDF(x$time[ix], x$lower[ix])
            ecdF.u = cenECDF(x$time[ix], x$upper[ix])
            
            obj[[i]] = 
            .predict.cenfit(newdata, ecdF, ecdF.l, ecdF.u, do.conf, x$conf.int)

            if (!do.conf) names(obj[[i]]) = as.character(newdata)
          }
        names(obj) = names(x$strata)
      }

    return(obj)
}

quantile.cenfit =
function(x, 
         probs   = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95), 
         do.conf = FALSE, 
         ...)
{
    obj = NULL

    # Local utility to assign names to quantile vectors or dataframes
    nameObj =
    function(obj) 
    {
      if (do.conf) names(obj)[c(1,2)] = c("prob", "quantile")
      else
        {
          names(obj) = 
            paste(formatC(100 * probs, format = "fg", wid = 1), "%", sep='')
        }
      return(obj)
    }

    #
    if (is.null(x$strata))
      {
        inv   = cenECDF(x$surv, x$time)
        inv.l = cenECDF(x$surv, x$lower)
        inv.u = cenECDF(x$surv, x$upper)
        
        obj = .predict.cenfit(probs, inv, inv.l, inv.u, do.conf, x$conf.int)

        obj = nameObj(obj)
      }
    else
      {
        nstrata = length(x$strata)
        for (i in 1:nstrata)
          {
            ix = (rep(1:nstrata, x$strata) == i)

            inv   = cenECDF(x$surv[ix], x$time[ix])
            inv.l = cenECDF(x$surv[ix], x$lower[ix])
            inv.u = cenECDF(x$surv[ix], x$upper[ix])
            
            obj[[i]] = 
            .predict.cenfit(probs, inv, inv.l, inv.u, do.conf, x$conf.int)
            
            obj[[i]] = nameObj(obj[[i]])
          }
        names(obj) = names(x$strata)
      }

    return(obj)
}

.mean.cenfit =
function (stime, surv, n.risk, n.event) 
{
    mean = NULL
    varmean = NULL

    min.stime = min(stime)
    min.time = min(0, min.stime)
    n = length(stime)

    hh = c(ifelse((n.risk[-n] - n.event[-n]) == 0, 
                   0, 
                   n.event[-n]/(n.risk[-n] * (n.risk[-n] - n.event[-n]))), 0)

    ndead = sum(n.event)
    dif.time = c(diff(c(min.time, stime)), 0)

    if (!is.matrix(surv)) 
      {
        mean = dif.time * c(1, surv)
        varmean = sum(rev(cumsum(rev(mean))^2)[-1] * hh)
      }
    else 
      {
        n = nrow(surv)
        mean = dif.time * rbind(1, surv)
        if (n == 1) 
            temp = mean[2, , drop = FALSE]
        else temp = (apply(mean[(n + 1):2, , drop = FALSE], 2, 
            cumsum))[n:1, , drop = FALSE]
        varmean = c(hh %*% temp^2)
      }

    obj = c(sum(mean), sqrt(varmean))
    names(obj) = c("rmean", "se(rmean)")

    return(obj)
}

mean.cenfit =
function (x, ...) 
{
    obj = NULL

    if (is.null(x$strata))
      {
        obj = .mean.cenfit(x$time, x$surv, x$n.risk, x$n.event)
      }
    else
      {
        nstrata = length(x$strata)
        for (i in 1:nstrata)
          {
            ix = (rep(1:nstrata, x$strata) == i)

            obj[[i]] = 
              .mean.cenfit(x$time[ix], x$surv[ix], x$n.risk[ix], x$n.event[ix])
          }
        names(obj) = names(x$strata)
      }

    return(obj)
}

median.cenfit =
function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    quantile.cenfit(x, 0.5)
}

sd.cenfit =
function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    obj = NULL

    if (is.null(x$strata)) obj = (sqrt(x$n) * as.numeric(mean.cenfit(x)[2]))
    else
      {
        nstrata = length(x$strata)
        for (i in 1:nstrata)
          {
            n = as.numeric(x$strata.all[i])
            obj[[i]] = (sqrt(n) * as.numeric(mean.cenfit(x[i])[2]))
          }
        names(obj) = names(x$strata)
      }

    return(obj)
}

# Wrapper for survival::survdiff 
# default of rho = 1 means Peto-Peto is used.
cendiff = 
function(formula, rho = 1, ...)
#function(formula, data, subset, na.action, rho = 1, ...)
{
    x = survival::survdiff(formula, rho, ...)
    # To do: modify the call object to reflect cenfit call
    # for now, just nullify it -- users typically remember the call.
    x$call = NULL
    return(x)
}

print.cenfit =
function(x, ...) 
{
    obj = NULL

    summaryVec =
    function(x)
    {
      n      = x$n
      events = sum(x$n.event)
      median = median(x)
      mean   = mean(x)

      return(c(n, events, median, mean))
    }

    if (is.null(x$strata))
      {
        obj = summaryVec(x)
        names(obj) = c("n", "events", "median", "mean", "se(mean)")
      }
    else
      {
        obj = summaryVec(x[1])
        for (i in 2:length(x$strata))
          {
            obj = rbind(obj, summaryVec(x[i]))
          }
        colnames(obj) = c("n", "events", "median", "mean", "se(mean)")
        rownames(obj) = names(x$strata)
      }

    print(obj, ...)
    invisible(obj)
}

#####
#####
# What follows is hydras code from the survival package.
# For now, the easiest thing to do is MODify it for our needs.

# MODified survival::summary.survfit

# MODified survival::print.summary.survfit
# The labeling of output vectors/columns is burried in here.  
# We want labeling to be consistent with our world view.
print.summary.cenfit =
function(x, 
         digits = max(options()$digits - 4, 3), 
         ...) 
{
    savedig <- options(digits=digits)
    on.exit(options(savedig))

## MOD: We can remember what the call was, thank you.
##    if (!is.null(cl<- x$call)) {
##	cat("Call: ")
##	dput(cl)
##	cat("\n")
##	}

    omit <- x$na.action
    if (length(omit)) 
	    cat(naprint(omit), "\n")
    if (x$type == 'right' || is.null(x$n.entered)) {
	mat <- cbind(x$time, x$n.risk, x$n.event, x$surv)
# MOD: time --> obs/observation
	cnames <- c("obs", "n.risk", "n.event")
        }

    else if (x$type == 'counting') {
	mat <- cbind(x$time, x$n.risk, x$n.event, x$n.entered,
		     x$n.exit.censored, x$surv)
# MOD: time --> obs/observation
	cnames <- c("obs", "n.risk", "n.event", 
		    "n.entered", "n.censored")
        }
    if (is.matrix(x$surv)) ncurve <- ncol(x$surv)
    else	           ncurve <- 1
    if (ncurve==1) {                 #only 1 curve
# MOD: survival --> prob/probability
	cnames <- c(cnames, "prob")
	if (!is.null(x$std.err)) {
	    if (is.null(x$lower)) {
		mat <- cbind(mat, x$std.err)
		cnames <- c(cnames, "std.err")
	        }
	    else {
		mat <- cbind(mat, x$std.err, x$lower, x$upper)
		cnames <- c(cnames, 'std.err',
			    paste("lower ", x$conf.int*100, "% CI", sep=''),
			    paste("upper ", x$conf.int*100, "% CI", sep=''))
	        }	
	    }
        }
# MOD: survival --> prob/probability
    else cnames <- c(cnames, paste("prob", seq(ncurve), sep=''))

    if (!is.null(x$new.start)) {
	mat.keep <- mat[,1] >= x$new.start
	mat <- mat[mat.keep,,drop=FALSE]
	if (is.null(dim(mat)))
		stop(paste("No information available using new.start =", x$new.start, "."))
        }
    if (!is.matrix(mat)) mat <- matrix(mat, nrow=1)
    if (!is.null(mat)) {
	dimnames(mat) <- list(NULL, cnames)
	if (is.null(x$strata))
		prmatrix(mat, rowlab=rep("", nrow(mat)))
	else  { #print it out one strata at a time
	    if (!is.null(x$times.strata))
		    strata <- x$times.strata
	    else
		    strata <- x$strata
	   
	    if (!is.null(x$new.start))
		    strata <- strata[mat.keep]
	    for (i in levels(strata)) {
		who <- (strata==i)
		cat("               ", i, "\n")
		if (sum(who) ==1)
			print(mat[who,])
	        else
		    prmatrix(mat[who,], rowlab=rep("", sum(who)))

		cat("\n")
 	        }
	    }
        }
    else 
	stop("There are no events to print.  Please use the option ",
	    "censored=TRUE with the summary function to see the censored ",
	    "observations.")
    invisible(x)
}


###-->> END Survival Analysis-based functions
