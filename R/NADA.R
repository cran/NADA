###
# NADA for R and S-PLUS by Lopaka Lee.
#
# Version 1.0-2
# Copyright (2004) Lopaka Lee
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


## Generic methods for ros objects

###
# Here we make some S3 generics that don't usually exist.
#
# This shuts up the warning we get for masking functions in the base package
.conflicts.OK = 1
#
# abline() as a S3 generic -- see abline.ros below 
#   Note that abline is in graphics:: in R > 1.9.x
abline.default =
    if (exists("abline", where = NULL, inherits = FALSE, mode = "function"))
        base::abline else graphics::abline

abline = 
function (a = NULL, b = NULL, h = NULL, v = NULL, reg = NULL,
          coef = NULL, ...) UseMethod("abline")


# sd() as a S3 generic -- see sd.ros below 
#   Note that sd is in stats:: in R > 1.9.x
sd.default = 
    if (exists("sd", where = NULL, inherits = FALSE, mode = "function"))
        base::sd else stats::sd 

sd = function(x, na.rm) UseMethod("sd")
#
# median() as a S3 generic -- see median.ros below
#   Note that median is in stats:: in R > 1.9.x
median.default = 
    if (exists("sd", where = NULL, inherits = FALSE, mode = "function"))
        base::median else stats::median 

median = function(x, na.rm) UseMethod("median")
#
#predict.default = predict
#
###

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
    ret = summary.lm(object)
    if (plot)
      {
        oldClass(object) = "lm"
        plot.lm(object)
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

# .cohN Calculates 'Cohn' numbers -- quantities described by 
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


## General utility functions

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
    return(100*(length(obs[censored])/length(obs)))
}
