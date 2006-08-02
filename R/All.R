###
# NADA for R by Lopaka Lee.
#
# Version 1.3-0
# Copyright (2004, 2005, 2006) Lopaka Lee
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

### Definitions common to All sections in this package

## Globals

.onLoad = function(lib, pkg) 
{
    require(methods)
    #library.dynam("NADA", pkg, lib)
}

NADAprobs = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)

## Generics

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setGeneric("print", function(x, ...) standardGeneric("print"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("mean", function(x, ...) standardGeneric("mean"))

setGeneric("sd", function(x, na.rm=FALSE) standardGeneric("sd"))

setGeneric("median", function(x, na.rm=FALSE) standardGeneric("median"))

setGeneric("min", function(..., na.rm=FALSE) standardGeneric("min"))
setGeneric("max", function(..., na.rm=FALSE) standardGeneric("max"))

setGeneric("quantile", function(x, ...) standardGeneric("quantile"))

setGeneric("predict", function(object, ...) standardGeneric("predict"))

setGeneric("pexceed", function(object, ...) standardGeneric("pexceed"))

setGeneric("lines", function(x, ...) standardGeneric("lines"))

setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))

setGeneric("cor", function(x, y = NULL, use = "all.obs",
          method = c("pearson", "kendall", "spearman")) standardGeneric("cor"))

# LCL and UCL return string representations of 
# the lower and upper conf limits of an object (e.g. "0.95LCL")
setGeneric("LCL", function(x) standardGeneric("LCL"))
setGeneric("UCL", function(x) standardGeneric("UCL"))

## Broken for the time being
#setGeneric("abline", 
#           function(a, b, h, v, reg, coef, untf, col, lty, lwd, ...) 
#           standardGeneric("abline"))

setGeneric("residuals", function(object, ...) standardGeneric("residuals"))

setGeneric("coef", function(object, ...) standardGeneric("coef"))

## Maybe in the future transform() could be a generic 
#  Remember that transform is different in R minor versions < 3
#if (as.numeric(version$minor) < 3) {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(x, ...) standardGeneric("transform"))
#} else {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(`_data`, ...)
#                 standardGeneric("transform"))
#}


## Classes

setClass("NADAList", "list")

## Methods

#setMethod("summary", signature(), function(x, ...))

setMethod("print", signature("NADAList"), function(x, ...)
{
    tag = names(x)
    for (i in 1:length(x))
      {
        cat(tag[i], "\n")
        print(x[[i]])
        cat("\n")
      }
})

# Dennis' censored boxplots 
cenboxplot =
function(obs, censored, group, log=TRUE, range=0, ...) 
{
  if (log) log="y"
  else     log=""

  if (missing(group)) ret = boxplot(obs, log=log, range=range, ...)
  else                 ret = boxplot(obs~group, log=log, range=range, ...)

  # Draw horiz line at max censored value
  abline(h=max(obs[censored])) 

  invisible(ret)
}

# Dennis' censored xy plots -- need to fix this
cenxyplot =
function(x, xcen, y, ycen, log="xy", ...) 
{
    # Setup plot
    plot(x, y, log=log, type="n", ...)
    # Plot uncensored values
    points(x[!ycen], y[!ycen], ...)
    # Plot censored values
    points(x[ycen], y[ycen], ...)
}

# A first-cut at a summay function for censored data.  To do: groups.
censummary =
function(obs, censored) 
{
    smry = 
    function(obs, cen)
    {
        ret = cohn(obs, censored)

        ret$n = length(obs)
        ret$n.cen = length(obs[censored])

        cat("Summary:\n")
        props = c(ret$n, ret$n.cen, pctCen(obs, censored), min(obs), max(obs))
        names(props) = c("n", "n.cen", "pct.cen", "min", "max")
        print(props)

        limits = t(data.frame(ret$C, ret$A, ret$P))
        colnames(limits) = ret$limit
        rownames(limits) = c("C", "A", "P")

        cat("\nThresholds and counts:\n")
        print(limits["C",])

        cat("\nUncensored between each threshold:\n")
        print(limits["A",])

        cat("\nROS probability of exceeding each threshold:\n")
        print(limits["P",])
    }
    ret = smry(obs, censored)

    invisible(ret)
}

# Broken until I can figure out the eval.parent problem in cenreg
censtats =
function(obs, censored) 
{
    #stop("This convenience function is currently broken ... sorry")

    skm  = cenfit(obs, censored)
    sros = cenros(obs, censored)
    smle = cenmle(obs, censored)

    med  = c(median(skm), median(sros), median(smle))
    sd   = c(sd(skm), sd(sros), sd(smle)[1])
    mean = c(mean(skm)[1], mean(sros), mean(smle)[1])

    len = c(length(obs), length(which(censored == T)), pctCen(obs, censored))
    names(len) = c("n", "n.cen", "pct.cen")
    print(len)

    data.frame(median=med, mean=mean, sd=sd, row.names=c("K-M", "ROS", "MLE"))
}


#-->> BEGIN general utility functions

##
# split_qual extracts qualifed and unqualifed vectors from a vector
# containing concatenated qualifiying characters and numeric values
# like "<0.5".  Only handles one kind of censoring character/symbol.
splitQual =
function(v, qual.symbol = "<")
{
    v = as.character(v)

    obs = as.numeric(sub(qual.symbol, "", v))
    cen = rep(FALSE, length(obs))
    cen[grep(qual.symbol, v)] = TRUE 

    return(list(obs=obs, cen=cen))
}

## pctCen -- Simple function to save some typing
pctCen =
function(obs, censored)
{
    if (!is.logical(censored)) stop("censored must be logical vector!\n")

    return(100*(length(obs[censored])/length(obs)))
}

#-->> END general utility functions

