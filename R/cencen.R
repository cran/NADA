#-->> BEGIN survival analysis based functions -- "cencen" section 

## Generics

setGeneric("cenfit", 
           function(obs, censored, groups, ...) standardGeneric("cenfit"))

setGeneric("cendiff", 
           function(obs, censored, groups, ...) standardGeneric("cendiff"))

setGeneric("flip", function(x) standardGeneric("flip"))

## Classes 

setOldClass("Surv")
setOldClass("survfit")

setClass("Cen", representation(Surv="Surv", flipFactor="numeric"))
setClass("cenfit", representation(survfit="survfit"))

## Methods

# Begin cencen utility functions

# LCL() and UCL() return string representations of 
# the lower and upper conf limits of cenfit object ("0.95LCL")
LCL =
function(x) { paste(x@survfit$conf.int, "LCL", sep='') }

UCL =
function(x) { paste(x@survfit$conf.int, "UCL", sep='') }

# stepfind() -- More applicable version of stats::stepfun(). 
# Used in predict.cenfit() and quantile.cenfit().
# Given decreasingly ordered x and y vectors of a step function, 
# find the y value associated with any given x value.  
# Optionally right or left looking on the number line. 
stepfind =
function(x, y, val, right=TRUE) 
{
    findStep =
    function(x, y, val, right)
    {
        i = length(x)
        if (val >= max(x)) return(max(y))
        if (val <= min(x)) return(min(y))

        while (as.logical(i)) 
          {
            if (x[i] > val || identical(all.equal(x[i], val), TRUE)) break
            i = i - 1
          }
        return(ifelse(right, y[i], y[i+1]))
    }

    sapply(val, findStep, x=x, y=y, right=right)
}

## End cencen utility functions

# Cen() is analgous to Surv() in survival package.
# However, Cen() can be used in the flip() routines below to
# rescale left-censored data into right-censored data 
# for use in the survival package routines.
Cen =
function(obs, censored, type="left")
{
    new("Cen", Surv=Surv(obs, !censored, type=type), flipFactor = max(obs))
}

setMethod("print", signature(x="Cen"), function(x, ...) print(x@Surv))

## flip() methods for rescaling input for use in survival routines.

# Note that flip()ing a Cen object results in a Surv object.
setMethod("flip", signature(x="Cen"), function(x)
{
    surv = x@Surv
    # Subtract the flipFactor from all observations
    surv[,1] = x@flipFactor - surv[,1]

    # Update the censoring type
    if (attr(surv, "type") == "right") attr(surv, "type") = "left"
    else if (attr(surv, "type") == "left") attr(surv, "type") = "right"
    else stop("Can only flip \"left\" or \"right\" censored Cen objects")

    return(surv)
})

# flip()ing a formula just symbolically updates the 
# response (which should be a Cen object). 
# Result is like: flip(Cen(obs, cen))~groups
setMethod("flip", signature(x="formula"), function(x) update(x, flip(.)~.))

## Begin cencen.* functions
# These routines are allow cenfit and cendiff methods
# to be used with Cen objects, vectors, or formulas as input.
# Note they all convert input to a formula and call the formula method.

cencen.Cen =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(a~1, list(a=cl[[2]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cencen.vectors =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(Cen(a, b)~1, list(a=cl[[2]], b=cl[[3]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cencen.vectors.groups =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = substitute(Cen(a, b)~g, list(a=cl[[2]], b=cl[[3]], g=cl[[4]]))
    f = as.formula(f)
    environment(f) = parent.frame()
    callGeneric(f, ...)
}
## End cencen.* routines

# cenfit for formulas
setMethod("cenfit", 
          signature(obs="formula", censored="missing", groups="missing"), 
          function(obs, censored, groups, ...)
{
    sf = survfit(flip(obs), ...)

    cenObj = eval(obs[[2]], environment(obs))
    sf$time = cenObj@flipFactor - sf$time 

    new("cenfit", survfit=sf)
})

setMethod("cenfit", 
          signature(obs="Cen", censored="missing", groups="missing"), 
          cencen.Cen)

setMethod("cenfit", 
          signature(obs="numeric", censored="logical", groups="missing"), 
          cencen.vectors)

setMethod("cenfit", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)

# cendiff for formulas
setMethod("cendiff", 
          signature(obs="formula", censored="missing", groups="missing"), 
          function(obs, censored, groups, rho=1, ...)
{
    x = survival::survdiff(flip(obs), rho=rho, ...)
    # To do: modify the call object to reflect cenfit call
    # for now, just nullify it -- users typically remember the call.
    x$call = NULL
    return(x)
})

#setMethod("cendiff", 
#          signature(obs="Cen", censored="missing", groups="missing"), 
#          cencen.Cen)

setMethod("cendiff", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)

# Indexing cenfit objects.
# When this happens, we throw away the strata information.
setMethod("[", signature(x="cenfit", i="numeric", j="missing"), 
          function(x, i, drop=FALSE) 
{
    s = x@survfit

    if (is.null(s$strata)) stop("can't index; object contains only one ECDF")

    s$n = s$strata.all[i]

    s = s[i]

    s$strata        = NULL
    s$strata.all    = NULL
    s$ntimes.strata = NULL

    new("cenfit", survfit=s)
})

setMethod("plot", signature(x="cenfit"),
          function(x, y, log  = 'x', axLimitFactor = 0.8, 
                   ylab = "Probability", xlab = "Value", 
                   lty  = seq(1,6), ...)
{
    s = x@survfit

    firstx = (min(s$time)*axLimitFactor)

    plot(s, log=log, firstx=firstx, ylab=ylab, xlab=xlab, lty=lty, ...)

    # Draw a vertical line at the minimum observation -- this is
    # the lower extent of the step curve or ECDF.
    abline(v=min(s$time))

    abline(h=1.0)
})

setMethod("lines", "cenfit", function(x, ...) lines(x@survfit, ...))

# Private predict method for cenfit objects -- does not handle strata!
.predict.cenfit =
function(x, newdata, conf.int=FALSE) 
{
    ret  = NULL
    s = x@survfit
    pred = stepfind(s$time, s$surv, newdata)

    if (!conf.int) ret = pred
    else
      {
        pred.l = stepfind(s$time, s$lower, newdata)
        pred.u = stepfind(s$time, s$upper, newdata)

        ret = data.frame(newdata, pred, pred.l, pred.u)

        names(ret) = c("obs", "prob", LCL(x), UCL(x))
      }

    return(ret)
}

# Public predict method for cenfit objects
setMethod("predict", signature(object="cenfit"),
          function(object, newdata, conf.int=FALSE) 
{
    ret  = NULL
    s = object@survfit

    if (is.null(s$strata)) ret = .predict.cenfit(object, newdata, conf.int)
    else
      {
        for (i in 1:length(s$strata))
          {
            ret[[i]] = .predict.cenfit(object[i], newdata, conf.int)
          }
        names(ret) = names(s$strata)
        class(ret) = "NADAlist"
      }
      
    return(ret)
})

# Private quantile method for cenfit objects -- does not handle strata!
.quantile.cenfit =
function(x, newdata, conf.int=FALSE) 
{
    ret  = NULL
    s = x@survfit
    quan = stepfind(s$surv, s$time, newdata, FALSE)

    if (!conf.int) 
      {
        ret = quan
        names(ret) = 
          paste(formatC(100 * newdata, format = "fg", wid = 1), "%", sep='')
      }
    else
      {
        quan.l = stepfind(s$surv, s$lower, newdata, FALSE)
        quan.u = stepfind(s$surv, s$upper, newdata, FALSE)

        ret = data.frame(newdata, quan, quan.l, quan.u)

        names(ret) = c("quantile", "obs", LCL(x), UCL(x))
      }

    return(ret)
}

# Public quantile method for cenfit objects
setMethod("quantile", signature(x="cenfit"),
          function(x, probs = NADAprobs, conf.int = FALSE, ...)
{
    ret = NULL
    s = x@survfit

    if (is.null(s$strata)) ret = .quantile.cenfit(x, probs, conf.int)
    else
      {
        for (i in 1:length(s$strata))
          {
            ret[[i]] = .quantile.cenfit(x[i], probs, conf.int)
          }
        names(ret) = names(s$strata)
        class(ret) = "NADAlist"
      }

    return(ret)
})

# Private mean method for cenfit objects -- does not handle strata!
# Given a cenfit or survfit object returns the calculated mean and se(mean)
.mean.cenfit =
function(x)
{
    x = x@survfit

    stime   = x$time
    surv    = x$surv
    n.risk  = x$n.risk
    n.event = x$n.event

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

    ret = c(sum(mean), sqrt(varmean))
    names(ret) = c("rmean", "se(rmean)")

    return(ret)
}

# Public mean method for cenfit objects
setMethod("mean", signature(x="cenfit"), function(x, ...)
{
    ret = NULL
    s = x@survfit

    if (is.null(s$strata)) ret = .mean.cenfit(x)
    else
      {
        for (i in 1:length(s$strata)) ret[[i]] = .mean.cenfit(x[i])
        names(ret) = names(s$strata)
        class(ret) = "NADAlist"
      }

    return(ret)
})

setMethod("median", signature(x="cenfit"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    ret = as.vector(quantile(x, 0.5))
    names(ret) = names(x@survfit$strata)
    return(ret)
})

setMethod("sd", signature(x="cenfit"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    ret = NULL
    s = x@survfit

    if (is.null(s$strata)) ret = (sqrt(s$n) * as.numeric(mean(x)[2]))
    else
      {
        for (i in 1:length(s$strata))
          {
            n = as.numeric(s$strata.all[i])
            ret[[i]] = (sqrt(n) * as.numeric(mean(x[i])[2]))
          }
        names(ret) = names(s$strata)
      }

    return(ret)
})

setMethod("print", signature(x="cenfit"), function(x, ...)
{
    s = x@survfit

    summaryVec =
    function(x)
    {
      s = x@survfit

      n      = s$n
      cen    = s$n - sum(s$n.event)
      median = median(x)
      mean   = mean(x)

      return(c(n, cen, median, mean))
    }

    ret = NULL
    tag = c("n", "n.cen", "median", "mean", "se(mean)")

    if (is.null(s$strata))
      {
        ret = summaryVec(x)
        names(ret) = tag
      }
    else
      {
        ret = summaryVec(x[1])
        for (i in 2:length(s$strata))
          {
            ret = rbind(ret, summaryVec(x[i]))
          }
        colnames(ret) = tag
        rownames(ret) = names(s$strata)
      }

    print(ret, ...)
    invisible(ret)
})

setMethod("summary", signature(object="cenfit"),
          function(object, times, censored=FALSE, scale=1, ...)
{
    strataSummary =
    function(x)
    {
      s = x@survfit

      # Reverse vectors -- we have leftist views ;)
      s = lapply(s, rev)
      ret = data.frame(s$time, s$n.risk, s$n.event, 
                       s$surv, s$std.err, s$lower, s$upper)
      names(ret) = c("obs", "n.risk", "n.event", 
                     "prob", "std.err", LCL(x), UCL(x))
      return(ret)
    }

    s   = object@survfit
    ret = NULL

    if (is.null(s$strata)) ret = strataSummary(object)
    else
      {
        for (i in 1:length(s$strata)) ret[[i]] = strataSummary(object[i])
        names(ret) = names(s$strata)
        class(ret) = "NADAlist"
      }

    return(ret)
})

#-->> END survival analysis based methods -- "cencen" section
