#-->> BEGIN Maximum Likelihood Estimation (MLE) Regression section

## Generics

setGeneric("cenmle",
  function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setOldClass("survreg")

setClass("cenmle", representation(survreg="survreg"))
setClass("cenmle-gaussian", representation("cenmle"))
setClass("cenmle-lognormal", representation("cenmle"))

## Core Methods

# This keeps things orthogonal with the survival package
cenreg = cenmle

setMethod("cenmle",
          signature(obs="formula", censored="missing", groups="missing"),
                    function(obs, censored, groups, dist, ...)
{
    dist = ifelse(missing(dist), "lognormal", dist)
    switch(dist,
        gaussian  = new_cenmle_gaussian(obs, dist, ...),
        lognormal = new_cenmle_lognormal(obs, dist, ...),
        survreg(asSurv(obs), dist=dist, ...)
    )
})

setMethod("cenmle", 
          signature(obs="Cen", censored="missing", groups="missing"), 
          cencen.Cen)

setMethod("cenmle",
          signature(obs="numeric", censored="logical", groups="missing"),
                    cencen.vectors)

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)

setMethod("summary", signature(object="cenmle"), function(object, ...)
{
    # To do: modify the call object to reflect the NADA call
    # for now, just nullify it -- users typically remember the call.
    object@survreg$call = NULL
    summary(object@survreg, ...)
})

setMethod("print", signature(x="cenmle"), function(x, ...)
{
    print(summary(x, ...))

    ## This is nice, but it doesn't work with strata
    #ret = c(mean(x), median(x), sd(x))
    #names(ret) = c("mean", "median", "sd")
    #print(ret)
    #invisible(ret)
})


setMethod("predict", signature(object="cenmle"), summary) 

#setMethod("predict", signature(object="cenmle"), 
#          function(object, newdata, conf.int=FALSE, ...)
#{
#    predict(object@survreg, newdata, ...)
#})

setMethod("residuals", signature(object="cenmle"), function(object, ...)
{
    residuals(object@survreg, ...)
})

setMethod("coef", signature(object="cenmle"), function(object, ...)
{
    coef(object@survreg, ...)
})

setMethod("median", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef))
})

setMethod("sd", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    ret = exp(2*x@survreg$coef + x@survreg$scale^2)*(exp(x@survreg$scale^2)-1)
    as.vector(sqrt(ret))
})

setMethod("mean", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef + 0.5*(x@survreg$scale)^2))
})

setMethod("median", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$coef)
})

setMethod("sd", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$scale)
})

setMethod("mean", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$coef)
})


## Supporting Functions 

# cenmle for lognormal distributions

new_cenmle_lognormal =
function(formula, dist, ...)
{
    new("cenmle-lognormal", survreg=survreg(asSurv(formula), dist=dist, ...))
}

# cenmle for gaussian, or normal, distributions

# If a normal distribution is assumed the input data must be expressed
# as an interval between zero and the DL.  They cannot simply be stated as
# 'left' censored, because that allows some probability of going below 0,
# and estimates will be biased low and wrong.  So with the normal option
# and left censoring, internally we must use interval censoring.  The end of
# the interval are the detected values.  The start of the interval will
# have identical numbers in it for the detects, and a 0 for the 
# nondetects (a simple trick is: start = obs - obs * censored).

new_cenmle_gaussian =
function(formula, dist, ...)
{
    obs      = formula[[2]][[2]]
    censored = formula[[2]][[3]]
    groups   = formula[[3]]

    f = as.formula(substitute(Surv(o - o * c, o, type="interval2")~g, 
                                   list(o=obs, c=censored, g=groups)))
    environment(f) = parent.frame()
    new("cenmle-gaussian", survreg=survreg(f, dist=dist, ...))
}

#-->> END Regression on Maximum Likelihood Estimation (MLE) section

