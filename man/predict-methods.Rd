\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,cenfit-method}
\alias{predict,ros-method}
\title{Methods for function predict in package NADA}
\description{
  Functions that perform predictions using NADA model objects.

  For \code{ros} models, predict the normal quantile of a value.

  For \code{cenfit} objects, predict the probabilities of new observations.
}
\usage{
\S4method{predict}{ros}(object, newdata, ...)

\S4method{predict}{cenfit}(object, newdata, conf.int=FALSE, ...)
} % End usage
\arguments{

\item{object}{
    An output object from a NADA function such as \code{\link{ros}}.
}
\item{newdata}{
        Numeric vector of data for which to predict model values.
        For \code{ros} objects this will be new normalized quantiles 
        of plotting positions.  For \code{cenfit} objects this will be
        new observations for which you desire the modeled probabilities.
}
\item{conf.int}{
        A logical indicating if confidence intervals should be computed.  For
        \code{cenfit} objects, the confidence interval is set during the call
        to \code{\link{cenfit}}.  Currently not supported for \code{ros} 
        objects.
}
\item{...}{
        Additional arguments passed to the generic method.
}
} % End arguments
\keyword{methods}
