\name{quantile-methods}
\docType{methods}
\alias{quantile-methods}
\alias{quantile,ANY-method}
\alias{quantile,cenfit-method}
\alias{quantile,ros-method}
\alias{NADAprobs}
\title{Methods for function quantile in Package NADA}
\description{
  Methods for the function \code{quantile} in package \pkg{NADA} 

  Compute the modeled values of quantiles or probabilities using a model object.
}
\usage{

\S4method{quantile}{ros}(object, probs=NADAprobs, ...)

\S4method{quantile}{cenfit}(object, probs=NADAprobs, conf.int=FALSE, ...)
} % End usage
\arguments{

\item{object}{
    An output object from a NADA fuction such as \code{\link{ros}}.
}
\item{probs}{
    Numeric vector of probabilites for which to calculate model values.
    The default is the global variable 
    NADAprobs = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95).
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
\examples{
    
    data(Cadmium)

    mymodel = cenfit(Cadmium$Cd, Cadmium$CdCen, Cadmium$Region)

    quantile(mymodel, conf.int=TRUE)
}
\keyword{methods}
