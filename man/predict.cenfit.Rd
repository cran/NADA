\name{predict.cenfit} 

\alias{predict.cenfit} 
\alias{quantile.cenfit} 
\alias{mean.cenfit}
\alias{median.cenfit}
\alias{sd.cenfit}

\title{Query and Prediction with ECDFs}
\description{
    Functions that perform query and/or prediction with empirical
    cumulative distribution function objects (ECDF).  ECDF objects are
    constructed using \code{\link{cenfit}}.

    \code{mean} returns the modeled mean of an ECDF.
    \code{median} returns the modeled median of an ECDF.
    \code{sd} returns the modeled standard deviation of an ECDF.
 
    \code{quantile} produces sample quantiles corresponding to the given
    probabilities using an ECDF. The smallest observation corresponds
    to a probability of 0 and the largest to a probability of 1.

    \code{predict} Predicted values based on an ECDF.
} 
\usage{
    \method{mean}{cenfit}(x, ...)
    \method{median}{cenfit}(x, na.rm=FALSE)
    \method{sd}{cenfit}(x, na.rm=FALSE)
    \method{quantile}{cenfit}(x, 
            probs=c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95),
            do.conf=FALSE, ...)
    \method{predict}{cenfit}(object, newdata, do.conf=FALSE, ...)
}
\arguments{
    \item{x, object}{
        A ECDF object constructed using \code{\link{cenfit}}
    }
    \item{na.rm}{
        logical. Should missing values be removed?
    }
    \item{probs}{
        Numeric vector with values in [0,1] -- the quantiles to predict.
    }
    \item{newdata}{
        Numeric vector of normalized quantiles of plotting positions.
    }
    \item{do.conf}{
        A logical indicating if confidence intervals should be computed.
    }
    \item{...}{
        Additional arguments passed to the generic method.
    }

}
\details{
    Note that computed confidence intervals are set using \code{cenfit}.
    Thus, the confidence intervals returned in functions with the
    \code{do.conf} option are fixed at the value originally passed to
    \code{cenfit}.
}

\value{
    All functions return a numeric vector of results by default.
    Functions with \code{do.conf} option will return a data frame if
    this option is true.
} 
\references{
    Dennis R. Helsel (2005).  Nondectects and Data Analysis; Statistics for
    censored environmental data.  John Wiley and Sons, USA, NJ.

    Dennis R. Helsel (1990), Less Than Obvious: Statistical Methods for,
    Environmental Science and Technology,
    vol.24, no. 12, pp. 1767-1774

    Dennis R. Helsel and Timothy A. Cohn (1988), Estimation of descriptive
    statistics for multiply censored water quality data, Water Resources
    Research vol. 24, no. 12, pp.1997-2004 
} 
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
    \code{\link{cenfit}}
}
\note{
    Some of these fuctions mask the original functions in the base
    package.  For the documentation on the original functions use:
    \code{help("foo", package="base")}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    mycenfit = cenfit(Cen(obs, censored)) 

    mean(mycenfit)
    median(mycenfit) 
    sd(mycenfit)
 
    quantile(mycenfit, probs=c(0.90, 0.95))
    predict(mycenfit, 1.5, do.conf=TRUE)
} 
\keyword{survival} 
