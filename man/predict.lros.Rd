\name{predict.lros} 

\alias{predict.lros} 
\alias{quantile.lros} 

\alias{mean.lros}
\alias{median}
\alias{median.default}
\alias{median.lros}
\alias{sd}
\alias{sd.default}
\alias{sd.lros}

\title{Query and Prediction with LROS models} 
\description{
    Functions that perform query and/or prediction with Linear Regression
    on Order Statistics (\code{\link{lros}}) objects.

    \code{mean} returns the modeled mean of a LROS model.
    \code{median} returns the modeled median of a LROS model.
    \code{sd} returns the modeled standard deviation of a LROS model.
 
    \code{quantile} produces sample quantiles corresponding to
    the given probabilities. The smallest observation corresponds to a
    probability of 0 and the largest to a probability of 1.

    \code{predict} Predicted values based on a LROS model object.
} 
\usage{
    \method{mean}{lros}(x, ...)
    \method{median}{lros}(x, na.rm=FALSE)
    \method{sd}{lros}(x, na.rm=FALSE)
    \method{quantile}{lros}(x, probs=c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95), ...)
    \method{predict}{lros}(object, newdata, ...)
}
\arguments{
    \item{x, object}{
        A LROS model constructed using \code{\link{lros}}
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
    \item{...}{
        Additional arguments passed to the generic method.
    }
}
\details{
}

\value{
    All functions return a numeric vector of results.
} 
\references{

    Lee and Helsel (in press) Statistical analysis of environmental data
    containing multiple detection limits: S-language software for linear
    regression on order statistics, Computers in Geoscience vol. X, pp. X-X

    Lee and Helsel (in press) Baseline models of trace elements in
    drinking water of the United States. Applied Geochemistry vol. X, pp. X-X.

    Dennis R. Helsel (2004), Nondetects And Data Analasis:
    John Wiley and Sons, New York.

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
    \code{\link{lros}}
    \code{\link{hc.ppoints}}
}
\note{
    Some of these fuctions mask the original functions in the base
    package.  For the documentation on the original functions use:
    \code{help("foo", package="base")}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    mylros = lros(obs, censored) 

    mean(mylros)
    median(mylros) 
    sd(mylros)
 
    quantile(mylros, probs=c(0.90, 0.95))
    predict(mylros, 1.5)
} 
\keyword{regression} 
