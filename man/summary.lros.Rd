\name{summary.lros} 
\alias{summary.lros} 
\alias{plot.lros} 
\alias{abline.lros} 
\title{Summary and Plotting of Linear Regression on Order Statistics models} 
\description{
    Generic implementations of summary and plotting functions for LROS models.
    
    \code{summary.lros} provides the standard summary of the linear model.
    Has the option of producing the standard linear-model quality plots.

    \code{plot.lros} constructs a probability plot representation of
    the LROS model.

    \code{abline.lros} adds a line representing the linear regression model
    to an existing LROS prob-plot -- For comparing two or more LROS models.
} 
\usage{
    \method{summary}{lros}(object, plot=FALSE, ...) 
    \method{plot}{lros}(x, plot.censored=FALSE, lm.line=TRUE, grid=TRUE, ylab="Value", pch=16, ...)
    \method{abline}{lros}(a, ...)
}
\arguments{

%% summary.lros arguments
    \item{object}{
        A LROS model constructed using \code{\link{lros}}
    }
    \item{plot}{
        Logical indicating whether the standard \code{\link{lm}} quality plots
        should be drawn.
    }

%% plot.lros arguments
    \item{x}{
        A LROS model constructed using \code{\link{lros}}
    }
    \item{plot.censored}{
        A logical indicating if the modeled censored values should
        be plotted.
    }
    \item{lm.line}{
        A logical indicating if the a line (abline) representing the
        linear regression model should be drawn on the plot.
    }
    \item{grid}{
        A logical indicating if dotted lines are to be drawn
        at major axis division.
    }
    \item{ylab}{
        A character string for the y-axis lable.
    }
    \item{pch}{
        The numerical code of the plotting symbols for the uncensored data.
        See \code{\link{par}}
    }

    \item{...}{
        Additional arguments passed to the generic method.
    }

    \item{a}{
        A LROS model constructed using \code{\link{lros}}
    }
 
}
\details{
}

%% Function return values 
%\value{
%} 
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
    \code{\link{lros}},
    \code{\link{predict.lros}},
    \code{\link{quantile.lros}},
    \code{\link{sd.lros}},
    \code{\link{mean.lros}},
    \code{\link{median.lros}}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    mylros = lros(obs, censored) 
    summary(mylros, plot=TRUE) # Cycles through standard lm() plots

    plot(mylros, plot.cen=TRUE) # Forms a prob-plot of the LROS model

    obs2      = c(0.5,    1.0,  1.5,   1.5,    10,   100,   200)
    censored2 = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

    abline(lros(obs2, censored2)) # Draws this regression line on existing plot
} 
\keyword{regression} 
