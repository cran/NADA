\name{summary.ros} 
\alias{summary.ros} 
\alias{plot.ros} 
\alias{abline} 
\alias{abline.default} 
\alias{abline.ros} 
\title{Summary and Plotting of Regression on Order Statistics models} 
\description{
    Generic implementations of summary and plotting functions for ROS models.
    
    \code{summary.ros} provides the standard summary of the linear model.
    Has the option of producing the standard linear-model quality plots.

    \code{plot.ros} constructs a probability plot representation of
    the ROS model.

    \code{abline.ros} adds a line representing the linear regression model
    to an existing ROS prob-plot -- For comparing two or more ROS models.
} 
\usage{
    \method{summary}{ros}(object, plot=FALSE, ...) 
    \method{plot}{ros}(x, plot.censored=FALSE, lm.line=TRUE, grid=TRUE, ylab="Value", pch=16, ...)
    \method{abline}{ros}(a, ...)
}
\arguments{

%% summary.ros arguments
    \item{object}{
        A ROS model constructed using \code{\link{ros}}
    }
    \item{plot}{
        Logical indicating whether the standard \code{\link{lm}} quality plots
        should be drawn.
    }

%% plot.ros arguments
    \item{x}{
        A ROS model constructed using \code{\link{ros}}
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
        A ROS model constructed using \code{\link{ros}}
    }
 
}
\details{
}

%% Function return values 
%\value{
%} 
\references{

    Lee and Helsel (in press) Statistical analysis of environmental data
    containing multiple detection limits: S-language software for 
    regression on order statistics, Computers in Geoscience vol. X, pp. X-X

    Lee and Helsel (in press) Baseline models of trace elements in
    major aquifers of the United States. Applied Geochemistry vol. X, pp. X-X.

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
    \code{\link{ros}},
    \code{\link{predict.ros}},
    \code{\link{quantile.ros}},
    \code{\link{sd.ros}},
    \code{\link{mean.ros}},
    \code{\link{median.ros}}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    myros = ros(obs, censored) 
    summary(myros, plot=TRUE) # Cycles through standard lm() plots

    plot(myros, plot.cen=TRUE) # Forms a prob-plot of the ROS model

    obs2      = c(0.5,    1.0,  1.5,   1.5,    10,   100,   200)
    censored2 = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

    abline(ros(obs2, censored2)) # Draws this regression line on existing plot
} 
\keyword{regression} 
