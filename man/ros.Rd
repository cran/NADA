\name{ros} 
\alias{ros} 
\alias{lros} 
\alias{print.ros}
\alias{as.data.frame.ros}
\title{Regression on Order Statistics} 
\description{
    \code{ros} is an implementation of a Regression on Order Statistics
    (ROS) designed for multiply censored analytical chemistry data.

    The method assumes data contains zero to many left censored (less-than)
    values.

    By default, \code{ros} performs a log transformation prior to, 
    and after operations over the data.
} 
\usage{ ros(obs, censored, forwardT="log", reverseT="exp")} 
\arguments{
        \item{obs}{
            A numeric vector of observations.  This includes both censored
	        and uncensored observations.
        } 
        \item{censored}{
            A logical vector indicating TRUE where an observation in v is
	        censored (a less-than value) and FALSE otherwise.
        } 
        \item{forwardT}{
	        A name of a function to use for transformation prior to performing
	        the ROS fit.  Defaults to \code{log}.
        } 
        \item{reverseT}{
	        A name of a function to use for reversing the transformation
	        after performing the ROS fit. Defaults to \code{exp}.
        }

} 
\details{

    The proceedure first computes the Weibull-type plotting positions of
    the combined uncensored and censored observations using a formula
    designed for multiply-censored data (see \code{hc.ppoints}).
    A linear regression is formed using the plotting positions of the
    uncensored observations and their normal quantiles.  This model is
    then used to estimate the concentration of the censored observations
    as a function of their normal quantiles.  Finally, the observed
    uncensored values are combined with modeled censored values to
    corporately estimate summary statistics of the entire population. By
    combining the uncensored values with modeled censored values, this
    method is more resistant of any non-normality of errors, and reduces
    any transformation errors that may be incurred.

}

%% Function return values 
\value{
    \code{ros} returns an object of class c("ros", "lm").

    \code{print} displays a simple summary of the ROS model.
    \code{as.data.frame} converts the modeled data in a ROS model to
    a data frame. Note that this discards all linear-model information
    from the object.
} 
\references{

    Lee and Helsel (in press) Statistical analysis of environmental data
    containing multiple detection limits: S-language software for linear
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
    \code{\link{splitQual}}, \code{\link{predict.ros}}, 
    \code{\link{plot.ros}}, \code{\link{summary.ros}}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    myros = ros(obs, censored) 

    summary(myros)
    as.data.frame(myros)
} 
\keyword{regression} 
