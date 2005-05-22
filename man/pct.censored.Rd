\name{pct.censored} 
\alias{pct.censored} 
\title{Calculate the percentage of values censored} 
\description{
    \code{pct.censored} is a simple, but convenient, function that calculates
    the percentage of censored values.  
} 
\usage{ pct.censored(obs, censored)} 
\arguments{
        \item{obs}{
            A numeric vector of observations.  This includes both censored
	        and uncensored observations.
        } 
        \item{censored}{
            A logical vector indicating TRUE where an observation in v is
	        censored (a less-than value) and FALSE otherwise.
        } 
} 
\details{
    100*(length(obs[censored])/length(obs))
}

%% Function return values 
\value{
    \code{pct.censored} returns a single numeric value representing
    the percentage of values censored in the ``obs" vector.
} 
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
    \code{\link{splitQual}}, \code{\link{ros}}, 
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    pct.censored(obs, censored) 
} 
\keyword{utilities} 
