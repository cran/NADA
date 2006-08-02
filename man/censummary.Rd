\name{censummary}
\alias{censummary}
\title{
    Produces basic summary statistics on censored data
}
\description{
    Produces basic, and hopefully useful, summary statistics on censored data.
}
\synopsis{censummary(obs, censored)}
\usage{
    censummary(obs, censored)
}
\arguments{
\item{obs}{
    A numeric vector of observations.
} 
\item{censored}{
    A  logical vector indicating TRUE where an observation in x is
    censored (a less-than value) and FALSE otherwise.  
}
} % arguments
\value{
    A list with the following elements ...
}
\references{
    Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
    censored environmental data.  John Wiley and Sons, USA, NJ.
}
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>

    Dennis Helsel <dhelsel@usgs.gov>
} 
%\seealso{
%    \code{\link{NADA}},
%}
\examples{
    data(DFe)
    with(DFe, censummary(Summer, SummerCen))
}
\keyword{hplot}
