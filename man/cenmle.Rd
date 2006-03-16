\name{cenmle}
\alias{cenmle}
\alias{cenreg}
\title{
    Regression by Maximum Likelihood Estimation for Left-censored Data
}
\description{

     Regression by Maximum Likelihood (ML) Estimation for left-censored
     ("nondetect" or "less-than") data.  This routine computes regression
     estimates of slope(s) and intercept by maximum likelihood
     when data are left-censored.  It will compute ML estimates of
     descriptive statistics when explanatory variables following the ~
     are left blank.  It will compute ML tests similar in function and
     assumptions to two-sample t-tests and analysis of variance when
     groups are specified following the ~.  It will compute regression
     equations, including multiple regression, when continuos explanatory
     variables are included following the ~.  It will compute the ML
     equivalent of analysis of covariance when both group and continuous
     explanatory variables are specified following the ~.  To avoid
     an appreciable loss of power with regression and group hypothesis
     tests, a probability plot of residuals should be checked to ensure
     that residuals from the regression model are approximately gaussian.
}
\synopsis{cenmle(obs, censored, groups, ...)}
\usage{
    cenmle(obs, censored, groups, ...)
}
\arguments{
\item{obs}{
    Either a numeric vector of observations or a formula.
    See examples below.
} 
\item{censored}{
    A logical vector indicating TRUE where an observation in `obs' is
    censored (a less-than value) and FALSE otherwise.
}
\item{groups}{
    A factor vector used for grouping `obs' into subsets.
}
\item{...}{
    Additional items that are common to this function and the \code{survreg}
    function from the `survival' package.  The most important of which
    is `dist' (distribution to use).  See Details below.
}
} %arguments
\details{
    This routine is a front end to the \code{survreg} routine in the
    \code{survival} package.  
    
    There are many additional options that are supported and documented
    in \code{survfit}.  Only a few have application to the geosciences.

    The most important option is `dist' which specifies the distributional
    model to use in the regression.  The default is `lognormal'.  

    Also supported is a `gaussian' or a normal distribution.  The use of
    a gaussian distribution requires an interval censoring context for
    left-censored data.  Luckily, this routine automatically does this
    for you -- simply specify `gaussian' and the correct manipulations
    are done.

    If any other distribution is specified besides lognormal or gaussian, 
    the return object is a raw survreg object -- it is up to the user to
    `do the right thing' with the output (and input for that matter).

    If you are using the formula interface:  The \code{censored} and
    \code{groups} parameters are not specified -- all information is provided
    via a formula as the \code{obs} parameter.  The formula must have a
    \code{Cen} object as the response on the left of the \code{~} operator and,
    if desired, terms separated by + operators on the right.  
    See Examples below.  
} %details
\value{
    a \code{cenmle} object.
    Methods defined for \code{cenmle} objects are provided for
    \code{mean}, \code{median}, \code{sd}.
}
\references{
Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
censored environmental data.  John Wiley and Sons, USA, NJ.
}
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
%    \code{\link{cenmle-methods}},
    \code{\link{Cen}}, 
    \code{\link{mean-methods}}, 
    \code{\link{sd-methods}}, 
    \code{\link{median-methods}}, 
%    \code{\link{quantile-methods}}, 
%    \code{\link{predict-methods}}, 
    \code{\link{summary-methods}} 
}

\examples{

    # Create a MLE regression object 

    data(TCEReg)

    obs      = TCEReg$TCEConc 
    censored = TCEReg$TCECen

    mycenmle = cenmle(obs, censored) 

    summary(mycenmle)
    median(mycenmle)
    mean(mycenmle)
    sd(mycenmle)

    # With groupings
    groups = as.factor(TCEReg$LandUse)
    cenmle(obs, censored, groups)
    
    # Formula interface -- no groups
    cenmle(Cen(obs, censored)) 
    # assume a gaussian distribution
    #cenmle(Cen(obs, censored), dist="gaussian") 

    # Formula interface -- a complicated regression
    attach(TCEReg)
    cenmle(Cen(TCEConc, TCECen)~PopDensity+PctIndLU+Depth)
}
\keyword{survival}
