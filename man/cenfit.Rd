\name{cenfit}
\alias{cenfit}
\alias{print.cenfit}
\alias{summary.cenfit}
\alias{print.summary.cenfit}
\alias{plot.cenfit}
\alias{print.cenfit}
\alias{[.cenfit}
%\alias{cenfit.km}
%\alias{cenfit.coxph.null}
%\alias{cenfit.coxph}
%\alias{basehaz}
\title{
    Compute an ECDF for Censored Data
}
\description{
    Computes an estimate of an empirical cumulative distribution function
    (ECDF) for censored data using the Kaplan-Meier method.
}
\synopsis{cenfit(formula, data=NULL, ...)}
\usage{
%    survfit(formula, data, weights, subset, na.action, 
%            newdata, individual=F, conf.int=.95, se.fit=T, 
%            type=c("kaplan-meier","fleming-harrington", "fh2"),
%            error=c("greenwood","tsiatis"),
%            conf.type=c("log","log-log","plain","none"),
%            conf.lower=c("usual", "peto", "modified"))

cenfit(formula, 
       na.action, 
       se.fit=T, 
       conf.int=.95, 
       conf.type=c("log","log-log","plain","none"),
       conf.lower=c("usual", "peto", "modified"))

%\method{[}{cenfit}(x,...,drop=FALSE)
%basehaz(fit,centered=TRUE)

}
\arguments{
\item{formula}{
    The formula must have a \code{Cen} object as the response on the
    left of the \code{~} operator and, if desired, terms separated by +
    operators on the right.  One of the terms may be a \code{strata}
    object.  For a single survival curve the \code{"~ 1"} part of the
    formula is not required.
}
%\item{data}{
%    a data frame in which to interpret the variables named in the formula,
%    or in the \code{subset} and the \code{weights} argument.
%}
%\item{weights}{
%    The weights must be nonnegative and it is strongly recommended that 
%    they be strictly positive, since zero weights are ambiguous, compared
%    to use of the \code{subset} argument.
%}
%\item{subset}{
%    expression saying that only a subset of the rows of the data
%    should be used in the fit.
%}
\item{na.action}{
    a missing-data filter function, applied to the model frame.
%    , after any \code{subset} argument has been used.
    Default is \code{options()$na.action}.
}
%\item{newdata}{
%    a data frame with the same variable names as those that appear in
%    the \code{coxph} formula.  Only applicable when \code{formula} is a
%    \code{coxph} object.  The curve(s) produced will be representative of
%    a cohort who's covariates correspond to the values in \code{newdata}.
%    Default is the mean of the covariates used in the \code{coxph} fit.
%}
%\item{individual}{
%    a logical value indicating whether the data frame represents different
%    time epochs for only one individual (T), or whether multiple rows
%    indicate multiple individuals (F, the default).  If the former only
%    one curve will be produced; if the latter there will be one curve
%    per row in \code{newdata}.
%}
\item{se.fit}{
    a logical value indicating whether standard errors should be computed.
    Default is \code{TRUE}.
}
\item{conf.int}{
    the level for a two-sided confidence interval on  the survival
    curve(s).  Default is 0.95.
}
%\item{type}{
%    a character string specifying the type of empirical cumulative
%    distribution function to generate.  Possible values are
%    \code{"kaplan-meier"}, \code{"fleming-harrington"} or \code{"fh2"}
%    if a formula is given and \code{"aalen"} or \code{"kaplan-meier"}
%    if the first argument is a \code{coxph} object, (only the first
%    two characters are necessary).  The default is \code{"aalen"} when
%    a \code{coxph} object is given, and it is \code{"kaplan-meier"}
%    otherwise.
%}
%\item{error}{
%    either the string \code{"greenwood"} for the Greenwood formula or
%    \code{"tsiatis"} for the Tsiatis formula, (only the first character
%    is necessary).  The default is \code{"tsiatis"} when a \code{coxph}
%    object is given, and it is \code{"greenwood"} otherwise.
%}
\item{conf.type}{
    One of \code{"none"}, \code{"plain"}, \code{"log"} (the default),
    or \code{"log-log"}.  Only enough of the string to uniquely identify
    it is necessary.  The first option causes confidence intervals not to
    be generated.  The second causes the standard intervals \code{curve
    +- k *se(curve)}, where k is determined from \code{conf.int}.
    The log option calculates intervals based on the cumulative hazard
    or log(probability). The last option bases intervals on the log hazard
    or log(-log(probability)).  These last will never extend past 0 or 1.
}
\item{conf.lower}{
    controls modified lower limits to the curve, the upper limit remains
    unchanged.  The modified lower limit is based on an 'effective n'
    argument.  The confidence bands will agree with the usual calculation
    at each detected value, but unlike the usual bands the
    confidence interval becomes wider at each censored observation.
    The extra width is obtained by multiplying the usual variance by
    a factor m/n, where n is the number currently at risk and m is the
    number at risk at the last censored observation.  (The bands thus
    agree with the un-modified bands at each detected observation.)
    This is especially useful for ECDF curves with a long flat tail.

    The Peto lower limit is based on the same 'effective n' argument as
    the modified limit, but also replaces the usual Greenwood variance
    term with a simple approximation.  It is known to be conservative.
}
%\item{x}{a \code{cenfit} object}
%\item{fit}{a \code{coxph} object}
%\item{centered}{
%    Compute the baseline hazard at the covariate mean rather than at zero?
%}
} %%END arguments
\value{
    a \code{cenfit} object.
    Methods defined for \code{cenfit} objects are provided for
    \code{print}, \code{plot}, \code{lines}, \code{predict},
    \code{mean}, \code{median}, \code{sd}, \code{quantile}.

%    For \code{basehaz}, a dataframe with the baseline hazard, times,
%    and strata.

    If the input formula contained factoring groups 
    (ie., \code{cenfit(Cen(obs, censored)~groups))}, individual ECDFs can be
    obtained by indexing (eg., \code{model[1]}, etc.).
}
\details{

%Actually, the estimates used are the Kalbfleisch-Prentice (Kalbfleisch
%and Prentice, 1980, p.86) and the Tsiatis/Link/Breslow, which reduce to
%the Kaplan-Meier and Fleming-Harrington estimates, respectively, when
%the weights are unity.  When curves are fit for a Cox model, subject
%weights of \code{exp(sum(coef*(x-center)))} are used, ignoring any value
%for \code{weights} input by the user.  There is also an extra term in
%the variance of the curve, due to the variance ofthe coefficients and
%hence variance in the computed weights.


%The Greenwood formula for the variance is a sum of terms d/(n*(n-m)),
%where d is the number of deaths at a given time point, n is the sum
%of \code{weights} for all individuals still at risk at that time,
%and m is the sum of \code{weights} for the deaths at that time.
%The justification is based on a binomial argument when weights are all
%equal to one; extension to the weighted case is ad hoc.  Tsiatis (1981)
%proposes a sum of terms d/(n*n), based on a counting process argument
%which includes the weighted case.


%The two variants of the F-H estimate have to do with how ties are handled.
%If there were 3 deaths out of 10 at risk, then the first would increment
%the hazard by 3/10 and the second by 1/10 + 1/9 + 1/8.  For curves
%created after a Cox model these correspond to the Breslow and Efron
%estimates, respectively, and the proper choice is made automatically.
%The \code{fh2} method will give results closer to the Kaplan-Meier.


%Based on the work of Link (1984), the log transform is expected to produce
%the most accurate confidence intervals.  If there is heavy censoring,
%then based on the work of Dorey and Korn (1987) the modified estimate
%will give a more reliable confidence band for the tails of the curve.

This, and related routines, are front ends to routines in the
\code{survivial} package.  Since the survival routines can not handle
left-censored data, these routines tranparently handle ``flipping" input
data and resultant calculations.  Additionally provided are query and
prediction methods for \code{cenfit} objects (see below).

The ECDF estimates used are the Kalbfleisch-Prentice variety (Kalbfleisch
and Prentice, 1980, p.86), which reduces to the Kaplan-Meier, when
weights are unity.

The Greenwood formula for the variance is used to calculate errors. The
Greenwood formula uses a sum of terms d/(n*(n-m)), where d is the number
of detections at a given concentration, n is the sum of 'weights' for
all individuals at or below that concentration,, and m is the sum of
'weights' for the detected obs at that time.  The justification is based
on a binomial argument when weights are all equal to one; extension to the
weighted case is ad hoc.  Tsiatis (1981) proposes a sum of terms d/(n*n),
based on a counting process argument which includes the weighted case.
}

\references{
Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
censored environmental data.  John Wiley and Sons, USA, NJ.


Dorey, F. J. and Korn, E. L. (1987).  Effective sample sizes for
confidence intervals for survival probabilities.  \emph{Statistics in
Medicine} 6, 679-87.


Fleming, T. H. and Harrington, D.P. (1984).  Nonparametric estimation of
the survival distribution in censored data.  \emph{Comm. in Statistics}
13, 2469-86.


Kalbfleisch, J. D. and Prentice, R. L. (1980).  \emph{The Statistical
Analysis of Failure Time Data.} Wiley, New York.


Link, C. L. (1984). Confidence intervals for the survival function using
Cox's proportional hazards model with covariates.  \emph{Biometrics}
40, 601-610.


Tsiatis, A. (1981). A large sample study of the estimate for the
integrated hazard function in Cox's regression model for survival
data. \emph{Annals of Statistics} 9, 93-108.
}
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
    \code{\link{Cen}}, 
    \code{\link{predict.cenfit}}, 
    \code{\link{quantile.cenfit}}, 
    \code{\link{median.cenfit}}, 
    \code{\link{mean.cenfit}}, 
    \code{\link{sd.cenfit}}, 
    \code{\link{cendiff}}, 
    \code{\link{survfit}}
}

\examples{

    # fit a Kaplan-Meier ECDF, plot and summarize it.
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    mycenfit = cenfit(Cen(obs, censored)) 

    plot(mycenfit)
    summary(mycenfit)
}
\keyword{survival}
