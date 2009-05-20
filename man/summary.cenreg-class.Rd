\name{summary.cenreg-class}
\docType{class}
\alias{summary.cenreg-class}
\alias{summary.cenreg}

\title{Class "summary.cenreg"}
\description{
    A "summary.cenreg" object is returned from \code{summary}.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form 
    \code{summary(cenreg(obs, censored, groups))}.
}
\section{Slots}{
  \describe{
    No slots, it just extends a list.
    %\item{\code{survreg}:}{Object of class \code{"survreg"}  }
  }
}
\section{Extends}{
    Class \code{"list"}.
    Class \code{"vector"}, by class \code{"list"}.
}
\section{Methods}{
  \describe{
    %\item{summary}{\code{signature(object = "cenreg")}: ... }
  }
}
\author{
    Lopaka Lee <rclee@usgs.gov>

    Dennis Helsel <dhelsel@usgs.gov>
} 
\seealso{
  \code{\link{cenreg}}
}
\keyword{classes}
