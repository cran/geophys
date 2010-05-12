\name{randFRY}
\Rdversion{1.1}
\alias{randFRY}
\title{Random data for Fry analysis
}
\description{Random distribution of
  points that have a specified radial spacing.
}
\usage{
randFRY(N = 200, LIM = c(0, 0, 100, 100), rlen = 5, PLOT=FALSE)
}
\arguments{
  \item{N}{integer, number of points
}
  \item{LIM}{Limit in x-y direction
}
\item{rlen}{distance length.
  points are selected if they are greater than this distance
  from each other.
}
 \item{PLOT}{logical, TRUE=make a plot of data }


}
\details{
  A spatial data set is generated
  that has random distribution
  with points no closer than a specified
  distance apart.
}
\value{
  \item{x}{x-coordinate}
  \item{y}{y-coordinate}
}
\references{

  Fry, N., (1979)
Random point distributions and strain measurement in rocks
Tectonophysics, 60:89-105.
 
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{dofry
}
\examples{

RDAT = randFRY(100, LIM=c(0,0, 100, 100) , rlen=5   )
length(RDAT$x)
plot(RDAT$x, RDAT$y, asp=1, pch=".", cex=2)


}
\keyword{misc}
