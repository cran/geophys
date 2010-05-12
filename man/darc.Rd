\name{darc}
\alias{darc}
\title{Circular Arc
}
\description{Draw acircular arc
  from angle 1 to angle 2 at a given location.
}
\usage{
darc(rad = 1, ang1 = 0, ang2 = 360, x1 = 0, y1 = 0, n = 1)
}
\arguments{
  \item{rad}{radius
}
  \item{ang1}{angle 1, degrees
}
  \item{ang2}{angle 2, degrees
}
  \item{x1}{x location, plot coordinates
}
  \item{y1}{y location, plot coordinates
}
  \item{n}{increment for number of segments, degrees
}
}
\details{If angle1 > angle2 arc is drawn in opposite direction
}
\value{list(x,y)
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{

plot(c(0,1), c(0,1), type='n', ann=FALSE, asp=1)

A = darc(.3, 23, 47, .5, .5, n=1)

lines(A$x, A$y)



}
\keyword{misc}
