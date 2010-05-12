\name{randpoles}
\alias{randpoles}
\title{Random Poles
}
\description{Create a set of random poles distributed
  around a given pole.
}
\usage{
randpoles(az, iang, alphadeg, opt = "unif", BALL.radius = 1, N = 10, add = TRUE, ...)
}
\arguments{
  \item{az}{azimuth of given direction
}
  \item{iang}{inclination angle of given direction
}
  \item{alphadeg}{standard error for distribution
}
  \item{opt}{Optional distribution, "unif" or "norm"
}
  \item{BALL.radius}{radius of small circle to draw
}
  \item{N}{Number of points to create
}
  \item{add}{logical, TRUE=add points to net plot
}
  \item{\dots}{graphical parameters for plotting
}
}
\details{A given geographic pole is supplied and
  N normally or uniformly distributed poles are
  created.
}
\value{
  \item{az}{azimuths of new poles}
  \item{dip}{dips of new poles}
  \item{x}{x position on plot}
  \item{y}{y position on plot}
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{net
}
\examples{


net();
GG = randpoles(30, 40, 10, opt="norm", N=20)

addsmallcirc(30, 40, 10)




}
\keyword{misc}
