\name{DGzx}
\alias{DGzx}
\title{Gravity anomaly in 2.5D
}
\description{Gravity anomaly in 2.5-Dimensions from
  an arbitrary polynomial at many stations.
}
\usage{
DGzx(xs, zs, xv, zv, den)
}
\arguments{
  \item{xs}{station locations in X
}
  \item{zs}{station locations in Z
}
  \item{xv}{x-vertices
}
  \item{zv}{z-vertices
}
  \item{den}{density contrast
}
}
\details{calculate the
  2.5D solution to gravity.
  Orientation of the vertices
  should be right handed.
}
\value{vector of Delta-Gz and Delta-Gx at each station
}
\references{
  Won and Bevis
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\examples{

nstn = 10
xstart = -10000
xend = 10000
xcen = 0
zcen = 5000
RAD = 2000


xs = seq(from=xstart, by=(xend-xstart)/nstn , length=nstn)
zs = rep(0, length=length(xs))

den = 0.2

Np = 6
theta = seq(from=0, to=2*pi, length=Np)
KZ = list(x=NA, y=NA)
KZ$x = xcen+RAD*cos(theta)
KZ$y = zcen+RAD*sin(theta)




Ngrav = DGzx(xs, zs,  KZ$x, KZ$y, den)




}
\keyword{misc}
