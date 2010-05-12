\name{get.heat2}
\alias{get.heat2}
\alias{ get.heat}
\title{Heat Equation Solution
}
\description{Solution of the diffusion equyation when
  contact is made instananeously with an infinite half space.
  Used for sudden exposure of magma at the surface, as in an oceanic ridge.
  
}
\usage{
get.heat2(x, T0, k, t)
}
\arguments{
  \item{x}{numeric vector, distance
}
  \item{T0}{Temperature at infinity
}
  \item{k}{diffusion constant
}
  \item{t}{time
}
}
\details{
  get.heat2 is used when temperature at the boundary (atmosphere)
  remains constant for all time.
  get.heat is used when temperature in half space changes (warms up).
  
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{
 k = 1*10^(-6)
dt = 3600
 dz = 20*10^(-2)

T0=25

T1 = 1200



x =  seq(from=0, to=80, by=0.5)/100

i = 1
Tx =  get.heat2(x, T1-T0, k, i*dt)

plot(Tx, max(x)-x, type='n', xlim=c(700, 1200) , axes=FALSE, xlab="Temp", ylab="Depth, cm")
axis(3)
axis(2, at=pretty(x), labels=100*(max(x)-pretty(x)))

for(j in 1:5)
{
Tx =  get.heat2(x, T1-T0, k, j*dt)
Tx = Tx+T0
lines(Tx, max(x)-x, lty=2, col=j)
}




}
\keyword{misc}
