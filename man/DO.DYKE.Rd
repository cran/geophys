\name{DO.DYKE}
\alias{DO.DYKE}
\title{Heat for Dyke Intrusion
}
\description{Solve heat equation for
  dyke intrusion.
}
\usage{
DO.DYKE(a = a, x = x, t = t, k = k, T0 = T0, NDIM = TRUE)
}
\arguments{
  \item{a}{Thickness of dyke (meters)
}
  \item{x}{vector of distances across dyke (m)
}
  \item{t}{time in seconds
}
  \item{k}{diffusion constant
}
  \item{T0}{Temperature at infinity
}
  \item{NDIM}{logical, TRUE=use nondimensional temperature for text labels
}
}

\value{graphical side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{get.heat2, DO.HALFSPACE
}
\examples{

 a = 10;
k=10^6


  DO.DYKE(NDIM=TRUE)

####  dimensional variable plot
   DO.DYKE(NDIM=FALSE)


}
\keyword{misc}
