\name{stress}
\alias{stress}
\title{Stress Box
}
\description{Display interactive Stress Box
  illustrating normal and shear stress on
  an arbitrary plane given a stress tensor.
}
\usage{
stress(PPs = matrix(ncol = 4, nrow = 3), Rview = c(-130, -50), xscale = 100, Stensor = matrix(ncol = 3, nrow = 3))
}
\arguments{
  \item{PPs}{3Points in 3 by 4 matrix
}
  \item{Rview}{Viewing matrix
}
  \item{xscale}{arbitrary scale
}
  \item{Stensor}{stress tensor
}
}
\details{Interactive program to examine stress.
  To change the plane orientation,
  click on one of the points and then on
  a new position along the axis selected.

  If a stress tensor is provided, then a second graphics device
  is  required to show the Mohr's circle.

  
  
}
\value{Graphical Side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{
  Uses RPMG for interaction.
}


\seealso{pstart, PLOTbox, PLOTplane, NORMvec
}
\examples{

\dontrun{

stress()

Stensor = matrix(c(
15, 0, 0,
0, 10, 0,
0,  0, 5), ncol=3)


stress(Stensor=Stensor)

}


}
\keyword{misc}
