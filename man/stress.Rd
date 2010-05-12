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
#####   no stress tensor provided => no mohr's circle
stress()

#########################  run program with a stress tensor
Stensor = matrix(c(
15, 0, 0,
0, 10, 0,
0,  0, 5), ncol=3)


stress(Stensor=Stensor)

P1 = c(0.2, 1, 1, 0)
P2 = c(1, 0.1, 1, 0)
P3 = c(1, 1, 0.4, 0)

  S = stressSETUP(P1, P2, P3, xscale=30   )

stress(PPs = S$PPs, Rview =S$Rview, xscale = S$xscale, Stensor=Stensor )




}


}
\keyword{misc}
