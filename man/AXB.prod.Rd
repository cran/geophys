\name{AXB.prod}
\Rdversion{1.1}
\alias{AXB.prod}
\title{Vector Cross Product
}
\description{Vector Cross Product
}
\usage{
AXB.prod(A, B)
}
\arguments{
  \item{A}{vector, 2 or 3 elements
}
  \item{B}{vector, 2 or 3 elements
}
}
\details{Can be 2D or 3D vectors.  If lists are input they are converted
  to vectors.
}
\value{vector
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\examples{
  AXB.prod(c(1,0,0), c(0,1,0)) ###  =  c(0, 0 , 1)
  AXB.prod(c(0,1,0), c(1,0,0)) ###  =  c(0, 0 ,-1)



}
\keyword{misc}
