\name{PolarDecomp}
\alias{PolarDecomp}
\title{Polar Decomposition
}
\description{Polar Decomposition for Strain
}
\usage{
PolarDecomp(A)
}
\arguments{
  \item{A}{Strain Matrix
}
}
\details{
  Polar decomposition uses the svd
  to extract 2 matrices that represent the
  stretch and rotation of a strain: A = UP.
  U is orthogonal rotation matrix and P is the stretch tensor.
  These are extracted from the singular value decomposition.
  
}
\value{
  \item{P}{matrix, stretch tensor}
  \item{U}{matrix, orthogonal rotation matrix}
}
\references{
  <http://en.wikipedia.org/wiki/Finite_strain_theory>
  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{svd
}
\examples{


A = matrix(runif(4, -1, 1), ncol=2)
PD = PolarDecomp(A)

E = svd(A)

###  W         S           V
 E$u %*% diag(E$d)%*%  t(E$v)

 ###t(E$v) %*% diag(E$d)%*%  (E$u)

P = E$v%*% diag(E$d)%*%  t(E$v)
U =  E$u %*% t(E$v)


U %*% P 


}
\keyword{misc}
