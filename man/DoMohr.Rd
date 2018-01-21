\name{DoMohr}
\alias{DoMohr}
\title{Mohr's Circle
}
\description{Mohr's Circle
}
\usage{
DoMohr(Stensor=diag(c(3,2,1)), axis=NULL)
}
\arguments{
  \item{Stensor}{Stress tensor
  }
  \item{axis}{vector of  axes to plot, c(1,2,3,4) }

  
}
\details{Stress tensor is a 2 by 2 or 3 by 3 symmetric matrix
}
\value{Graphical Side effects
}
\references{

  W.D. Means,  Stress and Strain: Basic Concepts of Continuum Mechanics
  for Geologists,  Springer; 1976. 4th ed., 339 p. (ISBN-10:0387075569).
  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{DoMohrFig1, stress
}
\examples{
Stensor = matrix(c(2,.3, .3, 1), ncol=2)

Stensor = matrix(c(50, 40, 40, 10), ncol=2)

DoMohr(Stensor) 


Stensor = matrix(c(
-50, 40, 20,
40, 10, 10,
20,  10, 5), ncol=3)

DoMohr(Stensor) 


}
\keyword{misc}
