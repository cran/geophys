\name{mohrleg}
\alias{mohrleg}
\title{Legend for Mohr
}
\description{Legend for Mohr
}
\usage{
mohrleg(ES)
}
\arguments{
  \item{ES}{Eigen Value Decomposition, output of function eigen
}
}
\details{Add notes to plots of Mohr's circles.  Uses the eigenvalues of the decomposition.
}
\value{Graphical Side Effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{DoMohr
}
\examples{
Stensor = matrix(c(50, 40, 40, 10), ncol=2)

DoMohr(Stensor) 


}
\keyword{misc}
