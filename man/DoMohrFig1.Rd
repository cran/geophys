\name{DoMohrFig1}
\alias{DoMohrFig1}
\title{Annotated Stress Tensor
}
\description{Annotated plot of a 2D stress tensor
}
\usage{
DoMohrFig1(Stensor=matrix(c(5,1, 1, 3), ncol=2), rot1=NULL)
}
\arguments{
  \item{Stensor}{ stress tensor
}
  \item{rot1}{ rotation matrix
}
}

\value{Graphical Side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{DoMohr
}
\examples{

Stensor =matrix(
       c(50, 40,
         40, 10), ncol=2)

DoMohrFig1(Stensor)

}
\keyword{misc}
