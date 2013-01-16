\name{rev2RH}
\alias{rev2RH}
\title{Reverse to Right Hand
}
\description{Reverse orientation of the vectors
  to a right handed polygon.
}
\usage{
rev2RH(pol)
}
\arguments{
  \item{pol}{List with x and y components
}
}
\details{
  program reverses the order of the x and y components
  so they traverse in the oposite direction.
}
\value{list
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{dircheck, flipZEE
}
\examples{
L=list()
L$x=c( 0.2188,-0.6668,-1.0992,-0.5105, 0.5106)
L$y=c( 0.76822, 0.75292,-0.01741,-0.74692,-0.48165)
rev2RH(L)

}
\keyword{misc}
