\name{flipZEE}
\alias{flipZEE}
\title{Flip Z direction
}
\description{Flip the zee direction
  so negative is down.
}
\usage{
flipZEE(pol)
}
\arguments{
  \item{pol}{list of x,y
}
}
\details{
  Reverse the sign of the y-direction
  to achieve negative down Z values.
}
\value{same as input with -y
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{Uses the y component as if it were Z.
}


\seealso{rev2RH, dircheck
}
\examples{
L=list()
L$x=c( 0.2188,-0.6668,-1.0992,-0.5105, 0.5106)
L$y=c( 0.76822, 0.75292,-0.01741,-0.74692,-0.48165)
flipZEE(L)



}
\keyword{misc}
