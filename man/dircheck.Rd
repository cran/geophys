\name{dircheck}
\alias{dircheck}
\title{Check direction of polygon
}
\description{Check direction of polygon for Right Handed orietnation
}
\usage{
dircheck(pol)
}
\arguments{
  \item{pol}{list of x,y
}
}
\details{
  Gets cross product of consecutive points
  in polygon and returns the sign of each.
  If most are positive, should try to reverse.
}
\value{
  vector of signs
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{
  Co-linear values should return zero.
}


\seealso{rev2RH, flipZEE
}
\examples{

L=list()
L$x=c( 0.2188,-0.6668,-1.0992,-0.5105, 0.5106)
L$y=c( 0.76822, 0.75292,-0.01741,-0.74692,-0.48165)
dircheck(L)

K = rev2RH(L)
###  after reversing, 
dircheck(K)




}
\keyword{misc}
