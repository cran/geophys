\name{centroid}
\alias{centroid}
\title{centroid
}
\description{calculate the centroid of a polygon
}
\usage{
centroid(p)
}
\arguments{
  \item{p}{list: x,y
}
}
\details{same code as centroid.polygon
}
\value{vector x, y
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{

L=list()
L$x=c( 0.2188,-0.6668,-1.0992,-0.5105, 0.5106)
L$y=c( 0.76822, 0.75292,-0.01741,-0.74692,-0.48165)
centroid(L)



}
\keyword{misc}
