\name{Sect2vex}
\Rdversion{1.1}
\alias{Sect2vex}
\title{Intersect 2 Vectors
}
\description{Finder intersection point of two vectors in a plane
}
\usage{
Sect2vex(V1, V2)
}
\arguments{
  \item{V1}{vector of 2 points
}
  \item{V2}{vector of 2 points
}
}
\details{if vectors are parallel return NULL
}
\value{x-y location of intersection
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{AXB.prod
}
\examples{

v1 = matrix(runif(4), ncol=2)
v2 =  matrix(runif(4), ncol=2)


pt = Sect2vex(v1, v2)

plot(c(v1[,1],v2[,1],pt[1] ) , c(v1[,2],v2[,2],pt[2]) , asp=1)

arrows(v1[1,1], v1[1,2], v1[2,1], v1[2,2])
arrows(v2[1,1], v2[1,2], v2[2,1], v2[2,2])

points(pt[1], pt[2], col='red', cex=2, pch=3)






}
\keyword{misc}
