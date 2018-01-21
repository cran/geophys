\name{NinePointCircle}
\Rdversion{1.1}
\alias{NinePointCircle}
\title{Nine Point Circle
}
\description{Nine Point cirlce of a triangle
}
\usage{
NinePointCircle(P1, P2 = c(0, 1), P3 = c(1, 0), add = FALSE, SHOW = TRUE)
}
\arguments{
  \item{P1}{vector, Point 1
}
  \item{P2}{vector, Point 1
}
  \item{P3}{vector, Point 1
}
  \item{add}{add to existing plot
}
  \item{SHOW}{create a new plot and add
}
}
\details{circle passes through nine points that
  can be calculated for any triangle.
  Also known as Feuerbach's circle, Euler's circle,
  Terquem's circle, the six-point circle, the twelve-point circle, the
  n-point circle,
  the medioscribed circle,
  the mid circle,
  the  circum-midcircle.

  

}
\value{
  list of essential points:
  \item{A}{2-vector, vertex point 1}
  \item{B}{2-vector, vertex point 2}
  \item{C}{2-vector, vertex point 3}
  \item{D}{2-vector, mid-point opposite A}
  \item{E}{2-vector, mid-point opposite B}
  \item{F}{2-vector, mid-point opposite C}
    \item{G}{2-vector, foot altitude point opposite A }
  \item{H}{2-vector, foot altitude point opposite B}
  \item{I}{2-vector, foot altitude point opposite C}
  \item{J}{2-vector, mid point from S-A}
  \item{K}{2-vector, mid point from S-B}
  \item{L}{2-vector, mid point from S-C}
  \item{S}{2-vector, Intersection point of altitudes}
  \item{CEN}{2-vector, center of nine point circle}
  \item{R}{radius of nine point circle}
}
\references{
<http://en.wikipedia.org/wiki/Nine-point_circle>
  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{TriangleInfo, TriangleCenter
}
\examples{

  P1 = 10*runif(2)
      P2 =  10*runif(2)
      P3 =  10*runif(2)

      TRI =  NinePointCircle(P1, P2, P3, add=TRUE, SHOW=TRUE)



}
\keyword{misc}
