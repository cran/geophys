\name{TriangleInfo}
\Rdversion{1.1}
\alias{TriangleInfo}
\title{Triangle Information
}
\description{Given a triangle (three non-colinear points, return info that is useful
}
\usage{
TriangleInfo(P1, P2 = c(0, 1), P3 = c(1, 0), add = FALSE)
}
\arguments{
  \item{P1}{vector or list(x,y)  of 2D point 1
}
  \item{P2}{vector or list(x,y)  of 2D point 2
}
  \item{P3}{vector or list(x,y)  of 2D point 3
}
  \item{add}{logical, if TRUE add a variety of lines to the plot
}
}
\details{
  If either P1 P2 or P3 are an x,y list, they are converted to vectors.
  
  If P2 and P3 are missing, P1 contains all the points.
  Input can be list of 3 x,y values, a 3 by 2 matrix .

  
}
\value{
  list
  \item{BI}{x,y location of center of inscribed circle}
  \item{CIRCUM}{x,y location of center of circumscribed circle}
  \item{IH}{x,y location of Intersection of bisectors}
  \item{CEN}{centroid location}
  \item{r}{radius of inscribed circle}
  \item{R}{radius of  circumscribed circle}
  \item{AngBis}{loci of vectors bisecting each corner on opposite leg}
  \item{H}{lengths of altitudes}
  \item{M}{lengthsof medians}
  \item{TEE}{lengths of perpendiculars}
  \item{Area}{Area of triangle}


}
\references{
http://en.wikipedia.org/wiki/Triangle_center

http://en.wikipedia.org/wiki/Encyclopedia_of_Triangle_Centers

  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{Intersect2points, Sect2vex, AXB.prod, TriangleCenter, rot2Zplane
}
\examples{



P1 = runif(2)
P2 = runif(2)
P3 = runif(2)

  

theX = c(P1[1],P2[1], P3[1])
theY = c(P1[2],P2[2], P3[2])


        plot(theX , theY, asp=1)
        
        
        points(theX , theY)
        text(theX , theY, labels=1:3, pos=3)
        lines(c(theX, theX[1]), c(theY, theY[1]) )

TriangleInfo(P1, P2, P3, add=TRUE)
#####    TriangleInfo(rbind(P1, P2, P3),   add=TRUE)




}
\keyword{misc}
