\name{vecproj}
\alias{vecproj}
\title{Vector Projection 
}
\description{Vector Projection information, such as angle
  and distances between points
}
\usage{
vecproj(P1, P2)
}
\arguments{
  \item{P1}{Point 1
}
  \item{P2}{Point 2
}
}
\details{
  The distances returned are the legs of right right triangles
  where the cosine of the angle is used to get the projection
  distance of the opposite side on the specified direction.
  See the example for an illustration.
}
\value{
  cang=cang, angrad=angrad, angdeg=angdeg, dis1=d1, dis2=d2
  \item{cang}{cosine of angle between points}
  \item{angrad}{angle, radians}
  \item{angdeg}{angle, degrees}
  \item{dis1}{distance}
  \item{dis2}{distance}
  
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{perpproj
}
\examples{

P1 = c(2, 3)
P2 = c(5, 2)


I = vecproj(P1, P2)


plot(c(0, P1[1], P2[1]), c(0, P1[2], P2[2]), asp=1, ann=FALSE)

arrows(0, 0, P1[1],  P1[2], length=.1)
arrows(0, 0,  P2[1],  P2[2], length=.1)

text(P1[1],  P1[2], "Point 1", pos=3)
text(P2[1],  P2[2], "Point 2", pos=3)

j1 = atan2(P1[2], P1[1])*180/pi
j2 = atan2(P2[2], P2[1])*180/pi

L1 = vlength(P1)
L2 = vlength(P2)

A = GEOmap::darc(L1*.2, j1, j2, 0, 0, n=1)

lines(A)
an = length(A$x)
arrows(A$x[an-1] , A$y[an-1] ,A$x[an] , A$y[an] ,  length=.08 )

text(A$x[an/2] , A$y[an/2], labels=format(I$angdeg, digits=4)  , pos=4)


 V1 = c( 0,P1[1], 0,  P1[2])
   V2 = c( 0,P2[1], 0,  P2[2])
  

    PP = perpproj( V1, V2, add=FALSE  )

arrows(P1[1],P1[2],PP$P2[1], PP$P2[2],  length=.07, lty=2, col='red') 
arrows(P2[1],P2[2],PP$P1[1], PP$P1[2],  length=.07, lty=2, col='blue') 

 labelLine( c(0, 0) , PP$P2 , lab="dis1", dinch = .25, aty=1,
acol='blue', above=FALSE  )


 labelLine( c(0, 0) , PP$P1 , lab="dis2", dinch = .25, aty=1,
acol='blue'  )


}
\keyword{misc}
