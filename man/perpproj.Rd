\name{perpproj}
\alias{perpproj}
\title{Perpendicular Projection
}
\description{Perpendicular Projection in two-dimensions
}
\usage{
perpproj(V1, V2, add = FALSE)
}
\arguments{
  \item{V1}{vector 1 
}
  \item{V2}{vector 1 
}
  \item{add}{logical, TRUE add to plot
}
}
\details{
  V1 and V2 are arrays with, for example, V1=c(x1 , x2 , y1, y2)
  or lists with V2=list(x=c(x1, x2), y-c(y1, y2)) which gets converted
  into the first format.

  The points are from the head of each vector projected
  onto the the line formed by the other vector.
  If the the x1 y1 of each vector do not coincide,
  the interesections of two lines is determined
  and the tail of each vector is moved to that point
  for determining the projection positions.

  
}
\value{
  list:
  \item{P1}{x, y}
  \item{P2}{x, y}
   
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{vecproj
}
\examples{


  V1 = c( 2,  6, 1, 9)
    V2 = c( 0, 5,  1, 2)

  
    PP = perpproj( V1, V2, add=FALSE  )


    R = range(c(V1, V2, unlist(PP) ))

    plot(R, R, type='n', asp=1)

    arrows(V1[1], V1[3], V1[2], V1[4], length=.08 ) 
   arrows(V2[1], V2[3], V2[2], V2[4], length=.08 ) 

    points(PP$P1[1],PP$P1[2], col='red')
       points(PP$P2[1],PP$P2[2], col='blue')



    arrows(V2[2], V2[4],PP$P1[1],PP$P1[2] , length=.08, col='red' )
    
    arrows(V1[2], V1[4], PP$P2[1],PP$P2[2], length=.08, col='blue' ) 
 


}
\keyword{misc}
