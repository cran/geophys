\name{points2line}
\alias{points2line}
\title{Points to Line
}
\description{Given a point and a matrix of lines,
  calculate the projection of the
  vector of the point to the first coordinate of
  each line and the perpendicular distance.
  
}
\usage{
points2line(Lp, VL)
}
\arguments{
  \item{Lp}{Point list (x,y)
}
\item{VL}{Matrix of points, N by 4
  as X1, Y1, X2, Y2
}
}
\details{
  The first point in the
  matrix is considered the base.
}
\value{
  list:
  \item{rat}{cosine projection }
  \item{srat}{sine projection }
   
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{REplane
}
\examples{


S= stressSETUP()
      pstart()
PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')
      PLOTplane(S$Rp, planecol="brown")

 basepoint = 3

        legpoints = c(7,4,2)

   

        VL =   cbind( rep(S$Rbox[basepoint,1] , length(legpoints)),
     rep(S$Rbox[basepoint,2] , length(legpoints)),
          S$Rbox[legpoints,1], S$Rbox[legpoints,2])

      Lp=list() 
Lp$x=c(40.0180732557)
Lp$y=c(40.4847345741)

 G = points2line(Lp, VL )

}
\keyword{misc}
