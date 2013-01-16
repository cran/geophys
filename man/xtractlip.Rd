\name{xtractlip}
\Rdversion{1.1}
\alias{xtractlip}
\title{Extract Convex Hull Ellipse
  from fry
}
\description{Extract Convex Hull Ellipse
  from fry-diagram
}
\usage{
xtractlip(AF)
}
\arguments{
  \item{AF}{x-y List out put of plotfry
    
}
}
\details{To get closest points to center,
  program cycles around in a circular manner pulling out the closest
  points.
  These points are used to estimate the chull ellipse.
  
}
\value{List: see output of ellipsoidhull:
     an object of class '"ellipsoid"', basically a 'list' with several
     components
  
}
\references{

   Fry, N., (1979)
Random point distributions and strain measurement in rocks
Tectonophysics, 60:89-105.


  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{ellipsoidhull
}
\examples{

RDAT = randFRY(400, LIM=c(0,0, 200, 200) , rlen=5   )
length(RDAT$x)
plot(RDAT$x, RDAT$y, asp=1, pch=".", cex=2)
flag = sqrt( (RDAT$x-mean(RDAT$x))^2  + (RDAT$y-mean(RDAT$y))^2)<75

DAT = list(x=RDAT$x[flag], y=RDAT$y[flag])
 x = DAT$x
    y = DAT$y
 plot(x,y, asp=1, pch=21, col='red', bg='gold', ann=FALSE, axes=FALSE)

FF = dofry(DAT$x, DAT$y )
AF = plotfry(FF, dis=30)

Z = xtractlip(AF)

lines(cluster::predict.ellipsoid(Z$hull), col='red')  
   

}
\keyword{misc}
