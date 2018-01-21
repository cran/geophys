\name{plotfry}
\Rdversion{1.1}
\alias{plotfry}
\title{Plot Fry Output
}
\description{Plot Fry Output limited by a specified distance.
}
\usage{
plotfry(fry, dis, col = grey(0.7), ann = FALSE, axes = FALSE)
}
\arguments{
  \item{fry}{list output of dofry
}
  \item{dis}{distance to limit plotting to the center
}
  \item{col}{color for points
}
  \item{ann}{logical, annotation of axes, see par
}
  \item{axes}{logical, whether to plot the axes, see par
}
}
\details{Used for plotting the points
  attained through the fry method
}
\value{
  \item{x}{x-coordinates}
  \item{y}{y-coordinates}
  \item{mx}{x center of plot}
  \item{my}{y center of plot}
  \item{dis}{distance for limiting plot}
}
\references{
Fry, N., (1979) <http://dx.doi.org/10.1016/0040-1951(79)90135-5>
Random point distributions and strain measurement in rocks
Tectonophysics, 60:89-105.

}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{dofry, Showfry
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
