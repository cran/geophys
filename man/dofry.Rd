\name{dofry}
\Rdversion{1.1}
\alias{dofry}
\title{Fry Estimation
}
\description{Calculate the Fry diagram for estimating shear on
  random point patterns
}
\usage{
dofry(x, y, PLOT = FALSE)
}
\arguments{
  \item{x}{x-coordinates
}
  \item{y}{y-coordinates
}
  \item{PLOT}{logical, TRUE=plot
}
}
\details{The number of points returned
  can be large and plotting can take a long time, so
  when plotting can restrict the points by using plotfry.
}
\value{
  List:
  \item{X}{x-coordinates}
  \item{Y}{y}
  \item{mx}{mean X}
  \item{my}{mean Y}
  
}
\references{
Fry, N., (1979) <http://dx.doi.org/10.1016/0040-1951(79)90135-5>
Random point distributions and strain measurement in rocks
Tectonophysics, 60:89-105.
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{Showfry, plotfry, xtractlip
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
