\name{Showfry}
\Rdversion{1.1}
\alias{Showfry}
\title{Show Fry Diagrams
}
\description{Show Fry Diagrams with random data
}
\usage{
Showfry(RDAT, shear = matrix(c(1, 1.2, 0, 1)), rad = 75)
}
\arguments{
  \item{RDAT}{data list
}
  \item{shear}{2D shearing matrix
}
\item{rad}{radius of points from the center
  from which to select the points for analysis.
}
}
\details{
  The original data is
  deformed and a circular subset is extracted
  for analysis.
  This is to show that the
  shape of the data does not affect the
  assessment of the fry ellipse.
  Plots input data and fry diagram, side by side.
  
}
\value{graphical side effects
}
\references{
  Fry, N., (1979)
Random point distributions and strain measurement in rocks
Tectonophysics, 60:89-105.

}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{dofry, plotfry, xtractlip
}
\examples{

RDAT = randFRY(400, LIM=c(0,0, 200, 200) , rlen=5   )
length(RDAT$x)
plot(RDAT$x, RDAT$y, asp=1, pch=".", cex=2)

u = par( no.readonly = TRUE)


par(mfrow=c(3,2))
    shr = 0.0
simpleshear = matrix(c(1, shr, 0,  1), ncol=2)

Showfry(RDAT, simpleshear, 75)
    shr = 1.2
simpleshear = matrix(c(1, shr, 0,  1), ncol=2)

Showfry(RDAT, simpleshear, 75)

epsilon1 = 0.4
H = matrix(c(1+epsilon1, 0, 0,  1/(1+epsilon1) ), ncol=2)

Showfry(RDAT, H, 75)

par(u)


}
\keyword{misc}
