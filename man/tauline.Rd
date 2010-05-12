\name{tauline}
\alias{tauline}
\alias{tauplane}

\title{Shear Stress along Line
}
\description{Calculate the shear stress along an arbitrary line in a
  plane with stress orientation
}
\usage{
tauline(Rp, P1, P2, Rview, ES, NN)
tauplane(Rp, L, Rview, ES, NN)


}
\arguments{
  \item{Rp}{ rotated points describing plane
}
  \item{P1}{point 1 extracted from screen (locator)
}
  \item{P2}{point 2 extracted from screen
}
  \item{Rview}{rotation matrix for viewing
}
  \item{ES}{eigen value decomposition from eigen
}
  \item{NN}{normal vector to plan in unrotated coordinates
  }
\item{L}{list locations (x,y) in the figure, projected to the plane }
  
}
\details{
  Used internally in stress.
  When the plan is plotted, if two points are
  located on the figure,
  the points are positions on the plan and un-rotated
  using the Rview matrix.
  Then the shear stress in the plan along that line
  is calculated and returned.
}
\value{shear stress along the line indicated
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{stress,NORMvec 
}
\examples{


S= stressSETUP()


pstart()

     PLOTplane(S$Rp, planecol="brown")
     PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')




##  L = locator(2)

L=list()
L$x=c(-13.6305297057, 52.6412739525)
L$y=c(26.2697350325,32.4501696158)


Stensor = matrix(c(
15, 0, 0,
0, 10, 0,
0,  0, 5), ncol=3)


P1 = list(x=L$x[1], y=L$y[1])
P2 = list(x=L$x[2], y=L$y[2])


ES = eigen(Stensor) 
NN = NORMvec(S$PPs, S$xscale, S$Rview, aglyph=S$aglyph, add=FALSE)

tauline(S$Rp, P1, P2, S$Rview, ES, NN)





}
\keyword{misc}
