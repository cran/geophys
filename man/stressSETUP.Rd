\name{stressSETUP}
\alias{stressSETUP}
\title{Setup Stress Box
}
\description{Setup default stress box, plane information and normal
  vector 3d plotting glyph
  
}
\usage{
stressSETUP(P1=c(.2, 1,1,0 ), P2=c(1, .1,1,0), P3=c(1, 1,.4,0), xscale=30)
}
\arguments{
  Arguments are optional
  \item{P1}{4-vector, point 1}
  \item{P2}{4-vector,point 2}
  \item{P3}{4-vector,point 3}
  \item{xscale}{scale factor for plotting}

  }
\details{Initialize the stress set up.
  
}
\value{
List:
\item{xscale}{arbitrary scale}
\item{Rview}{viewing rotation matrix}
\item{BOX}{Box coordinates, scaled}
\item{AX}{XYZ Axis coordinates}
\item{Rbox}{rotated box}
\item{Rax}{XYZ rotated axes}
\item{PPs}{Points along the Box axes}
\item{Rp}{rotated points}
\item{aglyph}{3D glyph for plotted norm}
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{stress, PLOTbox, PLOTplane, NORMvec
}
\examples{


S = stressSETUP()
pstart()

PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')

pstart()

PLOTplane(S$Rp, planecol="brown")
PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')


NORMvec(S$PPs, S$xscale, S$Rview, S$aglyph, add=TRUE)


}
\keyword{misc}
