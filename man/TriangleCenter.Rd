\name{TriangleCenter}
\Rdversion{1.1}
\alias{TriangleCenter}
\title{Triangle Center
}
\description{Extract Triangle center in 3D
}
\usage{
TriangleCenter(P1, P2, P3, A1= 0, A2= 360, KNum=10)
}
\arguments{
  \item{P1}{3-vector, point(x,y,z)
}
  \item{P2}{3-vector, point(x,y,z)
}
  \item{P3}{3-vector, point(x,y,z)
  }
  \item{A1}{degrees, initial angle in plane(default=0)
  }
  \item{A2}{degrees, final angle in plane(default=360)
  }
  \item{KNum}{Divisor Number to divide range by (default=10)
    
  }



  
}
\details{
  Program rotates the object to the X-Y plane
  and does calculations in 2D, then rotates back.
}
\value{
  \item{Center}{x-y of center of the inscribed circle}
  \item{r}{radius of inscribed}
  \item{Cinscribed}{inscribed circle points around center}
    \item{CIRCUM}{x-y of center of the circumscribed circle}


  
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{TriangleInfo
}
\examples{

S = stressSETUP()
pstart()

PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')

pstart()

PLOTplane(S$Rp, planecol="brown")
PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')


NORMvec(S$PPs, S$xscale, S$Rview, S$aglyph, add=TRUE)
   P1 = S$PPs[1, 1:3] 
    P2 = S$PPs[2, 1:3] 
    P3 = S$PPs[3, 1:3]

BV = TriangleCenter(S$PPs[1,1:3],S$PPs[2,1:3], S$PPs[3,1:3] )

 CIRCview =   BV$Cinscribed  %*% S$Rview

    lines(CIRCview[,1], CIRCview[,2], col='purple')

cview =   BV$Center  %*% S$Rview

    points(cview[1,1], cview[1,2])


}
\keyword{misc}
