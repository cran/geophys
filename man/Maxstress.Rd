\name{Maxstress}
\alias{Maxstress}
\title{Maximum Stress 
}
\description{Maximum Normal and Shear stress in a plane
}
\usage{
Maxstress(NN, Stensor)
}
\arguments{
  \item{NN}{Normal Vector
}
  \item{Stensor}{Stress tensor
}
}
\details{Stress tensor can be
  entered as a vector of 3 or six values
  or as a 3 by 3 matrix.
  If three values are provided
  provided they are assumed to be the
  principal stresses.
  If six values are provided,
  they are c(x11, x12, x13, x22, x33).
  In either case eigen values and eigen vectors are
  calculated and used in the following
  analysis.
  
  
  
}
\value{
  \item{NN}{Normal vector to plane is returned}
  \item{sigNORMmax}{maximum normal stress}
  \item{tauSHEARmax}{maximum shear stress in the plane}
  
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{stress
}
\examples{
  S = stressSETUP( )

Nvec = NORMvec(S$PPs, S$xscale, S$Rview, S$aglyph , add = FALSE)

Stensor = matrix(c(
15, 0, 0,
0, 8, 0,
0,  0, 5), ncol=3)

 

Mstress  = Maxstress(Nvec, Stensor)

DoMohr(Stensor)
 axis(1)
axis(2)


points(Mstress$sigNORMmax , Mstress$tauSHEARmax, pch=21, col='blue'  , bg='gold' )

u=par('usr')

segments(0, Mstress$tauSHEARmax, Mstress$sigNORMmax ,
Mstress$tauSHEARmax, lty=2, col='green'  )

text(mean(c(0, Mstress$tauSHEARmax)),  Mstress$tauSHEARmax,
"MaxShear in Plane", pos=3)


segments(Mstress$sigNORMmax , u[3] , Mstress$sigNORMmax ,
Mstress$tauSHEARmax, lty=2, col='purple'  )

text(Mstress$sigNORMmax , u[3], "MaxNormal stress", adj=c(0,-1) )


}
\keyword{misc}
