\name{okada85}
\alias{okada85}
\title{Okada Fault Deformation Model
}
\description{Calculate the
  deformation due to a fault
  based on Okada 1985 code.
}
\usage{
okada85(e = 0, n = 0, depth = 0, strike = 20,
dip = 20, L = 5, W = 3, rake = 20, slip = 1, U3 = 1, nu = 0.25)
}
\arguments{
  \item{e}{  vector or matrix, Easting Observation, km
}
  \item{n}{  vector or matrix, Northing Observation, km
}
  \item{depth}{depth of fault
}
  \item{strike}{degrees, strike of fault, measured from north
}
  \item{dip}{degrees, strike of fault, measured from horizontal
}
  \item{L}{ Length, km
}
  \item{W}{ Width, km
}
  \item{rake}{ rake (degrees)
}
  \item{slip}{ slip
}
  \item{U3}{ dilitation (m)
}
  \item{nu}{ Poisson's ration, default=0.25
}
}
\details{
       Computes displacements, tilts and strains at the surface of an
       elastic half-space, due to a dislocation defined by RAKE, SLIP,
       and OPEN on a rectangular fault defined by orientation STRIKE and
       DIP, and size LENGTH and WIDTH. The fault centroid is located
       (0,0,-DEPTH).

       R-Code is a translation of original Matlab code by Francois
       Beauducel.  UNC student Lan Ngo worked on the first version and
       J. M. Lees tested, verified and finalized the package.


  
}
\value{
list:
\item{uE}{East deformation}

  \item{uN}{North deformation }
  \item{uZ}{Vertical deformation}
  \item{uZE}{tilt vertical east}
  \item{uZN}{tilt vertical north}
  \item{uNN}{horizontal strains North-North}
  \item{uNE}{horizontal strains North-East}
  \item{uEN}{horizontal strains East-North}
  \item{uEE}{horizontal strains East-East}
}
\references{
         Aki K., and P. G. Richards (1980) (ISBN-13: 978-1891389634) Quantitative seismology, Freemann & Co,
         New York.
	 
          Okada Y., Surface deformation due to shear and tensile faults in a
             half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.

	     Acknowledgments: Francois Beauducel, Dmitry Nicolsky,
	     University of Alaska

	     

  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{
  Units should be in km.
  
  Note that vertical strain components can be obtained with following equations:
          uNZ = -uZN;
          uEZ = -uZE;
          uZZ = -(uEE + uNN)*NU/(1-NU);

Be careful about the definition of strike: seems that here the
strike may be measured positive north of east, not east of north as stated
in the documentation above.


  
}


\seealso{mogi1, mogiM    
}
\examples{

library(RPMG)

EX = seq(from=-10, to=10, length=50)
WHY = seq(from=-10, to=10, length=50)
M = RPMG::meshgrid(EX, WHY )

OKout  = okada85(M$x, M$y,2,30,70,5,3,-45,1,1) 
old.par <- par(no.readonly = TRUE) 

par(mfrow=c(2, 2))

image(EX, WHY, OKout$uN, col=terrain.colors(100) ,
asp=1, main="uN" );contour(EX, WHY, OKout$uN, add=TRUE)

image(EX, WHY, OKout$uE, col=terrain.colors(100) ,
asp=1 , main="uE" );contour(EX, WHY, OKout$uE, add=TRUE)

image(EX, WHY, OKout$uZ, col=terrain.colors(100) ,
asp=1, main="uZ" );contour(EX, WHY, OKout$uZ,  add=TRUE)

plot(EX, OKout$uZ[, 25 ])


par(old.par)


}
\keyword{misc}
