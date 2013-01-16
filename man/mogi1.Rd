\name{mogi1}
\alias{mogi1}
\title{Mogi's model
}
\description{Mogi's model (point source in elastic half-space)
}
\usage{
mogi1(d = 1, f = 1, a = 0.1, P = 1e+05, mu = 4e+09, nu = 0.25)
}
\arguments{
  \item{d}{distance along surface, m

}
  \item{f}{depth to source, m

}
  \item{a}{radius of sphere injected, m

}
  \item{P}{hydrostatic pressure of injection, Pa

}
  \item{mu}{shear modulus

}
  \item{nu}{Poisson's ratio
}
}
\details{
Units should be constistent, e.g.: R, F, A, Ur and Uz in m imply
V in m3; E, mu and P in Pa; Dt in rad, Er, Et and nu dimensionless.

  
}
\value{
  list
  \item{ur}{radial displacement}
  \item{uz}{vertical displacement}
  
}
\references{
     
 Mogi, K., Relations between the eruptions of various volcanoes and the
 deformations of the ground surfaces around them, Bull. Earthquake Res.
 Inst. Univ. Tokyo, 36, 99-134, 1958.


  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{

P=1e5
a = 0.1
    f = 1

 nu=0.25
    E = 10e9

mu = E/(2*(1+nu));
  
#############
rho = seq(from=-3, to=3, by=0.1)
o1 = mogi1(a, P, rho, mu, f)
plot(rho, o1$uz)

plot(rho, o1$ud)






EX = seq(from=-3, by=.1, to= 3)

mm = RPMG::meshgrid(EX, EX)

rho = sqrt(mm$x^2 +  mm$y^2)

o1 = mogi1(a, P, rho, mu, f)

#############
image(EX, EX, o1$uz, asp=1)


contour(    EX, EX, o1$uz, add=TRUE)
#############
image(EX, EX, o1$ur, asp=1)


contour(    EX, EX, o1$ur, add=TRUE)

#############
image(EX, EX, o1$uz, asp=1)
contour(    EX, EX, o1$ur, add=TRUE)


}
\keyword{misc}
