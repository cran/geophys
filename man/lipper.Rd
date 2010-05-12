\name{lipper}
\alias{lipper}
\title{Ellipse calculator
}
\description{Return parametric version of ellipse from equation
  of an ellipse.
}
\usage{
lipper(a, b, c, d = 0, f = 0, g = -1)
}
\arguments{
  \item{a}{numeric
}
  \item{b}{numeric
}
  \item{c}{numeric
}
  \item{d}{numeric
}
  \item{f}{numeric
}
  \item{g}{numeric
}
}
\details{
  General quadratic equation is,
  \eqn{a x^2 + 2 b x y + c y^2 + 2 d x + 2 f y + g = 0}
  
  
}
\value{
  vector=c(ap, bp, phi)
}
\references{
http://mathworld.wolfram.com/Ellipse.html
  
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\examples{

LIP = lipper(10, 3, 9)

theta  = seq(from=0, to=2*pi, length=360)

        phi = LIP[3]
        px = LIP[1]*cos(theta)*cos(phi)-  LIP[2]*sin(theta)*sin(phi)
        py = LIP[1]*cos(theta)*sin(phi) + LIP[2]*sin(theta)*cos(phi)


plot(px, py, type='p', asp=1)
lines(px, py)


}
\keyword{misc}
