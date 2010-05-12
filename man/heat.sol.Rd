\name{heat.sol}
\alias{heat.sol}
\title{Heat Equation Solution
}
\description{Heat Equation Solution
}
\usage{
heat.sol(x, T0, k, t)
}
\arguments{
  \item{x}{distance
}
  \item{T0}{temperature
}
  \item{k}{diffusion constant
}
  \item{t}{time vector
}
}
\details{INteractive code to see evolution
}
\value{Graphical Side Effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{get.heat, get.heat2,DO.HALFSPACE 
}
\examples{


\dontrun{
k = 1e-6
 t  = 1*365*24*60*60
T0 = 1000
x = seq(0,20, length=1000)
T =  heat.sol(x, T0, k, t)

}

}
\keyword{misc}
