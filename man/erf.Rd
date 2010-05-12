\name{erf}
\alias{erf}
\alias{erfc}
\alias{erfinv}

\alias{philpotts.erf}

\title{Error Function
}
\description{Error Function
}
\usage{
erf(x)
philpotts.erf(x)
erfc(x)
erfinv(x)

}
\arguments{
  \item{x}{
}
}
\details{Uses pnorm to calculate the error function.
  erfinv returns the inverse of the erf function.
  
}
\value{error function
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{philpotts.erf
}
\examples{

j = seq(from=0, to=5.5, by=0.01)
 y = erf(j)   
phil = philpotts.erf(j)

 plot(j,y)
points( j, phil, pch=3, col='red'  )


#### plot the difference
plot(j, y-phil)
lines(j, y-phil)


#####  inverse of erf function:

 j = seq(from=0, to=5.5, by=0.1)
 y = erf(j)

for(i in 1:length(j))
   {

     z = erfinv(y[i])
     print(paste(i, y[i], j[i], z, sep=" "))
   }


}
\keyword{misc}
