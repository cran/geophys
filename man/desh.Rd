\name{desh}
\alias{desh}
\title{Display Mesh
}
\description{Display Mesh
}
\usage{
desh(M, add = TRUE, PTS = TRUE, colmesh = grey(0.8), colpts = grey(0.5),
...)
}
\arguments{
  \item{M}{Mesh output of meshgrid
}
  \item{add}{logical, add to plot
}
  \item{PTS}{logical, add points
}
  \item{colmesh}{color for mesh
}
  \item{colpts}{color for points
  }
  \item{\dots}{graphical parameters from par
  }
  

  
}
\details{Used for showing strain
  mesh plots.
}
\value{graphical side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{meshgrid
}
\examples{

###  warping of the ground from a Mogi source
 P=1e5
 a = 0.01
         f = .1
      nu=0.25
         E = 10e13
     mu = E/(2*(1+nu));
    
EX = seq(from=-3, by=.4, to= 3)
     mm = meshgrid(EX, EX)
  rho = sqrt(mm$x^2 +  mm$y^2)

cen = list(x=mean(EX), y=mean(EX))

are = sqrt( (mm$x-cen$x)^2  + (mm$y-cen$y)^2 )


 o1 = mogi1(a, P, rho, mu, f)


shiftx =  o1$ur*(mm$x-cen$x)/are
shifty =  o1$ur*(mm$y-cen$y)/are 

newmesh = list(x=mm$x+shiftx , y=mm$y+shifty)

desh(mm, add=FALSE, PTS=FALSE, colmesh=grey(.8) )


desh(newmesh, add=TRUE, PTS=FALSE, colmesh=rgb(.6,.8, .6) )






}
\keyword{misc}
