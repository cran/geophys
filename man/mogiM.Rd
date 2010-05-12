\name{mogiM}
\alias{mogiM}
\title{Mogi MOdel
}
\description{Mogi model deformation returns the
  deformation from a point source
  presurized inflation in an elastic medium.

Mogi's model (point source in elastic half-space).
       computes radial
       and vertical displacements Ur and Uz, ground tilt Dt, radial and
       tangential strain Er and Et on surface, at a radial distance R
       from the top of the source due to a hydrostatic pressure inside a
       sphere of radius A at depth F, in a homogeneous, semi-infinite elastic
      body and approximation for A << F (center of dilatation). Formula by
      Anderson [1936] and Mogi [1958].



  
  
}
\usage{
mogiM(R = 1, F = 1, A = 0.1, P = 1e+05, E = 1e+10, nu = 0.25)
}
\arguments{
  \item{R}{Hoirizontal Distance frm source, m
}
  \item{F}{Depth below surface, m, positive down
}
  \item{A}{radius of magma chamber
}
  \item{P}{hydrostatic pressure change in the sphere
}
  \item{E}{elasticity (Young's modulus)
}
  \item{nu}{Poisson's ratio
}
}
\details{
  Original paper by Mogi used poisson's ratio equale to
  0.25, i.e. lame parameters lambda and nu were equal.
}
\value{
  list:
  \item{ur}{radial  displacements Ur}
  \item{uz}{vertical displacements Uz, Uz > 0 = UP}
  \item{dt}{ground tilt Dt}
  \item{er}{radial strain Er}
  \item{et}{tangential strain  Et on surface}
}
\references{

        Anderson, E.M., Dynamics of the formation of cone-sheets, ring-dikes,
             and cauldron-subsidences, Proc. R. Soc. Edinburgh, 56, 128-157, 1936.

	     Mogi, K., Relations between the eruptions of various volcanoes and the
             deformations of the ground surfaces around them, Bull. Earthquake Res.
              Inst. Univ. Tokyo, 36, 99-134, 1958.

}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\examples{


data(PXY)

delV = 2.3E13/(100^3)    ##### (convert to meter^3 from cm^3)
F = 2.8E5/100      ##### (convert to meter from cm    )

EX = seq(from=0, by=100, to= 9000)


Atest = mogiM(R=EX,F=F,A=delV)



 plot(PXY, pch=6, col='purple', xlim=c(0,9), ylim=c(0, 1) )
    ###  model
    lines(EX/1000, Atest$uz/max(Atest$uz))


############ best fit   optimization

library(stats)

    fr<-function(x)
      {

        Atest = mogiM(R=PXY$x*1000 ,F=x[1],A=x[2])


        rms = sum ( (PXY$y - Atest$uz/max(Atest$uz))^2 )

        return(rms)
      }
xin = c(2600, 2.0e+07)

FOUT = optim(xin , fr)

  Btest = mogiM(R=EX,F=FOUT$par[1] ,A=FOUT$par[2])

   plot(PXY, pch=6, col='purple', xlim=c(0,9), ylim=c(0, 1) )
 
 lines(EX/1000, Btest$uz/max(Btest$uz))







    
}
\keyword{misc}
