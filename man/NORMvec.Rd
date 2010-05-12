\name{NORMvec}
\alias{NORMvec}
\title{Plot Normal Vector
}
\description{Calculate and plot a normal vecotor to a plane
}
\usage{
NORMvec(PPs, xscale, Rview, aglyph = list(), add = TRUE)
}
\arguments{
  \item{PPs}{Points for plotting, these define the plane
}
  \item{xscale}{scale for the figure
}
  \item{Rview}{viewing matrix
}
  \item{aglyph}{glyph for plotting the normal vector in 3D
}
  \item{add}{logical, whether to add to current plot
}
}

\value{Graphical Side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{Used internally
}


\seealso{pstart, PLOTbox, PLOTplane, stress
}
\examples{

#######    set colors
 axcol = 'black'
 boxcol = 'blue'
 planecol = 'brown'

####   view angle
 Rview  =    ROTZ(-130) %*% ROTX(-50) 

###   set arbitrary scale
    xscale = 100

###  set glphy for plotting 3D vectors
 headlen =xscale* .3/6
       len =xscale* .7/6
       basethick =xscale* 0.05/2
       headlip =xscale* .02/2
       aglyph = Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )

####   set three points and combine them with rbind
     P1 = xscale*c(.2, 1,1,0)
    P2 = xscale*c(1, .1,1,0)
    P3 = xscale*c(1, 1,.4,0)

    PPs = rbind(P1, P2, P3)
    
  ###  convert points for view
Rp = PPs  %*% Rview

###   create box  and scale it
   BOX <-matrix(c(0,0,0,0,
            0, 1, 0,0,
            0, 1, 1,0,
            0, 0, 1,0,
            1,0,0,0,
            1, 1, 0,0,
            1, 1, 1,0,
            1, 0, 1,0), ncol=4, byrow=TRUE)


    BOX = xscale*BOX

    ##############   create axes and scale them
AX = matrix(c(0,0,0,0,
            1, 0, 0,0,
            0, 0, 0,0,
            0, 1, 0,0,
            0,0,0,0,
            0, 0, 1,0), ncol=4, byrow=TRUE)

    AX = 1.5*xscale*AX

 ##############   rotate axes and box
    Rax =  AX %*% Rview

    
          Rbox =   BOX %*% Rview

 ##############  start the figure
pstart(xscale=xscale)
 ######  plot the box
PLOTbox(Rax, Rbox, axcol= 'black', boxcol= 'blue')
 ##############  plot the plane
PLOTplane(Rp)

#######  plot the normal
 NN = NORMvec(PPs, xscale, Rview,aglyph=aglyph, add=TRUE)

}
\keyword{misc}
