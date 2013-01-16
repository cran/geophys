\name{PLOTplane}
\alias{PLOTplane}
\title{Plot Plane on box
}
\description{Plot Plane on box
}
\usage{
PLOTplane(Rp, planecol = "brown")
}
\arguments{
  \item{Rp}{Rotated points
}
  \item{planecol}{color for plotting plane
}
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{pstart, PLOTbox, stress
}
\examples{

#######    set colors
 axcol = 'black'
 boxcol = 'blue'
 planecol = 'brown'

####   view angle
 Rview  =    RFOC::ROTZ(-130) %*% RFOC::ROTX(-50) 

###   set arbitrary scale
    xscale = 100


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

     P1 = xscale*c(0, 0,.4,0)
    P2 = xscale*c(0.8, 0,0,0)
    P3 = xscale*c(0, .7, 0 ,0)

    PPs = rbind(P1, P2, P3)
    
  ###  convert points for view
Rp = PPs  %*% Rview

PLOTplane(Rp, planecol = "gold" )


}
\keyword{misc}
