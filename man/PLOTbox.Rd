\name{PLOTbox}
\alias{PLOTbox}
\title{Plot 3D box
}
\description{Plot 3D box
}
\usage{
PLOTbox(Rax, Rbox, axcol = "black", boxcol = "blue")
}
\arguments{
  \item{Rax}{rotated axes
}
  \item{Rbox}{rotated box
}
  \item{axcol}{axes color
}
  \item{boxcol}{box color
}
}

\value{Graphical Side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\seealso{pstart,  PLOTplane, stress, NORMvec
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




}
\keyword{misc}
