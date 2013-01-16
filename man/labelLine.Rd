\name{labelLine}
\alias{labelLine}
\title{Label a Line Segment
}
\description{Label a Line Segment
}
\usage{
labelLine(P1, P2, above = TRUE, dinch = 0.2, lab = "text",
 acode = 3, alength = 0.06, aty = 1, acol = "black", bty = 1,
 bcol = "black", tcol = "black", font = 1, cex = 1)
}
\arguments{
  \item{P1}{Point 1 c(x,y)
}
  \item{P2}{Point 2  c(x,y)
}
  \item{above}{logical, TRUE=label above the line, else below
}
  \item{dinch}{length the legs
}
  \item{lab}{character, text label
}
  \item{acode}{code for arrows, see arrows
}
  \item{alength}{length for arrows
}
  \item{aty}{lty for arrows
}
  \item{acol}{color for arrows
}
  \item{bty}{style for legs
}
  \item{bcol}{color for legs
}
  \item{tcol}{color for text
}
  \item{font}{font for text
}
  \item{cex}{character expansion for text, see par
}
}
\details{
  Two short lines are drawn perpendicular to the
  line between the points, the length of this line is dinch.
  The arrow is drawn between the legs,
  using the parameters provided.
  
  
}
\value{graphical side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{par, arrows, text
}
\examples{


 V1 = c(runif(1) ,  runif(1),runif(1) , runif(1))

  
 P1 = c(V1[1], V1[3])
P2 =  c(V1[2], V1[4])
 plot(c(P1[1], P2[1]), c(P1[2], P2[2] ), asp=1, type='n' )
   arrows(P1[1], P1[2], P2[1], P2[2], length=.04, col='red')
    
   labelLine( P1, P2 , lab="ABOVE", dinch = .5,
           aty=2, acol='blue'  )
   labelLine( P1, P2 , above=FALSE, lab="below",
         dinch = .5, aty=2, acol='green', tcol="magenta"  )


 


}
\keyword{misc}
