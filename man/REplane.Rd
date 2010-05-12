\name{REplane}
\alias{REplane}
\title{Replot plane
}
\description{Replot the plane after points have moved
}
\usage{
REplane(m1, Lp, PPs, Rbox, Rview, xscale)
}
\arguments{
  \item{m1}{axis chosen
}
  \item{Lp}{x,y location to move to.
}
  \item{PPs}{Points matrix
}
  \item{Rbox}{Box matrix
}
  \item{Rview}{Viewing matrix
}
  \item{xscale}{arbitrary scale
}
}
\details{Used internally in stress program.
  The new locations are constrained to lie on
  associated axis.

  
}
\value{graphical side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{stress
}

\keyword{misc}
