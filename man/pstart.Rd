\name{pstart}
\alias{pstart}
\title{Start an empty plot
}
\description{Start an empty plot
}
\usage{
pstart(xscale = 30, expand=1.2)
}
\arguments{
  \item{xscale}{Arbitray scale for plot
  }
  \item{expand}{factor to expand the plotting region, default=1.2
  }


  
}
\details{Creates a plotting region with no drawing.
  Default region is defined by:
 xlim= c(-expand*xscale, expand*xscale) and
 ylim= c(-expand*xscale, expand*xscale).

  


  
}
\value{Graphical Side effects
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{

pstart()

axis(1)
axis(2)


}
\keyword{misc}
