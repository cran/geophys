\name{ModelG}
\alias{ModelG}
\title{Interactive 2.5D Gravity Modeling
}
\description{Interactive 2.5D Gravity Modeling
}
\usage{
ModelG(Data, labs = c("Done"), obs=NULL, ZCOLS=RPMG::pastel.colors(24, seed=2 ))
}
\arguments{
  \item{Data}{list: xmin, xmax, zmin, zmax
}
  \item{labs}{Labels
  }
   \item{obs}{Observed Profile Gravity (x,g)
  }
  \item{ZCOLS}{ Color Palette for plotting blocks
  }
   
}
\details{Forward modeling using the
  code of Won and Bevis.

  For input, at a minimum, one must provide the boundary of
  the model region, xmin, xmax, ymin, ymax (or depth).
  more comprehensive input,
  i.e. with a starting model,
  will provide a model in addition:
  xmin, xmax, zmin, zmax, mod,  cens, n.
  The mod is a list of polygons, each
  with xy coordinates and a delta-rho value.

  
}
\value{Model
}
\references{
  Won and Bevis (1987) Computing the gravitational and magnetic anomalies due to a polygon: Algorithms and Fortran subroutines <doi:https://doi.org/10.1190/1.1442298>
}
\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}
\note{
  You type in the polygons
  interactively and execute the gravity forward modeling.
}


\seealso{DGzx
}
\examples{

\dontrun{
blum = list(xmin=0, xmax=26000, zmin=0, zmax=-1000)
gill = ModelG(blum)


}

}
\keyword{misc}
