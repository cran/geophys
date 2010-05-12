\name{rot2Zplane}
\Rdversion{1.1}
\alias{rot2Zplane}
\title{Rotate to Z-plane
}
\description{Rotation matrix to project a set of 3D points to z-plane
}
\usage{
rot2Zplane(vec, p)
}
\arguments{
  \item{vec}{3D normal vector 
}
  \item{p}{translation from this point
}
}
\details{provides a matrix to rotate and translate a
  set of points to the X-Y plane.  Used for internal calculations
}
\value{Matrix for projection
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}


\seealso{gmat
}
\examples{


P1=runif(3)
P2=runif(3)
P3=runif(3)

PP = rbind(P1, P2, P3, rep(1, 3) )

    g1 = PP[1, 1:3] - PP[3, 1:3]
    g2 = PP[2, 1:3] - PP[3, 1:3]


 B = unlist(AXB.prod(list(x = g1[1], y = g1[2], z = g1[3]), 
        list(x = g2[1], y = g2[2], z = g2[3])))
    B = B/sqrt(sum(B * B))


MF = rot2Zplane(B, P3) 

UMAT = t(PP) %*% MF


}
\keyword{misc}
