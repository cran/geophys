\name{setSTRESS}
\alias{setSTRESS}
\title{Set Stress Tensor
}
\description{Set Stress Tensor: given either
  3-vector, or 6-vector, create the symmetric stress tensor matrix.
  
}
\usage{
setSTRESS(Stensor)
}
\arguments{
  \item{Stensor}{
    stress tensor as a vector or a 3 by 3 matrix,
    a 3-vector (principal stresses)  or a 6-vector (principal and off
    diagonal in row-order.
    
}
}
\details{
  Can be either three components of diagonal,
  or six components of the stress tensor.
  If the number of elements is neither 3 nor 6, returns NULL
  and an error message.
}
\value{ List: eigen value decomposition
}

\author{
Jonathan M. Lees<jonathan.lees@unc.edu>
}

\examples{

ES = setSTRESS(c(12, -6, 3) )

ES = setSTRESS(c(12, -6, 3, 3, 6, 1) )




}
\keyword{misc}
