\name{PPoints}
\Rdversion{1.1}
\alias{PPoints}
\docType{data}
\title{Igneous Petrology Points

}
\description{Scanned points from igneous petrology

}
\usage{data(PPoints)}
\format{
  The format is:
List of 4
 $ x   : num [1:16] 45.4 56.4 64.1 71.1 66.3 ...
 $ y   : num [1:16] 9.39 13.56 11.36 8.36 5.84 ...
 $ lab : chr [1:16] "minor varieties" "phonolites" "trachytes" "rhyolites" ...
 $ rots: num [1:16] 45 0 0 0 0 0 0 -90 0 0 ...
}
\details{
Can click on here with a mixture to get
  a location on the plot.

}
\source{scanned figure

}
\references{Principles of Igneous and Metamorphic Petrology, A. R. Philpotts, 1990, Prentice Hall 

}
\examples{

data(Glines)
data(PPoints)

plot(Glines[[1]], type='n', xlab="SiO2", ylab="NA2O+K2O")
for(i in 1:length(Glines))
{

lines(Glines[[i]]$x, Glines[[i]]$y)

}

for(i in 1:length(PPoints$x))
{
text(PPoints$x[i],PPoints$y[i], label=PPoints$lab[i], srt=PPoints$rots[i])
}

}
\keyword{datasets}
