\name{Glines}
\Rdversion{1.1}
\alias{Glines}
\docType{data}
\title{Igneous Petrology Lines

}
\description{Scanned lines from igneous petrology

}
\usage{data(Glines)}
\format{
  The format is:
List of 10
 $ :List of 2
  ..$ x: num [1:77] 55.2 56.1 57.2 58.5 59.5 ...
  ..$ y: num [1:77] 1.62 1.67 1.82 2.15 2.35 ...
 $ :List of 2
  ..$ x: num [1:6] 50.7 51.7 53.1 54.4 56.3 ...
  ..$ y: num [1:6] 9.21 9.16 9.24 9.21 9.09 ...
 $ :List of 2
  ..$ x: num [1:14] 41.7 41.9 42.8 44.9 48.4 ...
  ..$ y: num [1:14] 3.35 3.55 4.22 5.82 8.49 ...
 $ :List of 2
  ..$ x: num [1:13] 52.1 52.2 52.2 52.2 52.3 ...
  ..$ y: num [1:13] 1.72 2.7 3.77 5.27 5.67 ...
 $ :List of 2
  ..$ x: num [1:20] 45.5 46.8 47.8 48.8 49.7 ...
  ..$ y: num [1:20] 5.74 5.69 5.69 5.69 5.64 ...
 $ :List of 2
  ..$ x: num [1:2] 63.2 63.4
  ..$ y: num [1:2] 6.64 3.45
 $ :List of 2
  ..$ x: num [1:2] 55.2 55.2
  ..$ y: num [1:2] 5.34 1.75
 $ :List of 2
  ..$ x: num [1:2] 46.7 53.2
  ..$ y: num [1:2] 7.04 7.07
 $ :List of 2
  ..$ x: num [1:2] 49.7 55.1
  ..$ y: num [1:2] 14.7 11.2
 $ :List of 2
  ..$ x: num [1:5] 58.8 60 61.8 62.3 62.5
  ..$ y: num [1:5] 11.5 12.3 13.4 13.8 13.9
}
\details{
  Can click on here with a mixture to get
  a location on the plot.

}
\source{scanned figure

}
\references{Principles of Igneous and Metamorphic Petrology, A. R. Philpotts, 1990, Prentice Hall (ISBN-13:978-0521880060)

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
