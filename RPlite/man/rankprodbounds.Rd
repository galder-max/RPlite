\name{rankprodbounds}
\alias{rankprodbounds}

\title{
  rankprodbounds takes RP scores of paired comparisons as input and
  quickly outputs a vector of good approximate corresponding p-values. 
}
\description{
  This was adapted from http://www.ru.nl/ publish/pages/726696/rankprodbounds.zip
  
  Described in:
  
  A fast algorithm for determining bounds and accurate approximate p-values of the rank
  product statistic for replicate experiments
  Tom Heskes, Rob Eisinga and Rainer Breitling
}
\usage{
rankprodbounds(rho, n, k, Delta)
}

\arguments{
  \item{rho}{
    vector of RP scores
}
  \item{n}{
    number of samples
}
  \item{k}{
    number of molecules, genes, isoforms or probesets
}
  \item{Delta}{
    the type of output, only geometric is used inside RPlite and it
    outputs the geometric mean between the p-value's upper and lower bounds
}
}
\details{

}
\value{
  the vector of p-values
}
\references{
  A fast algorithm for determining bounds and accurate approximate p-values of the rank
  product statistic for replicate experiments
  Tom Heskes, Rob Eisinga and Rainer Breitling
}
\author{
  I am not the author, please refer to the publication :
  A fast algorithm for determining bounds and accurate approximate p-values of the rank
  product statistic for replicate experiments
  Tom Heskes, Rob Eisinga and Rainer Breitling
}
\note{

}


\seealso{

}
\examples{

}

\keyword{ ~kwd1 }

