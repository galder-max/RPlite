\name{RPpaired}
\alias{RPpaired}

\title{
RP paired 
}
\description{

}
\usage{
RPpaired(m, nperms = 200)
}

\arguments{
  \item{m}{
    matrix of paired ratios.
}
  \item{nperms}{
    number of permutations to perform to approximate the f(gamma distribution).
}
}
\details{

}
\value{

}
\references{

}
\author{
maxime.tarabichi@ulb.ac.be
}
\note{

}

\seealso{

  
}
\examples{
## samples
nc<-20
##genes
nr<-10000
## matrix of expression
m<-matrix(rnorm(nc*nr),nr,nc)
##factor of paired leves
o<-c(paste("level1",1:(nc/2),sep="-"),paste("level21",1:(nc/2),sep="-"))

##matrix of ratios
mr<-m[,1:(nc/2)]-m[,(nc/2+1):nc]


## paired Rank Product
resRP<-RPpaired(mr,nperms=100)
head(resRP$lStats)

## to be done
## use gamma distribution instead of permutations
## for the null distribution

}
\keyword{ ~kwd1 }
