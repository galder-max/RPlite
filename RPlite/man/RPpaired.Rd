\name{RPpaired}
\alias{RPpaired}

\title{
RP paired 
}
\description{

}
\usage{
RPpaired(m, nperms = 200, method = c("gamma", "permutation", "geometricLU"))
}

\arguments{
  \item{m}{
    matrix of paired ratios.
}
  \item{nperms}{
    number of permutations to perform to approximate the f(gamma distribution).
}
  \item{method}{
    defines the way p-values are computed: either with permutation,
    gamma approximation or using the method described by Heskes et
    al. (see description of the rankprodbounds function)
}
}
\details{
  
}
\value{
  returns a list including p-values and q-values for each
  gene/probesets/isoform 
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
nc<-10
##genes
nr<-10000
## matrix of expression
m<-matrix(rnorm(nc*nr),nr,nc)
##factor of paired leves
o<-c(paste("level1",1:(nc/2),sep="-"),paste("level21",1:(nc/2),sep="-"))

##matrix of ratios
mr<-m[,1:(nc/2)]-m[,(nc/2+1):nc]


## paired Rank Product
resRP1<-RPpaired(mr,nperms=100)
head(resRP1$lStats)

resRP2<-RPpaired(mr,method="gamma")
head(resRP2$lStats)

resRP3<-RPpaired(mr,method="geometricLU")
str(resRP3$lStats)

plot(resRP1$lStats$pvals.do,
resRP2$lStats$pvals.do,pch=".")

plot(resRP1$lStats$pvals.do,
resRP3$lStats$pvals.do,pch=".")

plot(resRP3$lStats$pvals.do,
resRP2$lStats$pvals.do,pch=".")


}
\keyword{ ~kwd1 }
