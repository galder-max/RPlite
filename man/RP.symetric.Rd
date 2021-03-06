\name{RP.symetric}
\alias{RP.symetric}

\title{
RP symetric & multicore
}
\description{

}
\usage{
RP.symetric(lMs, lOs, nperms = 100, size = 64, mc.cores = NA)
}

\arguments{
  \item{lMs}{
list of matrices in the meta-analysis. A list containing one matrix can
  be passed. Or a list containing multiple matrices, with same number of
  rows and same ordering of rows (the mapping is not performed, nor
  checked for!)
}
  \item{lOs}{
list of 2-level factors corresponding to the matrices in lMs. Attention
  must be paid to the order of the levels in each factor. Ratios are
  always level1/level2. Ideally, each vector should heave the two same
  unordered levels.
}
  \item{nperms}{
    number of permutations to perform. 100 is typically good. Whereas
    RP follows a f(gamma distribution) for paired expriments, it is not
    the case for unpaired experiments. The null distribution (ND) is
    approximated by permutations. The more permutations, the better the
    approximation, the more time it takes to compute the ND.
}
  \item{size}{
    Very useful. this tells the function to cut the matrix of possible ratios is
    smaller matrices with size columns. It reduces memory consumption by
    very little increase in time consumption. Memory usage can get
    critical, especially for big meta-analyses. 
}
  \item{mc.cores}{
    Number of cores to use if parallelized. Decrease time consumption
    proportionnaly but increases memory needs.
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
nc1<-20
nc2<-10
##genes
nr<-2000
## matrix of expression
m1<-matrix(rnorm(nc1*nr),nr,nc1)
m2<-matrix(rnorm(nc2*nr),nr,nc2)

##factor of paired leves
o1<-as.factor(c(rep("level1",nc1/2),rep("level2",nc1/2)))
o2<-as.factor(c(rep("level1",nc2/2),rep("level2",nc2/2)))

## Rank Product for first matrix
resRP<-RP.symetric(list(m1),list(o1),nperms=100)
head(resRP$lStats)


##Down regulated Genes
sum(resRP$lStats$qvals.do<0.05)
##Up regulated Genes
sum(resRP$lStats$qvals.do<0.05)




## meta-analysis with Rank Product (2 datasets/matrices) with same
## number of rows and ordered rows (mapping is not performed)!
resRPmeta<-RP.symetric(list(m1,m2),list(o1,o2),nperms=100)

head(resRPmeta$finalRes)
##Down regulated Genes
sum(resRPmeta$lfinalRes$qvals.do<0.05)
##Up regulated Genes
sum(resRPmeta$lfinalRes$qvals.up<0.05)



## ADDING signal (mean fold change=2,sd=.1) to 50 first genes
wsignal<-1:50
m1[wsignal,1:(nc1/2)]<-m1[wsignal,1:(nc1/2)]-rnorm(length(m1[wsignal,1:(nc1/2)]),mean=2,sd=.1)
m2[wsignal,1:(nc2/2)]<-m2[wsignal,1:(nc2/2)]-rnorm(length(m2[wsignal,1:(nc2/2)]),mean=2,sd=.1)


resRPmeta<-RP.symetric(list(m1,m2),list(o1,o2),nperms=100)
##Down regulated Genes
sum(resRPmeta$finalRes$qvals.do<0.05)
##Up regulated Genes
sum(resRPmeta$finalRes$qvals.up<0.05)


## plot scores for all genes
cols<-rep("grey",nr)
cols[resRPmeta$finalRes$qvals.up<0.05]<-"red"
barplot(resRPmeta$finalRes$scores,
border=NA,
col=cols,
ylim=c(-500,500),
ylab="RP scores")
## add thresholds of RP scores corresponding to bonferroni-corrected
## p-values 0.05
abline(h=quantile(unlist(resRPmeta$null.distrib),
probs=c(0.05/nr,1-0.05/nr)),
col="grey",lty=2)

}

