RPpaired<-function(m,nperms=200)
  {
    Ngenes<-nrow(m)
    Nsamplepairs<-ncol(m)
    ranks<-apply(m,2,function(x) rank(x))
    scores<--rowSums(log(ranks))+ncol(m)*log(Ngenes+1)
    nulld<-rgamma(Ngenes*nperms,shape=Nsamplepairs,scale=1)
    fr<-compute.Stats(nulld,
                      scores,
                      Ngenes,
                      nperms)
    return(list(lStats=fr,
                null.distrib=list(nulld),
                final.scores=list(scores),
                nperms=nperms))
  }


