find.pairs<-function(lOs)
  {
    if(any(sapply(lOs,function(x) length(levels(x)))!=2))
      stop("Class factor vectors must be of length 2")
    pairs<-lapply(1:length(lOs),function(x)
                  {
                    l1<-which(lOs[[x]]==levels(as.factor(lOs[[x]]))[1])
                    l2<-which(lOs[[x]]==levels(as.factor(lOs[[x]]))[2])
                    vec<-list()
                    count<-1
                    for (i in l1)
                      for (j in l2)
                        {
                          vec[[count]]<-c(i,j)
                          count<-count+1
                        }
                    return(vec)
                  })
    return(pairs)
  }
