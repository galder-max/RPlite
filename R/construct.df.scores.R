construct.df.scores <-
function(lMs,lOs,lIsPaireds,pairs)
  {    
    number.pairs<-sapply(pairs, length)
    df.scores<-matrix(NA,ncol=sum(number.pairs),nrow=max(sapply(lMs,nrow)))
    count<-1
    for(i.datasets in 1:length(lOs))
      {
        .pairs<-pairs[[i.datasets]]
        for (j.dataset in 1:number.pairs[i.datasets])
          {
            df.scores[,count]<-lMs[[i.datasets]][,.pairs[[j.dataset]][1]]-lMs[[i.datasets]][,.pairs[[j.dataset]][2]]
            count<-count+1
          }
      }
    return(df.scores)    
  }
