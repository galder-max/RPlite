integrate.Results.General <-function(listRes,vec=NA, weights=NA,nperms=NA)
  {
    if(is.na(vec[1])) vec<-1:length(listRes)
    if(is.na(weights[1])) weights<-rep(1/length(vec),length(listRes))
    nb.genes<-length(listRes[[1]]$final.scores[[1]])
    cnperms<-sapply(vec,function(x)
                    if(!is.null(listRes[[x]]$null.distribF))
                    length(listRes[[x]]$null.distribF)
                    else length(listRes[[x]]$null.distrib[[1]]))
    print(cnperms)
    if(!is.na(nperms))
      {
        nperms<-min(cnperms)[1]
        null.distribF<-lapply(1:nperms,function(y)
                              {
                                rowSums(sapply(lapply(vec,
                                                      function(x)
                                                      if(!is.null(listRes[[x]]$null.distribF))
                                                      (listRes[[x]]$null.distribF[[y]])*(weights[x])
                                                      else
                                                      (listRes[[x]]$null.distrib[[1]][[y]])*(weights[x])
                                                      ),unlist))
                              })
      }
    else
      {
        null.distribF<-lapply(1:nperms,function(y)
                              {
                                rowSums(sapply(lapply(vec,
                                                      function(x)
                                                      if(!is.null(listRes[[x]]$null.distribF))                                                      
                                                      (listRes[[x]]$null.distribF[[sample(1:cnperms[x],1)]])*(weights[x])
                                                      else
                                                      (listRes[[x]]$null.distrib[[1]][[sample(1:cnperms[x],1)]])*(weights[x])
                                                      ),unlist))
                              })        
      }
    final.scoresF<-rowSums(sapply(lapply(vec,function(x)
                                         if(!is.null(listRes[[x]]$final.scoresF))
                                         (listRes[[x]]$final.scoresF)*(weights[x])
                                         else
                                         (listRes[[x]]$final.scores[[1]])*(weights[x])   
                                         ),unlist))  
    finalRes<-compute.Stats(unlist(null.distribF),final.scoresF,nb.genes, nperms)
    return(list(null.distribF=null.distribF,
                finalRes=finalRes,
                final.scoresF=final.scoresF))
  }
