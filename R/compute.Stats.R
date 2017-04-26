compute.Stats <-
function(null.distrib,final.scores, nb.genes, nperms)
  {
    ord.finalscore<-order(final.scores,decreasing=FALSE)
    final.scores.ordered<-final.scores[ord.finalscore]
    ords<-order(c(final.scores.ordered,null.distrib),decreasing=FALSE)
    nb.better<-(match(1:nb.genes,ords)-1:nb.genes)
    pvals.up<-nb.better/nperms/nb.genes
    pvals.do<-1-pvals.up
    qvals.up<-(nb.better/1:nb.genes)/nperms
    qvals.up[qvals.up>1]<-1
    qvals.do<-((nperms*nb.genes)-nb.better)/nb.genes:1
    qvals.do[qvals.do>1]<-1
    revert.ind<-match(1:nb.genes,ord.finalscore)    
    qvals.do[nb.genes:1]<-correctlocalQval(qvals.do[nb.genes:1])
    qvals.up<-correctlocalQval(qvals.up)
    results<-data.frame(rows=1:nb.genes,
                        scores=final.scores,
                        rank.scores=revert.ind,
                        pvals.up=pvals.up[revert.ind],
                        qvals.up=qvals.up[revert.ind],
                        pvals.do=pvals.do[revert.ind],
                        qvals.do=qvals.do[revert.ind])
  }
