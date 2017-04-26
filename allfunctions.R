
setwd("/Users/Maxime/Desktop/RPlite/")
library(RPlite) ## old one

RPpaired<-RPlite:::RPpaired

rankprodbounds <- function(rho,n,k,Delta){ 
  ## I am not the author of this function
  ##
  ## This was adapted from http://www.ru.nl/ publish/pages/726696/rankprodbounds.zip
  ##
  ## Described in:
  ##
  ## "A fast algorithm for determining bounds and accurate approximate p-values of the rank
  ## product statistic for replicate experiments"
  ## Tom Heskes, Rob Eisinga and Rainer Breitling
  ##
  ## As stated in the original code: 
  ## <<This implementation closely follows the description in Heskes, Eisinga, Breitling:
  ## "A fast algorithm for determining bounds and accurate approximate p-values of the
  ## rank product statistic for replicate experiments", further referred to as HEB.
  ## More specifically, this R function corresponds to the recursive variant, sketched
  ## as pseudocode in the additional material of HEB.
  ##
  ## updated version August 2015: fixed a bug with the help of Vicenzo Lagani>> 
  updateparam <- function(param,n,k,j,Delta) {
    k1 <- 1+k
    j1 <- 1+j 
    if(length(param[[k1,j1]]$e) == 0) {  # apparently empty, so needs to be calculated   
      if(j == 0) {   # initializing G_{k0}      
        param[[k1,j1]]$e <- n^k
        param[[k1,j1]]$d <- 0
                                        # the 2 lines above implement equation (11) in HEB
      }
      else {
        k0 <- k1-1
        j0 <- j1-1
        param <- updateparam(param,n,k-1,j-1,Delta)
                                        # checking that the parameters for (k-1,j-1) that are needed to compute the
            # parameters for (k,j) are indeed available; if not, they are themselves computed
        param00 = param[[k0,j0]]
        newa0 = param00$a+1
        newb0 = param00$b
        newc0 = param00$c/newa0
        param11 = param00
                                        # the 5 lines above predefine some parameters common to equations (9) and (10) in HEB    
        if(k == j){ # updates for G_{kk}    
          param11$e <- (1-Delta)*(1-param00$e)
          param11$d <- Delta*param00$d+param00$e
          param11$a <- c(1,param00$a,newa0)
          param11$b <- c(0,param00$b,newb0)
          param11$c <- c(param00$d,Delta*param00$c,newc0)
                                        # the 5 lines above implement equation (10) in HEB
        }
        else {  # updates for G_{kj}, j < k
          param <- updateparam(param,n,k-1,j,Delta)
                                        # checking that the parameters for (k-1,j) that are needed to compute the
                                        # parameters for (k,j) are indeed available; if not, they are themselves computed
          param01 <- param[[k0,j1]]      
          logn <- log(n)
          lognnkj <- (k-j)*logn
          newa1 <- param01$a+1
          newa <- c(newa0,newa1)
          newb <- c(newb0,param01$b)
          newc <- c(newc0,-param01$c/newa1)
          param11$e <- n*param01$e + (Delta-1)*(param00$e-param01$e)
          lognminb <- c(-1*param00$b * logn,(1-param01$b)*logn)
          param11$d <- Delta*param00$d + (1-Delta)*param01$d/n + 
            (param00$e-param01$e)/exp(lognnkj) - 
              sum(newc*(lognminb^newa))
          param11$a <- c(1,1,param00$a,param01$a,newa)
          param11$b <- c(0,1,param00$b,param01$b,newb)
          param11$c <- c(param00$d,-param01$d,
                         Delta*param00$c,(1-Delta)*param01$c/n,newc)
                                        # the 15 lines above implement equation (9) in HEB
        }
        param[[k1,j1]] <- makeunique(param11)
                                        # although not strictly necessary, the a, b and c vectors can possibly be shortened by
                                        # restricting oneselves to unique combinations of a and b values
      }
    }
    return(param)
  }
  
  makeunique <- function(param) { 
    ab <- t(rbind(param$a,param$b))
    uniqueab <- unique(ab)
    nunique <- dim(uniqueab)[1]
    param$a <- t(uniqueab[,1])
    param$b <- t(uniqueab[,2])
    newc <- rep(0,nunique)
    for(i in 1:nunique) {
      iii <- intersect(which(ab[,1]==uniqueab[i,1]),which(ab[,2]==uniqueab[i,2]))
      newc[i] <- sum(param$c[iii])  
    }
    param$c <- newc  
    return(param)  
  }  
  if(any(rho > n^k) || any(rho < 1)) stop('rho out of bounds') 
  if(is.numeric(Delta) == FALSE) {
    if(Delta == 'geometric') {
      temp1 <- rankprodbounds(rho,n,k,'upper')
      temp2 <- rankprodbounds(rho,n,k,'lower')
      pvalue <- sqrt(temp1*temp2)   # geometric mean of upper and lower bound
      return(pvalue)
    }
    else {
      Delta <- switch(Delta,
        upper = 1,        # for computing upper bound
        lower = 0)        # for computing lower bound
    }
  }
  logn <- log(n)
  allj <- ceiling(-(log(rho)/logn)+k)   # index specifying the interval that contains rho 
  minj <- min(allj)                     # lowest interval index
  maxj <- max(allj)                     # highest interval index
  param <- matrix(list(), nrow=k+1, ncol=maxj+1)
  for(i in 1:(k+1)){
	  for(j in 1:(maxj+1)){
		  param[[i,j]] <- list(a=c(),b=c(),c=c(),d=c(),e=c())
	  }
  } 
  for(j in minj:maxj){
    param <- updateparam(param,n,k,j,Delta)
  }  
  k1 <- 1+k
  G <- rep(0,length(rho))   # G is a vector of the same length as rho,
                            # for each rho bounding the number of rank products 
  for(j in unique(allj)) {  # updated: thanks to Vicenzo Lagani for pointing this out
    j1 <- 1+j
    iii <- which(allj == j)         # indices of all rank products that fall in interval j:
                                    # bounds for these rank products can be computed with
                                    # the same set of parameters                                    
    thisrho <- rho[iii]
    thisparam <- param[[k1,j1]]
    thisG <- thisparam$e
    if(j != 0) {
      nrho <- length(thisrho)
      nterms <- length(thisparam$a)
      thisG <- thisG + thisparam$d*thisrho
      d1 <- matrix(thisparam$c) %*% thisrho
      d2 <- matrix(rep(log(thisrho),nterms),nrow=nterms,byrow=TRUE) -
        t(matrix(rep(logn*(k-j+thisparam$b),nrho),nrow=nrho,byrow=TRUE))
      d3 <- t(matrix(rep(thisparam$a,nrho),nrow=nrho,byrow=TRUE)) 
      thisG <- thisG + colSums(d1*(d2^d3))
    }
    # the 10 lines above implement equation (8) in HEB
    G[iii] <- thisG
  }
  pvalue <- G/n^k
  return(pvalue)
}

RPpaired<-function (m, nperms = 200, method=c("gamma","permutation","geometricLU")) {
    Ngenes <- nrow(m)
    Nsamplepairs <- ncol(m)
    ranks <- apply(m, 2, function(x) rank(x))
    if(method[1]=="permutation")
      {
#        ranks<--log(ranks/(1+Ngenes-ranks))
        ranks<--log(ranks)
        scores<-rowSums(ranks)
        nulld <- unlist(lapply(1:nperms,function(x)
                               {
                                 rowSums(apply(ranks,2,sample))
                               }))
        fr <- RPlite:::compute.Stats(nulld, scores, Ngenes, nperms)
      }
    if(method[1]=="gamma")
      {
        scores <- -rowSums(log(ranks)) + ncol(m) * log(Ngenes + 1)
        nulld <- rgamma(Ngenes * nperms, shape = Nsamplepairs, scale = 1)
        fr <- RPlite:::compute.Stats(nulld, scores, Ngenes, nperms)
      }
    if(method[1]=="geometricLU")
      {
        scores <- apply(ranks,1,prod)
        nulld<-NULL
        fr<-list()
        err<-try(fr$pvals.do<-rankprodbounds(scores,n=Ngenes,k=ncol(m),Delta="geometric"), silent=F)
        if(sum(grepl("Error",unlist(err)))>0) print("Please, consider using method permutation or gamma instead of geomeitrcLU") 
        fr$pvals.up<-1-fr$pvals.do
        fr$qvals.do<-p.adjust(fr$pvals.do,method="fdr")
        fr$qvals.up<-p.adjust(fr$pvals.up,method="fdr")
        fr$scores <- scores
      }
    return(list(lStats = fr, null.distrib = list(nulld), final.scores = list(scores), 
                nperms = nperms))
}

RPunpaired<-RPlite:::RP.symetric

RPunpaired<-function (lMs, lOs, nperms = 100, size = 64, mc.cores = NA) 
{
    require(parallel)
    if (is.na(mc.cores)) 
        mc.cores <- parallel:::detectCores(all.tests = TRUE) - 
            2
    nb.genes <- nrow(lMs[[1]])
    level1 <- levels(as.factor(lOs[[1]]))[1]
    level2 <- levels(as.factor(lOs[[1]]))[2]
    print(paste(level1, "/", level2))
    final.scores <- lapply(1:length(lOs), function(x) giveRPsymmetric.divideMatrix(list(lMs[[x]]), 
                                                                                   list(lOs[[x]]),
                                                                                   size = size,
                                                                                   doPrint = (x == 1),
                                                                                   mc.cores = mc.cores))
    gc()
    print("Null distribution of scores")
    null.distrib <- lapply(1:length(lMs), function(y) {
        if (y != 1) 
            cat("\n")
        print(paste("   Dataset", y))
        pb <- txtProgressBar(style = 3)
        lapply(1:nperms, function(perm) {
            setTxtProgressBar(pb, perm/nperms)
            yt <- apply(lMs[[y]], 2, sample)
            giveRPsymmetric.divideMatrix(list(yt), list(lOs[[y]]), 
                                         size = size,
                                         doPrint = FALSE,
                                         mc.cores = mc.cores)
        })
      })
    cat("\n")
    gc()
    print("Ordering and computing pvals")
    print("    Computing separate Stats")
    lStats <- lapply(1:length(lMs), function(x) compute.Stats(unlist(null.distrib[[x]]), 
        final.scores[[x]], nb.genes, nperms))
    if (length(lMs) > 1) {
        print("    Computing Null distribution")
        null.distribF <- lapply(1:nperms, function(y) {
            rowSums(sapply(lapply(1:length(lMs), function(x) null.distrib[[x]][[y]]), 
                unlist))
        })
        print("    Computing final scores")
        final.scoresF <- rowSums(as.data.frame(final.scores))
        print("    Computing final result")
        finalRes <- compute.Stats(unlist(null.distribF), final.scoresF, 
            nb.genes, nperms)
        return(list(lStats = lStats, final.scores = final.scores, 
            final.scoresF = final.scoresF, finalRes = finalRes, 
            null.distrib = null.distrib, null.distribF = null.distribF))
    }
    else {
        return(list(lStats = lStats[[1]], final.scores = final.scores, 
            final.scoresF = NULL, finalRes = NULL, null.distrib = null.distrib, 
            null.distribF = NULL))
    }
}

giveRPsymmetric<-RPlite:::giveRPsymetric
giveRPsymmetric.divideMatrix<-RPlite:::giveRPsymetric.divideMatrix

giveRPsymmetric.divideMatrix<-function (lMs, lOs, doPrint = FALSE, size = 64, mc.cores = 1) 
{
    require(parallel)
    pairs <- find.pairs(lOs)
    l <- length(pairs[[1]])
    nbs <- l%/%size
    rest <- l%%size
    if (nbs == 0) 
        return(giveRPsymmetric(lMs, lOs, pairs, doPrint))
    if (nbs > 1) {
        res <- sapply(mclapply(1:(nbs - 1), function(i) {
            seq. <- ((i - 1) * size + 1):(i * size)
            pairs. <- lapply(seq., function(x) pairs[[1]][[x]])
            return(giveRPsymmetric(lMs, lOs, list(pairs.), doPrint = FALSE))
        }, mc.cores = mc.cores), unlist)
        res <- rowSums(res)
    }
    else res <- rep(0, nrow(lMs[[1]]))
    if (rest >= size%/%2) {
        seq. <- ((nbs - 1) * size + 1):(nbs * size)
        pairs. <- lapply(seq., function(x) pairs[[1]][[x]])
        res <- res + giveRPsymmetric(lMs, lOs, list(pairs.), doPrint = FALSE)
        seq. <- (nbs * size + 1):(nbs * size + rest)
        pairs. <- lapply(seq., function(x) pairs[[1]][[x]])
        res <- res + giveRPsymmetric(lMs, lOs, list(pairs.), doPrint = doPrint)
    }
    else {
        seq. <- ((nbs - 1) * size + 1):(nbs * size + rest)
        pairs. <- lapply(seq., function(x) pairs[[1]][[x]])
        res <- res + giveRPsymmetric(lMs, lOs, list(pairs.), doPrint = doPrint)
    }
    return(res)
}

find.pairs<-RPlite:::find.pairs
construct.df.scores<-RPlite:::construct.df.scores
correctlocalQval<-RPlite:::correctlocalQval
compute.Stats<-RPlite:::compute.Stats
integrate.Results.General<-RPlite:::integrate.Results.General


package.skeleton(name="RPlite",list=ls())
