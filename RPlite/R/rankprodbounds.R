rankprodbounds <-
function(rho,n,k,Delta){ 
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
