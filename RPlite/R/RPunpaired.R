RPunpaired <-
function (lMs, lOs, nperms = 100, size = 64, mc.cores = NA) 
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
