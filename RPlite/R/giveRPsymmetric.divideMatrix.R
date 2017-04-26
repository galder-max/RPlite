giveRPsymmetric.divideMatrix <-
function (lMs, lOs, doPrint = FALSE, size = 64, mc.cores = 1) 
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
