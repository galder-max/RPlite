correctlocalQval <-
function (qval, N = NA) 
{
    if (is.na(N)) 
        N <- length(qval)
    ord <- order(qval, decreasing = F)
    init <- 1
    follow <- 1
    while (follow < N) {
        if (ord[follow] < ord[follow + 1]) {
            follow <- follow + 1
        }
        else {
            qval[init:ord[follow]] <- qval[ord[follow]]
            follow <- which(ord == ord[follow] + 1)
            if (length(follow) == 0) 
                return(qval)
            init <- ord[follow] + 1
        }
    }
    return(qval)
}
