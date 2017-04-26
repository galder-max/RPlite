giveRPsymmetric <-
function (lMs, lOs, pairs, doPrint = FALSE) 
{
    nb.genes <- nrow(lMs[[1]])
    if (doPrint) 
        print("Constructing Ratios")
    df.scores <- construct.df.scores(lMs, lOs, lIsPaireds, pairs)
    ncol.df.scores <- ncol(df.scores)
    if (doPrint) 
        print("Ranking")
    ranks <- apply(df.scores, 2, rank)
    if (doPrint) 
        print("Final Scores")
    final.scores <- rowSums(log(ranks/((nb.genes + 1) - ranks)))
    return(final.scores)
}
